// lz_tail_test.cpp  --  ER recombination tail with LZ's own EFT2024 config
//
// Uses LZ_Detector_2024::SetEFT2024Configuration() so that every yield, width,
// density and W-value convention comes from LZ's committed configuration
// rather than from hand-transcribed table values.
//
// Usage:  ./lz_tail_test [N_trials] [skewness]
//         N_trials = 0 (default) runs PART 1 only
//         skewness = -999 (default) uses the built-in LUX Run04 model

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include "NEST.hh"
#include "VDetector.hh"
#include "LZ_WS2024.hh"

using namespace std;
using namespace NEST;

// LZ Table S3 as printed, for cross-check only.
// NOTE: m6 is a SENTINEL in NEST (only its sign is read), so it is held at 0
// here; this function therefore omits the density correction NEST applies
// when ERYieldsParam[5] < 0.
static double Qy_paper(double E) {
  const double m1=12.31,m2=84.91,m3=0.5707,m4=2.804,m5=34.04,
               m7=82.57,m8=4.557,m9=0.2207,m10=0.1272;
  return m1 + (m2-m1)/pow(1.+pow(E/m3,m4),m9)
       + m5 + (0.0-m5)/pow(1.+pow(E/m7,m8),m10);
}

// Eq. 4 of the Supplement, both ways. printed=true evaluates the expression as
// printed with Table S4; printed=false applies the component-wise mode
// normalisation that RecombOmegaER carries. Both are computed here rather
// than printed as literals, so that the line reports what the build actually
// does rather than what it is expected to do.
static double omega_eq4(double Nq, bool printed) {
  static const double A[2]  = {0.03174, 0.06211};
  static const double MU[2] = {2.588,   5.200  };
  static const double SG[2] = {0.3627,  1.359  };
  static const double AL[2] = {-0.3583, -2.599 };
  const double x = log10(Nq), s2 = sqrt(2.), inv = 1./sqrt(2.*M_PI);
  double t = 0.;
  for (int i = 0; i < 2; ++i) {
    double nrm = 1.;
    if (!printed) {
      const double mo = MU[i] + 2.*inv*AL[i]*SG[i]/sqrt(1.+AL[i]*AL[i]);
      nrm = 1./(exp(-0.5*(mo-MU[i])*(mo-MU[i])/(SG[i]*SG[i]))
                *(1.+erf(AL[i]*(mo-MU[i])/(SG[i]*s2))));
    }
    t += nrm*A[i]*exp(-0.5*(x-MU[i])*(x-MU[i])/(SG[i]*SG[i]))
            *(1.+erf(AL[i]*(x-MU[i])/(SG[i]*s2)));
  }
  return t;
}

int main(int argc, char** argv) {
  const long   N       = (argc>1) ? atol(argv[1]) : 0L;
  const double skew_in = (argc>2) ? atof(argv[2]) : -999.;

  const double E_keV = 70.95;   // combined energy of the event of interest
  const double field = 96.5;    // V/cm, from LZ_Detector_2024::FitEF

  LZ_Detector_2024 det;
  det.SetEFT2024Configuration();
  NESTcalc n(&det);

  bool inGas = false;
  double rho = NESTcalc::GetDensity(175.9, 1.86, inGas);   // LZ T, P

  YieldResult y = n.GetYields(NEST::beta, E_keV, rho, field, 131.293, 54.,
                              det.get_henr_nr_yield_params(),
                              det.get_henr_er_yield_params());

  vector<double> P = det.get_henr_nr_er_width_params();

  double Ne  = y.ElectronYield;
  double Nph = y.PhotonYield;
  double Nq  = Nph + Ne;
  double Ni  = Nq/(1.+y.ExcitonRatio);

  printf("density  = %.4f g/cm^3   inGas = %d   field = %.1f V/cm\n", rho, (int)inGas, field);
  printf("E        = %.2f keVee\n", E_keV);
  printf("NEST  :  Nph = %.1f   Ne = %.1f   Nq = %.1f   Nex/Ni = %.4f\n",
         Nph, Ne, Nq, y.ExcitonRatio);
  printf("         implied W = %.3f eV\n", 1000.*E_keV/Nq);
  printf("paper :  Ne = %.1f  (Table S3 as printed, no density correction)\n",
         Qy_paper(E_keV)*E_keV);
  printf("         ratio NEST/paper = %.4f\n\n", Ne/(Qy_paper(E_keV)*E_keV));

  // ---------------- PART 1: extract omega ----------------
  QuantaResult q0 = n.GetQuanta(y, rho, P, skew_in);
  double r     = q0.recombProb;
  double binom = r*(1.-r)*Ni;
  double om    = sqrt(max(0.,(q0.Variance-binom)))/Ni;

  printf("--- PART 1 (skewness arg = %.1f) ---\n", skew_in);
  printf("  recombProb       = %.5f\n", r);
  printf("  Variance         = %.1f   (binomial part %.1f)\n", q0.Variance, binom);
  printf("  sd(Ne)           = %.1f\n", sqrt(q0.Variance));
  printf("  omega EXTRACTED  = %.5f\n", om);
  const double om_impl  = omega_eq4(Nq, false);
  const double om_print = omega_eq4(Nq, true);
  printf("        code form  = %.5f   Eq.4 as printed = %.5f\n", om_impl, om_print);
  printf("        extracted/code = %.4f   printed/code = %.4f\n",
         om/om_impl, om_print/om_impl);
  if (fabs(om/om_impl - 1.) > 0.01)
    printf("  ** extracted and code-form omega differ by more than 1%%:"
           " check the configuration\n");
  printf("\n");

  if (N <= 0) { printf("(pass N to run PART 2)\n"); return 0; }

  // ---------------- PART 2: tail probability ----------------
  const double NE_CEIL  = pow(10.,4.15)/34.5;   // 409.4312 e-; on an integer
                                                // count this selects N_e <= 409
  const double NE_EVENT = 9268./34.5;           // 269 e-, the event
  long nc=0, nv=0; double s1=0,s2=0,s3=0;
  for (long i=0;i<N;i++){
    QuantaResult q = n.GetQuanta(y, rho, P, skew_in);
    double e=(double)q.electrons; s1+=e; s2+=e*e; s3+=e*e*e;
    if (e<NE_CEIL ) nc++;
    if (e<NE_EVENT) nv++;
  }
  double m=s1/N, v=s2/N-m*m, sd=sqrt(v);
  printf("--- PART 2: %ld trials ---\n", N);
  printf("  mean(Ne) = %8.1f\n", m);
  printf("  sd(Ne)   = %8.1f\n", sd);
  printf("  skew(Ne) = %8.3f\n", (s3/N-3*m*v-m*m*m)/(sd*sd*sd));
  printf("  P(Ne<%3.0f) = %.3e +/- %.1e\n",
         NE_CEIL,  (double)nc/N, sqrt((double)max(1L,nc))/N);
  printf("  P(Ne<%3.0f) = %.3e +/- %.1e\n",
         NE_EVENT, (double)nv/N, sqrt((double)max(1L,nv))/N);
  return 0;
}
