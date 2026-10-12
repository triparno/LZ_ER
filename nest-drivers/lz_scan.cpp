// lz_scan.cpp -- one energy point per process, for parallel running.
//
//   ./lz_scan <E_keV> <N_trials> <seed>
//
// Prints one machine-readable line:
//   E  Nq  NexNi  meanNe  sdNe  skewNe  N  n_ceil  P_ceil  dP_ceil  n_win  P_win  dP_win
//
// Threshold  Ne <= 409 e-       (WS ROI ceiling S2c < 10^4.15 phd, g2 = 34.5)
// Window     4529 < Nph < 5435  (S1c 500-600 phd, g1 = 0.1104 from LZ config)
//
// IMPORTANT: NEST's RandomGen is a singleton seeded deterministically. Without
// SetSeed every process produces identical output. Pass a distinct seed.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cmath>
#include "NEST.hh"
#include "RandomGen.hh"
#include "VDetector.hh"
#include "LZ_WS2024.hh"     // declares class LZ_Detector_2024

using namespace std;
using namespace NEST;

int main(int argc, char** argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s <E_keV> <N_trials> <seed>\n", argv[0]);
    return 1;
  }
  const double   E    = atof(argv[1]);
  const long     N    = atol(argv[2]);
  const uint64_t seed = strtoull(argv[3], nullptr, 10);

  RandomGen::rndm()->SetSeed(seed);

  const double field   = 96.5;                 // V/cm, LZ_Detector_2024::FitEF
  // 409.4312 e-. q.electrons is an integer, so "e < NE_CEIL" admits every
  // count up to and including 409: what this measures is P(N_e <= 409). The
  // analytic side must therefore evaluate its CDF at the half-integer edge,
  // 409.5, which is what analysis/response.py does. Do not "tidy" either
  // value to 409 -- fitting K at one ceiling and applying it at another
  // costs 1.75% on every absolute number and does not cancel.
  const double NE_CEIL = pow(10., 4.15)/34.5;
  const double NPH_MIN = 500./0.1104;          // 4529
  const double NPH_MAX = 600./0.1104;          // 5435

  LZ_Detector_2024 det;
  det.SetEFT2024Configuration();
  NESTcalc n(&det);

  bool inGas = false;
  const double rho = NESTcalc::GetDensity(175.9, 1.86, inGas);

  YieldResult y = n.GetYields(NEST::beta, E, rho, field, 131.293, 54.,
                              det.get_henr_nr_yield_params(),
                              det.get_henr_er_yield_params());
  vector<double> P = det.get_henr_nr_er_width_params();

  const double Nq = y.PhotonYield + y.ElectronYield;

  long nc = 0, nw = 0;
  double s1 = 0., s2 = 0., s3 = 0.;
  for (long i = 0; i < N; ++i) {
    QuantaResult q = n.GetQuanta(y, rho, P, -999.);
    const double e = (double)q.electrons;
    s1 += e; s2 += e*e; s3 += e*e*e;
    if (e < NE_CEIL) {
      ++nc;
      if (q.photons > NPH_MIN && q.photons < NPH_MAX) ++nw;
    }
  }
  const double m  = s1/N;
  const double v  = s2/N - m*m;
  const double sd = sqrt(v);
  const double sk = (s3/N - 3.*m*v - m*m*m)/(sd*sd*sd);

  printf("%7.2f %8.1f %7.4f %9.2f %8.2f %8.4f %12ld "
         "%7ld %.5e %.2e %7ld %.5e %.2e\n",
         E, Nq, y.ExcitonRatio, m, sd, sk, N,
         nc, (double)nc/N, sqrt((double)(nc > 0 ? nc : 1))/N,
         nw, (double)nw/N, sqrt((double)(nw > 0 ? nw : 1))/N);
  return 0;
}
