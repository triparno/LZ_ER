// lz_scan_s1s2.cpp -- the leakage scan carried through to S1c and S2c, so that
// the photon-detection and electron-extraction fluctuations omitted by
// lz_scan.cpp are included. One energy point per process.
//
//   probe mode:  ./lz_scan_s1s2 probe <E_keV> [N]
//   scan  mode:  ./lz_scan_s1s2 <E_keV> <N_trials> <seed> <s1c_idx> <s2c_idx>
//
// RUN THE PROBE FIRST. The layout of the vectors returned by GetS1 and GetS2
// differs between NEST versions, and this driver will not guess: in scan mode
// both indices are required arguments. Probe mode prints every element of both
// vectors with its mean, its spread and its ratio to the expected corrected
// pulse area, so the right index is read off in seconds. Pick the element whose
// mean is within a few per cent of g1*<Nph> (resp. g2*<Ne>) AND which has a
// non-zero spread and non-integer values: that is the corrected pulse area in
// phd. Integer-valued entries are hit or spike counts, not areas.
//
// Scan mode prints one line with the quanta-level counters of lz_scan.cpp and
// the S1c/S2c counters side by side, on the same events:
//
//   E Nq NexNi meanNe sdNe skewNe N n_ceil P_ceil dP_ceil n_win P_win dP_win
//                              n_ceil_s P_ceil_s dP_ceil_s n_win_s P_win_s dP_win_s n_slow
//
// The paired comparison within one run is exact: the same simulated events are
// classified both ways, so P_win_s/P_win measures the smearing with no
// sampling error of its own beyond the counts. The quanta-level columns will
// agree with scan24_results.txt statistically but not bit-for-bit: GetS1 and
// GetS2 draw from the same RandomGen singleton, so the sequence of GetQuanta
// draws differs from the quanta-only run.
//
// Cost. A full S1/S2 treatment of every trial would take days. It is not
// needed: an event can only be reclassified by the smearing if it lies near a
// boundary, and at a leakage of ~1e-7 almost none do. An event whose charge
// exceeds the S2c ceiling by more than MARGIN_SIG standard deviations cannot
// satisfy either cut however the detector smears it, so it is classified by
// the hard cut, which is exact far below the statistics of the run. The
// number that took the slow path is reported as n_slow so the saving is
// auditable.
//
// The prefilter tests the charge condition ONLY, and must continue to. Adding
// the S1c window to it would make P_ceil_s a copy of P_ceil for any event
// that cannot reach the window: harmless for P_win_s, which is what the paper
// uses, but it leaves the ceiling-only column looking measured where it is
// not.
//
// The signatures used here are those of commit 541b643, verified against
// include/NEST/NEST.hh lines 335 (GetS1), 351 (GetS2) and 368
// (SetDriftVelocity). Both getters return a const reference to an internal
// buffer; the copies below are deliberate, since GetS2 is called while the
// S1 result is still needed. Neither takes an outputTiming flag, and
// SetDriftVelocity takes the pressure as a fourth argument.
//
// Positions are fixed at the detector axis and mid-drift, which removes the
// position dependence of the corrections rather than modelling it. The scan
// asks what the detection fluctuations do to a response, not what the position
// distribution does.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>
#include "NEST.hh"
#include "RandomGen.hh"
#include "VDetector.hh"
#include "LZ_WS2024.hh"

using namespace std;
using namespace NEST;

static const double FIELD   = 96.5;                  // V/cm, LZ_Detector_2024
static const double G1      = 0.1104;                // phd/photon
static const double G2      = 34.5;                  // phd/electron
static const double S1_LO   = 500.,  S1_HI = 600.;   // phd, comparison region
static const double S2_CEIL = pow(10., 4.15);        // phd, WS RoI ceiling
static const double S2_MIN  = 645.;                  // phd, WS RoI raw S2 floor
static const double NE_CEIL = S2_CEIL/G2;            // 409.4312 e-; on an
// integer count this selects N_e <= 409, matching lz_scan.cpp. The soft cut
// below is on the CONTINUOUS corrected area S2c, where no such correction
// applies -- the two are different quantities and both are right as written.
static const double NPH_MIN = S1_LO/G1, NPH_MAX = S1_HI/G1;
static const double MARGIN_SIG = 8.0;                // slow-path band, in sigma

struct Geom {
  double rho, driftV, zmid, dt;
  vector<double> g2par;
};

static Geom setup(NESTcalc& n, LZ_Detector_2024& det) {
  Geom g;
  bool inGas = false;
  g.rho    = NESTcalc::GetDensity(175.9, 1.86, inGas);
  g.driftV = n.SetDriftVelocity(175.9, g.rho, FIELD, 1.86);
  g.zmid   = 0.5*(det.get_TopDrift() + det.get_cathode());
  g.dt     = (det.get_TopDrift() - g.zmid)/g.driftV;
  g.g2par  = n.CalculateG2(0);
  return g;
}

// one event: quanta, then S1 and S2 if it is near a boundary
static inline void simulate(NESTcalc& n, const Geom& g, const YieldResult& y,
                            const vector<double>& P, double Nq, uint64_t ev,
                            double E, int i1, int i2,
                            int& hard_ceil, int& hard_win,
                            int& soft_ceil, int& soft_win, long& nslow,
                            double& Ne_out) {
  QuantaResult q = n.GetQuanta(y, g.rho, P, -999.);
  const double ne  = (double)q.electrons;
  const double nph = (double)q.photons;
  Ne_out = ne;

  hard_ceil = (ne < NE_CEIL);
  hard_win  = (hard_ceil && nph > NPH_MIN && nph < NPH_MAX);
  soft_ceil = hard_ceil; soft_win = hard_win;

  // can the smearing move this event across a boundary?
  const double sdS1 = sqrt(nph*G1) ;              // phd, counting term
  const double sdS2 = sqrt(ne > 1. ? ne : 1.)*G2*0.3;  // phd, generous
  const bool nearCeil = fabs(G2*ne - S2_CEIL) < MARGIN_SIG*sdS2;
  const bool nearWin  = fabs(G1*nph - S1_LO)  < MARGIN_SIG*sdS1 ||
                        fabs(G1*nph - S1_HI)  < MARGIN_SIG*sdS1;
  // The ceiling counter must not depend on the S1c condition, or P_ceil_s
  // would be a copy of the hard value for any event that cannot reach the
  // window -- a number that looks measured and is not. The S2 condition alone
  // decides whether the slow path runs; an event failing it has
  // N_e > NE_CEIL + 8 sigma, so hard_ceil and hard_win are both already zero
  // and copying them is correct.
  const bool couldPass = (G2*ne < S2_CEIL + MARGIN_SIG*sdS2);
  if (!couldPass) return;
  (void)nearCeil; (void)nearWin;

  ++nslow;
  vector<int64_t> wt; vector<double> wa;
  vector<double> s1 = n.GetS1(q, 0., 0., g.zmid, 0., 0., g.zmid,
                              g.driftV, g.driftV, NEST::beta, ev, FIELD, E,
                              S1CalculationMode::Parametric, 0, wt, wa);
  vector<double> s2 = n.GetS2(q.electrons, 0., 0., g.zmid, 0., 0., g.zmid,
                              g.dt, g.driftV, ev, FIELD,
                              S2CalculationMode::Full, 0,
                              wt, wa, g.g2par);
  if ((int)s1.size() <= i1 || (int)s2.size() <= i2) return;
  const double S1c = s1[i1], S2c = s2[i2];
  const double S2raw = (i2 > 0 ? s2[i2-1] : S2c);
  soft_ceil = (S2c < S2_CEIL && S2raw > S2_MIN);
  soft_win  = (soft_ceil && S1c > S1_LO && S1c < S1_HI);
}

static int probe(double E, long N) {
  LZ_Detector_2024 det; det.SetEFT2024Configuration();
  NESTcalc n(&det);
  Geom g = setup(n, det);
  YieldResult y = n.GetYields(NEST::beta, E, g.rho, FIELD, 131.293, 54.,
                              det.get_henr_nr_yield_params(),
                              det.get_henr_er_yield_params());
  vector<double> P = det.get_henr_nr_er_width_params();

  vector<double> m1s, m2s, m1i, m2i;
  vector<int> inti, intj;
  double mnph = 0., mne = 0.;
  for (long k = 0; k < N; ++k) {
    QuantaResult q = n.GetQuanta(y, g.rho, P, -999.);
    mnph += q.photons; mne += q.electrons;
    vector<int64_t> wt; vector<double> wa;
    vector<double> s1 = n.GetS1(q, 0., 0., g.zmid, 0., 0., g.zmid,
                                g.driftV, g.driftV, NEST::beta, (uint64_t)k,
                                FIELD, E,
                                S1CalculationMode::Parametric, 0, wt, wa);
    vector<double> s2 = n.GetS2(q.electrons, 0., 0., g.zmid, 0., 0., g.zmid,
                                g.dt, g.driftV, (uint64_t)k, FIELD,
                                S2CalculationMode::Full, 0,
                                wt, wa, g.g2par);
    if (m1s.empty()) { m1s.assign(s1.size(),0.); m2s.assign(s1.size(),0.);
                       inti.assign(s1.size(),1);
                       m1i.assign(s2.size(),0.); m2i.assign(s2.size(),0.);
                       intj.assign(s2.size(),1); }
    for (size_t j = 0; j < s1.size() && j < m1s.size(); ++j) {
      m1s[j] += s1[j]; m2s[j] += s1[j]*s1[j];
      if (fabs(s1[j]-llround(s1[j])) > 1e-9) inti[j] = 0;
    }
    for (size_t j = 0; j < s2.size() && j < m1i.size(); ++j) {
      m1i[j] += s2[j]; m2i[j] += s2[j]*s2[j];
      if (fabs(s2[j]-llround(s2[j])) > 1e-9) intj[j] = 0;
    }
  }
  mnph /= N; mne /= N;
  printf("\nprobe at E = %.2f keV, %ld events\n", E, N);
  printf("  <Nph> = %.1f   ->  g1*<Nph> = %.2f phd   (expected S1c)\n", mnph, G1*mnph);
  printf("  <Ne>  = %.1f   ->  g2*<Ne>  = %.2f phd   (expected S2c)\n\n", mne, G2*mne);
  printf("  GetS1 returns %zu elements:\n", m1s.size());
  printf("   %3s %14s %12s %10s %8s\n", "idx", "mean", "sd", "/g1<Nph>", "integer");
  for (size_t j = 0; j < m1s.size(); ++j) {
    double m = m1s[j]/N, v = m2s[j]/N - m*m;
    printf("   %3zu %14.4f %12.4f %10.4f %8s%s\n", j, m, v > 0 ? sqrt(v) : 0.,
           G1*mnph != 0. ? m/(G1*mnph) : 0., inti[j] ? "yes" : "no",
           (!inti[j] && fabs(m/(G1*mnph)-1.) < 0.10) ? "   <-- candidate s1c_idx" : "");
  }
  printf("\n  GetS2 returns %zu elements:\n", m1i.size());
  printf("   %3s %14s %12s %10s %8s\n", "idx", "mean", "sd", "/g2<Ne>", "integer");
  for (size_t j = 0; j < m1i.size(); ++j) {
    double m = m1i[j]/N, v = m2i[j]/N - m*m;
    printf("   %3zu %14.4f %12.4f %10.4f %8s%s\n", j, m, v > 0 ? sqrt(v) : 0.,
           G2*mne != 0. ? m/(G2*mne) : 0., intj[j] ? "yes" : "no",
           (!intj[j] && fabs(m/(G2*mne)-1.) < 0.10) ? "   <-- candidate s2c_idx" : "");
  }
  printf("\n  Pass the two chosen indices to scan mode. Where several are marked,\n"
         "  take the one whose ratio is CLOSEST TO ONE: g1 and g2 are defined so\n"
         "  that the corrected area equals g*N, and the raw area is the one that\n"
         "  sits a per cent or so off. Do not simply take the later index. In\n"
         "  Parametric S1 mode the spike counts are set equal to the phd values,\n"
         "  so a later index may be a copy of the right one here and a different\n"
         "  quantity under S1CalculationMode::Full.\n"
         "\n  The two vectors do not share a layout: GetS2 carries the extracted\n"
         "  electron count and the S2 photon count ahead of the areas, so its\n"
         "  indices run two higher than GetS1's.\n"
         "\n  At commit 541b643 the answer is s1c_idx = 5, s2c_idx = 7.\n\n");
  return 0;
}

int main(int argc, char** argv) {
  if (argc >= 2 && strcmp(argv[1], "probe") == 0) {
    if (argc < 3) { fprintf(stderr, "usage: %s probe <E_keV> [N]\n", argv[0]); return 1; }
    RandomGen::rndm()->SetSeed(1);
    return probe(atof(argv[2]), argc > 3 ? atol(argv[3]) : 2000L);
  }
  if (argc < 6) {
    fprintf(stderr,
      "usage: %s <E_keV> <N_trials> <seed> <s1c_idx> <s2c_idx>\n"
      "       %s probe <E_keV> [N]        <- run this first\n", argv[0], argv[0]);
    return 1;
  }
  const double   E    = atof(argv[1]);
  const long     N    = atol(argv[2]);
  const uint64_t seed = strtoull(argv[3], nullptr, 10);
  const int      i1   = atoi(argv[4]);
  const int      i2   = atoi(argv[5]);

  RandomGen::rndm()->SetSeed(seed);
  LZ_Detector_2024 det; det.SetEFT2024Configuration();
  NESTcalc n(&det);
  Geom g = setup(n, det);
  YieldResult y = n.GetYields(NEST::beta, E, g.rho, FIELD, 131.293, 54.,
                              det.get_henr_nr_yield_params(),
                              det.get_henr_er_yield_params());
  vector<double> P = det.get_henr_nr_er_width_params();
  const double Nq = y.PhotonYield + y.ElectronYield;

  long nc = 0, nw = 0, ncs = 0, nws = 0, nslow = 0;
  double s1 = 0., s2 = 0., s3 = 0.;
  for (long i = 0; i < N; ++i) {
    int hc, hw, sc, sw; double ne;
    simulate(n, g, y, P, Nq, (uint64_t)i, E, i1, i2, hc, hw, sc, sw, nslow, ne);
    s1 += ne; s2 += ne*ne; s3 += ne*ne*ne;
    nc += hc; nw += hw; ncs += sc; nws += sw;
  }
  const double m  = s1/N;
  const double v  = s2/N - m*m;
  const double sd = sqrt(v);
  const double sk = (s3/N - 3.*m*v - m*m*m)/(sd*sd*sd);
  const double e = 1e-30;
  printf("%7.2f %8.1f %7.4f %9.2f %8.2f %8.4f %12ld "
         "%7ld %.5e %.2e %7ld %.5e %.2e "
         "%9ld %.5e %.2e %8ld %.5e %.2e %10ld\n",
         E, Nq, y.ExcitonRatio, m, sd, sk, N,
         nc,  (double)nc/N,  sqrt((double)(nc  > 0 ? nc  : 1))/N,
         nw,  (double)nw/N,  sqrt((double)(nw  > 0 ? nw  : 1))/N,
         ncs, (double)ncs/N, sqrt((double)(ncs > 0 ? ncs : 1))/N,
         nws, (double)nws/N, sqrt((double)(nws > 0 ? nws : 1))/N,
         nslow);
  (void)e;
  return 0;
}
