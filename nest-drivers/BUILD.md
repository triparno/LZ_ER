# Building

Clone NEST and check out commit 541b643. The detector class is
`LZ_Detector_2024` in `include/Detectors/LZ_WS2024.hh`; the configuration is
selected by `SetEFT2024Configuration()`.

Place lz_tail_test.cpp and lz_scan.cpp in the NEST source tree, add them as
targets, then

    cd build && cmake .. && make lz_tail_test lz_scan

## Why this commit and not a release

The Letter's main text cites NEST v2.4.5, its supplement discusses v2.4.0, and
LZ_WS2024.hh recommends v2.4.0. None reproduces the EFT2024 numbers: the
parameter vector has a component centre above unity, which selects the
log10(Nq) two-component branch of RecombOmegaER that earlier releases lack.

Commit 307fb27 (27 Aug 2026) added the configuration; 541b643 is what was used.

## Seeds

lz_scan.cpp takes energy, trial count and seed. The seed argument is required:
NEST's RandomGen is a singleton with deterministic initial state, so parallel
processes without distinct seeds return identical output.
