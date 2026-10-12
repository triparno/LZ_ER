# The 24-point leakage scan

scan24_results.txt holds the scan of Sec. 2: 1.084e10 trials across 24
energies from 45 to 110 keVee, 935 leakage counts. Column format in
COLUMNS.md.

It is the only artefact here that is not regenerable in minutes. To rebuild it:

    cd ../nest-drivers && ./run_scan24.sh 8     # ~2 h on 8 cores

A rerun replaces this file only if it validates: the runner checks every
job's exit status, the energy grid, the field count, that every field is a
finite number, and that the trials each row reports are the trials it was
asked for. On any failure the committed file is left alone and the candidate
is kept beside it. nest-drivers/test_runner_guards.sh exercises those checks.

The run is deterministic -- run_scan24.sh fixes seeds 201-224 and NEST's
RandomGen is a singleton -- so a rerun reproduces the committed numbers
exactly rather than only statistically. A rebuild should return
K = 0.937 +/- 0.031 against the analytic representation, with a Poisson
deviance of 20.8/23 = 0.91 (p = 0.59); fit_K.py regenerates all three,
the values quoted in the abstract and in Sec. 2.

analysis/fig_leakage_scan.py reads this file; set SCAN_FILE at the top of that
script to its path.
