# scan24_results.txt

Whitespace separated, one header line, 24 data rows sorted by energy.

    0   E        parent energy                      [keVee]
    1   Nq       total quanta
    2   NexNi    exciton-to-ion ratio
    3   meanNe   mean electrons
    4   sdNe     sigma_Ne
    5   skewNe   event-level skewness (gamma_ER)
    6   N        trials at this point
    7   n_ceil   hits with Ne < 409
    8   P_ceil   n_ceil / N
    9   dP_ceil  uncertainty on P_ceil
    10  n_win    hits also passing 500 < S1c < 600 phd
    11  P_win    n_win / N
    12  dP_win   uncertainty on P_win

fig_leakage_scan.py reads columns 0, 8, 9, 10, 11, 12.

Total budget 1.084e10 trials; 935 leakage counts. The 79 keV point is
deliberately under-sampled, being at the window edge where the bound is
kinematic rather than measured.

---

## scan24_smeared.txt

Produced by `nest-drivers/run_scan24_smeared.sh`, which carries the same
energies, trial counts and seeds through `GetS1` and `GetS2` and applies the
cuts in phd rather than in quanta. Columns 0-12 are as above; the quanta-level
counters are recomputed on the same events, so this file is self-contained.

    13  n_ceil_s   hits with S2c < 10^4.15 phd
    14  P_ceil_s   n_ceil_s / N
    15  dP_ceil_s  uncertainty on P_ceil_s
    16  n_win_s    hits also passing 500 < S1c < 600 phd
    17  P_win_s    n_win_s / N
    18  dP_win_s   uncertainty on P_win_s
    19  n_slow     events that needed the full GetS1/GetS2 treatment

`read_smeared_scan.py` reads 0, 6, 7, 8, 10, 11, 13, 14, 16, 17, 19.

The comparison is paired: columns 10-12 and 16-18 classify the *same*
simulated events two ways, so P_win_s/P_win carries no sampling error beyond
the counts. Columns 0-12 agree with scan24_results.txt statistically but not
bit-for-bit, because GetS1 and GetS2 draw from the same RandomGen singleton
and shift the GetQuanta sequence. As a check that this does not matter, the
window integral of column 11 reproduces the archived value to four figures:
1.4968e-6 against 1.4970e-6 keV.

**Columns 13-15 of the committed file are not measurements on every row.**
The prefilter that produced it ran the full detector treatment only for events
near the S1c window as well as the S2c ceiling, so for events that cannot
reach the window -- every event at 45, 88, 94, 104 and 110 keVee, and some
elsewhere -- n_ceil_s is a copy of n_ceil rather than a measurement. Those
rows are recognisable by n_slow being zero or small beside n_ceil. P_win_s,
columns 16-18, is unaffected: an event more than eight standard deviations
outside the window cannot enter it however the detector smears, so its
soft_win is zero either way, and an event inside both cuts was already
counted. Nothing in the paper reads columns 13-15. The driver in this archive
tests the charge condition alone, so a rerun measures them
properly, at the cost of a larger n_slow and the same wall time.
