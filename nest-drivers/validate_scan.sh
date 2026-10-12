# Validation shared by the two scan runners. Sourced, not executed.
#
#   validate_scan OUTPUT NCOL EXPECTED
#
#     OUTPUT    the candidate scan file, header line first
#     NCOL      the number of fields every data row must carry
#     EXPECTED  a file of "energy trials" lines, one per job
#
# A scan is expensive and is committed to the archive, so a run that fails
# must not be able to replace a good file with a bad one. The row count alone
# does not establish that: a driver can exit non-zero and still print a line,
# or print the right number of rows carrying nan, the wrong energies, or
# fewer trials than it was asked for. Every one of those produces a file that
# looks structurally correct and is wrong.
#
# This checks, in order: the field count on every row, that every field is a
# finite number, that the energies are exactly the set the job list asked for
# with no duplicates and none missing, and that the trial count the driver
# reports in column 7 is the one it was given. It reports every failure it
# finds rather than stopping at the first, and returns non-zero if there was
# one. The caller is responsible for checking the exit status of the jobs
# themselves; this sees only their output.
validate_scan() {
    local out="$1" ncol="$2" expect="$3"
    awk -v NCOL="$ncol" '
    NR == FNR { want[$1+0] = $2+0; nwant++; next }
    /^[[:space:]]*#/ { next }
    /^[[:space:]]*$/ { next }
    {
        nrow++
        e = $1 + 0
        if (NF != NCOL)
            problem(sprintf("E=%s: %d fields, expected %d", $1, NF, NCOL))
        for (i = 1; i <= NF; i++)
            if ($i !~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/)
                problem(sprintf("E=%s: field %d is not a finite number: %s",
                                $1, i, $i))
        if (!(e in want))
            problem(sprintf("energy %s is not in the job list", $1))
        else if (seen[e]++)
            problem(sprintf("energy %s appears more than once", $1))
        else if (NF >= 7 && ($7 + 0) != want[e])
            problem(sprintf("E=%s: driver reports %s trials, the job list asked for %d",
                            $1, $7, want[e]))
    }
    END {
        if (nrow != nwant)
            problem(sprintf("%d data rows, expected %d", nrow, nwant))
        for (e in want)
            if (!(e in seen)) problem(sprintf("energy %g is missing", e))
        exit (bad > 0)
    }
    function problem(msg) { printf "  %s\n", msg > "/dev/stderr"; bad++ }
    ' "$expect" "$out"
}
