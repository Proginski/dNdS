#!/bin/bash

# Deal with LINEAR FASTA FILES ONLY
awk '

    BEGIN {
        OFS = "\t"
        print "FASTA", "isHomogeneous"
    }
    BEGINFILE {
        homo = "TRUE"
        seq  = ""
    }

    /^[^>]/ && $0 != seq {
        if (seq == "") {
            seq = $0
        } else {
            homo = "FALSE"
            nextfile
        }
    }
    ENDFILE {
        if (seq != "") {
            print FILENAME, homo
    }
}
' "$@"