#!/bin/bash

tree="NA"
if [ ! -z "$2" ]  # Check if $2 is not empty
then
    tree=$(cat $2)
fi

awk -v tree="${tree}" '
    BEGIN{
        OFS="\t"
#		print "model","lnL","np"
    }

    /^lnL\(ntime/ {
        $0=gensub(/.*np:/,"","g",$0)
        $0=gensub(/\):/,"","g",$0)
        print FILENAME,tree,$2,$1
    }
' $1
