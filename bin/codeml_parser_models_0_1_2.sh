#!/bin/bash

awk -v header="${2}" '

# Only get the name of the ORF
function basename(file) {
                sub(".*/", "", file)
                sub("\\.out", "", file)
                return file
        }

BEGIN {
	OFS="\t"
	if(header == "--header"){ print "ORF","omega","p_neg1","p_neu1","w_neg1","w_neu1","p_neg2","p_neu2","p_pos2","w_neg2","w_neu2","w_pos2","np_0","lnL_0","np_1","lnL_1","np_2","lnL_2" }
}

$0 ~ /^Model 0:/ { model="zero" }
$0 ~ /^Model 1:/ { model="one"  }
$0 ~ /^Model 2:/ { model="two"  }

# Store the number of parameters and the lnL of each model
/^lnL/ {
	np[model]=gensub(/.*np:[[:space:]]*([0-9]+)\):[[:space:]]*([^[:space:]]+).*/,"\\1","g"$0)
	lnL[model]=gensub(/.*np:[[:space:]]*([0-9]+)\):[[:space:]]*([^[:space:]]+).*/,"\\2","g"$0)
}

# Get the particular metrics of each model

model == "zero" && $1 == "omega" {omega=$NF}
# the percentages of codons under negative/neutral/positive selection,
model == "one" && /^p:/ {
	p_neg[model]=gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\1","g",$0)
        p_neu[model]=gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\2","g",$0)
}
# the omega values associated to these negative/neutral/positive selections.
model == "one" && /^w:/ {
	w_neg[model]=gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\1","g",$0)
        w_neu[model]=gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\2","g",$0)
}
# the percentages of codons under negative/neutral/positive selection,
model == "two" && /^p:/ {
        p_neg[model]=gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\1","g",$0)
        p_neu[model]=gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\2","g",$0)
        p_pos[model]=gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\3","g",$0)
}
# the omega values associated to these negative/neutral/positive selections.
model == "two" && /^w:/ {
        w_neg[model]=gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\1","g",$0)
        w_neu[model]=gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\2","g",$0)
        w_pos[model]=gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\3","g",$0)
}


END {
	print basename(FILENAME), omega, p_neg["one"], p_neu["one"], w_neg["one"], w_neu["one"], p_neg["two"], p_neu["two"], p_pos["two"], w_neg["two"], w_neu["two"], w_pos["two"], np["zero"], lnL["zero"], np["one"], lnL["one"], np["two"], lnL["two"]
}
' $1
