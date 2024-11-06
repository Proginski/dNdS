process COLLECT_RESULTS_LYSIN_STYLE {

	label 'local_job'

	publishDir "${params.outdir}/"

	input:
		path output_files

	output:
		path "results_lysin_style_*.tsv"

	"""
	awk '
        function basename(file) {
                sub(".*/", "", file)
                sub("\\.out", "", file)
                return file
        }

        BEGIN{	FS=OFS="\t"	}

		# The header.
		NR == 1 { print "ORF","p_neg","p_neu","p_pos","w_neg","w_neu","w_pos" }


		# For each file, on the same line print :

		# the ORF, 
        FNR==1 { printf basename(FILENAME)"\t" }

		# the percentages of codons under negative/neutral/positive selection,
        /^p:/ {printf gensub(/p:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\\\1\t\\\\2\t\\\\3","g",\$0)"\t"}

		# the omega values associated to these negative/neutral/positive selections.
        /^w:/ {print gensub(/w:[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+)[[:space:]]+([0-9.]+).*/,"\\\\1\t\\\\2\t\\\\3","g",\$0)}

	' $output_files > results_lysin_style_\$(date -d "today" +"%Y%m%d_%H%M").tsv
	"""

}