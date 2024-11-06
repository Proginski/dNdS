process COLLECT_RESULTS_ONE_OMEGA {

	label 'local_job'

	publishDir "${params.outdir}/"

	input:
		path output_files

	output:
		path "results_one_omega_*.tsv"

	"""
	awk '
		function basename(file) {
                sub(".*/", "", file)
                sub("\\.out", "", file)
                return file
        }

        BEGIN{ OFS="\t"	}

		# The header.
		NR == 1 { print "ORF","omega" }

		# For each file, on the same line print :

		# the ORF, 
        FNR==1 { printf basename(FILENAME)"\t" }
		
		\$1 == "omega" {print \$NF}
		
		' $output_files > results_one_omega_\$(date -d "today" +"%Y%m%d_%H%M").tsv
	"""

}