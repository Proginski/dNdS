process COLLECT_RESULTS_MODELS_0_1_2 {
	
	label 'local_job'

	publishDir "${params.outdir}/"

	input:
		path output_files

	output:
		path "results_*.tsv"

	"""
	single_out="results_models_0_1_2_\$(date -d "today" +"%Y%m%d_%H%M").tsv"

	echo "ORF	omega	p_neg1	p_neu1	w_neg1	w_neu1	p_neg2	p_neu2	p_pos2	w_neg2	w_neu2	w_pos2	np_0	lnL_0	np_1	lnL_1	np_2	lnL_2" > \$single_out

	for file in $output_files
	do
		codeml_parser_models_0_1_2.sh  \$file >> \$single_out
	done
	#codeml_parser_models_0_1_2.sh $output_files >> \$single_out
	"""

}