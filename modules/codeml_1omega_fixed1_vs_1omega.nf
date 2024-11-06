process CODEML_1OMEGA_FIXED1_VS_1OMEGA {
	
	label 'dnds'

	maxRetries 1
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

	publishDir "${params.outdir}/codeml_out/"

	input:
		path ctl_file
		tuple val(orf), path(tree), path(aln)

	output:
		path "${orf}_models.tsv", optional : true

	"""
	# Run several models of codeml for the given alignment

	# Store their likelihoods and number of parameters in a file
	tab=\$'\t'
	echo "model\${tab}tree\${tab}lnL\${tab}np" > ${orf}_models.tsv

	for ctl in 1omega 1omega_fixed1
	do
		echo \$ctl

		# Create a codeml control file for the alignment
		sed "s~__ALN__~${orf}.aln~ ; s~__NWK__~${orf}.nwk~ ; s~__OUT__~${orf}_\${ctl}.out~" ${projectDir}/ctl/\${ctl}.ctl > ${orf}_\${ctl}.ctl

		# Here we go
		codeml ${orf}_\${ctl}.ctl

		# Write the np and lnL values of the model
		codeml_stats.sh ${orf}_\${ctl}.out ${orf}.nwk >> ${orf}_models.tsv

	done

	compare_codeml_models.py ${orf}_models.tsv

	detect_selection_with_compared_models.py ${orf}_models.tsv --pvalue 0.05
	"""

}
