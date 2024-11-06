process CODEML {

	label 'dnds'
	
	maxRetries 1
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

	publishDir "${params.outdir}/codeml_out/"

	input:
		path ctl_file
		tuple val(orf), path(tree), path(aln)

	output:
		path "${orf}.out", optional : true

	"""
	# Create a codeml control file for the alignment
	sed "s~__ALN__~${orf}.aln~ ; s~__NWK__~${orf}.nwk~ ; s~__OUT__~${orf}.out~" ${ctl_file} > ${orf}.ctl

	# Here we go
	codeml ${orf}.ctl	
	"""

}
