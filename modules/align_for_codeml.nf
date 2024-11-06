process ALIGN_FOR_CODEML {

	label 'dnds'

	publishDir "${params.outdir}/alignments/", pattern: "*.aln"

	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
	maxRetries 1

	input:
		path orthologs_fnas

	output:
		path('*.aln'), optional : true

	"""
	
	"""

}