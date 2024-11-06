process TREE_FOR_CODEML {

	label 'local_job'

	publishDir "${params.outdir}/"

	input:
		path in_tree
		path phylip_names

	output:
		path "complete_tree.nwk"

	"""
	get_phylip_named_tree.py -tree $in_tree -names $phylip_names -out complete_tree.nwk
	sed -i -E "s/:[0-9.]+//g" complete_tree.nwk # Remove banch lengths
	"""
}