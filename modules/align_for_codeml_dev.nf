process ALIGN_FOR_CODEML {

	label 'dnds'

	publishDir "${params.outdir}/alignments/", pattern: "*.aln"
	publishDir "${params.outdir}/trees/", pattern: "*.nwk"

	maxRetries 1
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

	input:
		path A_fna,  stageAs : "A_fna/*"
		path B_fnas, stageAs : "B_fnas/*"
		path ortho_pairs_tsv,   stageAs : "ortho_pairs/*"
		path phylip_names
		path tree
		val orfs

	output:
		tuple env(orf), path('env(orf).nwk'), path('env(orf).aln'), optional : true //path("${orf}.nwk"), path("${orf}.aln"), optional : true

	shell:
	'''
	# Turn groovy variables to bash variables
	focal=$(echo !{A_fna} | sed 's~.*/\\(.*\\)\\..*~\\1~')
	phylip_names=!{phylip_names}
	tree=!{tree}
	tab=$'\t'

	echo "focal = ${focal}"
	echo "phylip_names = ${phylip_names}"
	echo "tree = ${tree}"

	get_alignment_fastas.py --fna_a !{A_fna} --fna_b !{B_fnas} --otho !{ortho_pairs_tsv} --name_mapping $phylip_names -tree $tree
	'''

}