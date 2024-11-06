process GET_ALIGNMENT_FASTAS {

	label 'dnds'

	publishDir "${params.outdir}/trees/", pattern: "*.nwk"

	//memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (4.GB) }
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
	maxRetries 1

	input:
		path A_fna,  stageAs : "A_fna/*"
		path B_fnas, stageAs : "B_fnas/*"
		path ortho_pairs_tsv,   stageAs : "ortho_pairs/*"
		path phylip_names
		path tree
		path orfs

	output:
		path("*_n_orthologs.fna"), emit: orthologs_fnas
        path("*.nwk"), emit: newicks

	"""
    # First only keep the ORFs that are in the list provided by the user to avoid unnecessary processing.
	mv ${A_fna} ${A_fna}.tmp
	faSomeRecords ${A_fna}.tmp $orfs ${A_fna}

    sed -i "s/\$/\\\\t/" $orfs
	for tsv in $ortho_pairs_tsv
	do
		mv \$tsv \${tsv}.tmp
		grep -f $orfs \${tsv}.tmp > \$tsv || true
	done

    # Then,
	get_alignment_fastas.py --fna_a $A_fna --fna_b ${B_fnas} --ortho ${ortho_pairs_tsv} --name_mapping $phylip_names --tree $tree
	"""

}