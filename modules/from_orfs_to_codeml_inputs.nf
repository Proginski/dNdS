process FROM_ORFS_TO_CODEML_INPUTS {

	publishDir "${params.outdir}/alignments/", pattern: "*.aln"
	publishDir "${params.outdir}/trees/", pattern: "*.nwk"

	errorStrategy 'retry'
    maxRetries 3

	input:
		val orf
        path aln_dir
        path tree_dir

	output:
		tuple val(orf), path("${orf}.nwk"), path("${orf}.aln"), optional : true

    """
    echo $orf
    if[ -f "${aln_dir}/${orf}.aln" ]
    then 
        echo "${aln_dir}/${orf}.aln exists"
    fi
    if [ -f "${tree_dir}/${orf}.nwk" ]
    then
        echo "${tree_dir}/${orf}.nwk exists"
    fi
	"""
}
