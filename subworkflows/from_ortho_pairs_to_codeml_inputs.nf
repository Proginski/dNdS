
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_FILES                 } from '../modules/check_files.nf'
include { NAMES_TO_PHYLIP             } from '../modules/names_to_phylip.nf'
include { TREE_FOR_CODEML             } from '../modules/tree_for_codeml.nf'
include { GET_SEQUENCES               } from '../modules/get_sequences.nf'
include { GET_ALIGNMENT_FASTAS        } from '../modules/get_alignment_fastas.nf'
include { ALIGN_FOR_CODEML            } from '../modules/align_for_codeml.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FROM_ORTHO_PAIRS_TO_CODEML_INPUTS {

	take:
	ortho_ch

	main:

	A_fna = file(params.A)
	B_fnas = Channel.fromPath(params.B)

	CHECK_FILES(
		A_fna,
		B_fnas.collect(), // .collect() turns the channel into one single list
		ortho_ch.collect()
		)
	B_fnas = CHECK_FILES.out.B_fnas
	ortho_ch = CHECK_FILES.out.ortho_pairs

	// Get a '.csv' file with orignal genome names as first column and their corresponding short_phylip name as second column.
	// FASTA files basename
	A_name_ch = Channel.fromPath(params.A)
		.map{ fasta -> fasta.getBaseName() }
	B_names_ch = Channel.fromPath(params.B)
		.map{ fasta -> fasta.getBaseName() }
	names_ch = A_name_ch.concat(B_names_ch).unique()
	NAMES_TO_PHYLIP(names_ch.collectFile(name: 'names.txt', newLine: true))

	// Convert the input newick to a tree without internal nodes and with PHYLIP compatible names.
	TREE_FOR_CODEML(
					file(params.tree),
					NAMES_TO_PHYLIP.out
					)
	// Each ORF from the list of provided by the user is processed separately.
	if (params.A_seqlist) {
		user_ORFs_list = file(params.A_seqlist)
	} else {
		user_ORFs_list = file('dummy')
	}
	ORFs_ch = GET_SEQUENCES( A_fna, user_ORFs_list)
		.splitText()
		.map{it -> it.trim()}

	ORFs_ch = ORFs_ch.collectFile(name: 'orfs.txt', newLine: true).splitText( by: 100, file: "orf_subset.txt")
	// ORFs_ch = ORFs_ch.buffer(size: params.alignment_batch_size, remainder: true)

	// Produce of proper alignment file for one ORF and its orthologs.
	GET_ALIGNMENT_FASTAS(
					 A_fna,
					 B_fnas.collect(),
					 ortho_ch.collect(),
					 NAMES_TO_PHYLIP.out.first(),
					 TREE_FOR_CODEML.out.first(),
					 ORFs_ch
					)

	// Produce of proper alignment file for each ORF and its orthologs.
	ALIGN_FOR_CODEML( GET_ALIGNMENT_FASTAS.out.orthologs_fnas )

	// OUTPUT
	emit : ALIGN_FOR_CODEML.out
}