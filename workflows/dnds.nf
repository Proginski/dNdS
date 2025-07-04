/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	WELCOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info ""
log.info ">>> dNdS workflow <<<"
log.info ""
// log.info "This workflow runs codeml with the required ctl file as model."
// log.info "It is currently intended to receive a custom tree for each run."
// log.info ""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	WELCOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
	// log.info ">>> dnds <<<"
	log.info "Perfoms a dN/dS analysis using CODEML from PAML package."
	log.info "Takes as input at least a FASTA A, some FASTA B(s) and a newick tree."
	log.info "If no ortholog files are provided, it will use reciprocal best hits (RBH) as orthologs."
	log.info paramsHelp("nextflow run proginski/dnds --A <A NUCL FASTA> --B <B NUCL FASTA(S)> --tree <NEWICK WITH A AND B(s)> --ctl_file <CODEML CTL FILE> --outdir <OUTDIR>")
	exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
// validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FROM_ORTHO_PAIRS_TO_CODEML_INPUTS } from '../subworkflows/from_ortho_pairs_to_codeml_inputs.nf'
include { CODEML                            } from '../modules/codeml.nf'
include { CODEML_1OMEGA_FIXED1_VS_1OMEGA	} from '../modules/codeml_1omega_fixed1_vs_1omega.nf'
include { CODEML_1OMEGA_VS_2OMEGA       	} from '../modules/codeml_1omega_vs_2omega.nf'
include { COLLECT_RESULTS_ONE_OMEGA         } from '../modules/collect_results_one_omega.nf'
include { COLLECT_RESULTS_LYSIN_STYLE       } from '../modules/collect_results_lysin_style.nf'
include { COLLECT_RESULTS_MODELS_0_1_2      } from '../modules/collect_results_models_0_1_2.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DNDS {

	take :
	ortho_ch

	main:

	// Produce of proper alignment file for each ORF and its orthologs.
	FROM_ORTHO_PAIRS_TO_CODEML_INPUTS( ortho_ch )

	// Provided a channel of paths with both trees and alignments, but ensuring matching tree and alignment are send to the same task.
	codeml_inputs = FROM_ORTHO_PAIRS_TO_CODEML_INPUTS.out
		.map( tuple -> [ tuple[1], tuple[2] ] ) // only keep tree and aln
		.buffer( size: params.codeml_batch_size, remainder: true )
		.map{ tuple -> tuple.flatten() }
	
	// Use the deisred codeml config file.
	if (params.style == "branch_models"){
		// CODEML_1OMEGA_FIXED1_VS_1OMEGA(
		CODEML_1OMEGA_VS_2OMEGA(
			file(params.ctl_file),
			codeml_inputs
		)
	}
	else {
		CODEML(
			file(params.ctl_file),
			codeml_inputs
		)
	}

	// Build a summary '.tsv' file.
	// Collect the results according to the required codeml 'style'.
	if (params.style == "one_omega"){
		COLLECT_RESULTS_ONE_OMEGA  ( CODEML.out.collect() )
		// COLLECT_RESULTS_MODELS_0_1_2 ( CODEML.out.collect() )
	}
	if (params.style == "lysin"){
		COLLECT_RESULTS_LYSIN_STYLE( CODEML.out.collect() )
	}
}
