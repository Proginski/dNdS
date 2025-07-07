process CODEML_1OMEGA_VS_2OMEGA {
	
	label 'dnds'

	maxRetries 1
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

	publishDir "${params.outdir}/codeml_out/models_raw"

	input:
		path ctl_file
		path trees_and_alns

	output:
		path "*_models.tsv", optional : true

	script:
	"""
	ls *.nwk | sed "s/.nwk\$//" > orfs.txt

	cat orfs.txt | while read orf
	do
		echo \$orf

		aln=\${orf}_for_codeml.aln
		tree=\${orf}.nwk

		# Run several models of codeml for the given alignment

		# Store their likelihoods and number of parameters in a file
		tab=\$'\t'
		echo "model\${tab}tree\${tab}lnL\${tab}np" > \${orf}_models.tsv

		for ctl in 1omega 1omega_fixed1
		do
			echo \$ctl

			# Create a codeml control file for the alignment
			sed "s~__ALN__~\${aln}~ ; s~__NWK__~\${orf}.nwk~ ; s~__OUT__~\${orf}_\${ctl}.out~" ${projectDir}/ctl/\${ctl}.ctl > \${orf}_\${ctl}.ctl

			# run codeml
			codeml \${orf}_\${ctl}.ctl

			# Write the np and lnL values of the model
			codeml_stats.sh \${orf}_\${ctl}.out \${orf}.nwk >> \${orf}_models.tsv
		done

		# From the alignment tree, generate mutiple trees.
		i=0
		# For each branch, a tree where this particular branch (foreground) has an omega different from the rest of the tree
		# For each subtree, a tree where the subtree has an omega different from the rest of the tree
		tree_two_omegas_combinaisons.py \$tree | while read tree_string
		do
			echo \$tree_string
			i=\$((i+1))

			# Write the tree to a file
			echo \$tree_string > \${orf}_2omega_\${i}.nwk

			# For each custom tree, run two models with 2 omegas 
			# one where the foureground omega is fixed to 1
			# one where the foreground omega is estimated

			for ctl in 2omega 2omega_fixed1
			do
				echo \$ctl

				# Create a codeml control file for the alignment
				sed "s~__ALN__~\${aln}~ ; s~__NWK__~\${orf}_2omega_\${i}.nwk~ ; s~__OUT__~\${orf}_\${ctl}_\${i}.out~ ; s~__TREE__~\$tree_string~" ${projectDir}/ctl/\${ctl}.ctl > \${orf}_\${ctl}_\${i}.ctl

				# run codeml
				codeml \${orf}_\${ctl}_\${i}.ctl

				# Write the np and lnL values of the model
				codeml_stats.sh \${orf}_\${ctl}_\${i}.out \${orf}_2omega_\${i}.nwk >> \${orf}_models.tsv
			done
		done

		# Finally run a free omega model
		echo "free omega"

		# Create a codeml control file for the alignment
		sed "s~__ALN__~\${aln}~ ; s~__NWK__~\${orf}.nwk~ ; s~__OUT__~\${orf}_free_omega.out~" ${projectDir}/ctl/free_omega.ctl > \${orf}_free_omega.ctl

		# run codeml
		codeml \${orf}_free_omega.ctl

		# Write the np and lnL values of the model
		codeml_stats.sh \${orf}_free_omega.out \${orf}.nwk >> \${orf}_models.tsv


		compare_codeml_models.py \${orf}_models.tsv

		detect_selection_with_compared_models.py \${orf}_models.tsv --pvalue 0.05

	done
	"""

}
