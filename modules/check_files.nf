process CHECK_FILES {

	label 'local_job'

	input:
		path A_fna, stageAs : "A_fna/*"
		path B_fnas, stageAs : "B_fnas/*"
		path(ortho_pairs, stageAs : "ortho_pairs/*")

	output:
		path "B_fnas_checked/*",   emit: B_fnas
		path "ortho_pairs_checked/*_vs_*", emit: ortho_pairs

	shell = ['/bin/bash','-u']
	"""
	#### FASTA files ####
	echo -n "Checking FASTA files..."
	A_name=\$(echo ${A_fna} | sed 's~.*/\\(.*\\)\\..*~\\1~')
	if [ ! -s ${A_fna} ] 
		then
			echo ""
			echo "ERROR : ${A_fna} does not exists or is empty."
			exit 1
	fi
	for fasta in ${B_fnas}
	do
		if [ ! -s \${fasta} ] 
		then
			echo ""
			echo "ERROR : \${fasta} does not exists or is empty."
			exit 1
		fi

		name=\$(echo \${fasta} | sed 's~.*/\\(.*\\)\\..*~\\1~')
		if [[ \${name} == \${A_name} ]]
		then
			rm -f \$fasta
			continue
		fi
		echo "\${A_name}_vs_\${name}_orthologs.tsv" >> expected_ortho_pairs.txt
	done
	echo " Done."
	#####################


	#### Ortholog pairs files ####
	echo -n "Checking ortho pairs files..."

	# Test is ortho pairs tsv have OrthoFinder style names
	ls ortho_pairs/ | if grep -q ".__v__.*\\.tsv"
	then
		echo ""
		echo "OrthoFinder style names detected among ortholog pairs files."
		echo -n "Renaming ortholog pairs files..."
		ls ortho_pairs/ | while read tsv
		do
			new_tsv=\$(echo \${tsv} | sed "s~__v__~_vs_~ ; s~\\.tsv\$~_orthologs.tsv~ ")
			OrthoFinder_orthologues_to_2columns.py ortho_pairs/\${tsv} ortho_pairs/\${new_tsv}
		done
		echo " Done."
	fi

	# Test if we have all expected ortholog pairs files
	ls ortho_pairs/ > received_ortho_pairs.txt
	grep -vxf received_ortho_pairs.txt expected_ortho_pairs.txt > missing_ortho_pairs.txt
	if [ -s missing_ortho_pairs.txt ]
	then
		echo "ERROR : missing ortholog pairs tsv files (see 'missing_ortho_pairs.txt'):"
		cat missing_ortho_pairs.txt
		exit 1
	fi
	echo " Done."
	###############################


	mv B_fnas/ B_fnas_checked/
	mv ortho_pairs/ ortho_pairs_checked/
	"""
}