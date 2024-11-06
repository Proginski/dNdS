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

	for orf in !{orfs.join(' ')}
	do
		fa="\${orf}_n_orthologs.fna"
		echo "orf = \${orf}"
		
		# First make a fasta with the ORF and its orthologs

		# Begin with the focal's ORF
		# As name of each sequence in the aligment has to be the genome's name (PHYLIP compliant), add the focal's phylip name as the header
		phylip_focal=$(awk -v name="${focal}" -F"," '$1==name {print $2}' $phylip_names)
		echo ">${phylip_focal}" > $fa
		# Add the ORF's sequence
		faOneRecord !{A_fna} $orf | tail -n +2 >> $fa

		# Then for each neighbor,
		for B_fna in !{B_fnas}
		do
			name=$(echo $B_fna | sed 's~.*/\\(.*\\)\\..*~\\1~')
			phylip_name=$(awk -v name="${name}" -F"," '$1==name {print $2}' $phylip_names)
			echo "current neighbor phylip name : ${phylip_name=}"

			# If there is an ortholog,
			if grep -q "^${orf}${tab}" ortho_pairs/${focal}_vs_${name}_orthologs.tsv
			then
				# Get it
				seq_to_add=$(grep "^${orf}${tab}" ortho_pairs/${focal}_vs_${name}_orthologs.tsv | awk -F"\t" '{print $2}')
				echo "ortholog found : ${seq_to_add}"

				# Add the neighbor phylip_name as header
				echo ">${phylip_name=}" >> $fa

				# Find its sequence and add it to the alignment fasta
				#faOneRecord ${B_fna} $seq_to_add | tail -n +2 >> $fa
				getonlyseq.sh ${B_fna} $seq_to_add >> $fa
			fi
		done

		seq_nb=$(grep -c ">" $fa)

		# Then, if there is more than one entry in the alignment fasta,
		if [ $seq_nb -gt 1 ]
		then
			# Remove terminal stop codons
			linearizefasta.sh $fa > ${fa}.tmp
			sed -E "/>/! s/...$//" ${fa}.tmp > $fa
			rm ${fa}.tmp

			# Align both nucleotides and AAs with MAFFT (see '-p F')
			translatorX.pl -i $fa -p F -o $orf

			# Get codons aligned with respect to amino acids 
			pal2nal.pl ${orf}.aa_ali.fasta ${orf}.nt_ali.fasta -output fasta > ${orf}_pal2nal.aln

			# If the below line is removed, you must linearize before format conversion (linearizefasta.sh ${orf}_pal2nal.aln > ${orf}_pal2nal.aln.oneline)
			remove_mostly_gapped_positions.py ${orf}_pal2nal.aln ${orf}_clean.aln

			# Turn the FASTA alignment to PHYLIP format (input_file, number of sequences, number of columns)
			FASTAtoPHYL.pl ${orf}_clean.aln $seq_nb $(awk 'FNR==2 {print length}' ${orf}_clean.aln)
			mv ${orf}_clean.aln.phy ${orf}.aln

			# Generate the phylogenetic tree with only the species that have a match
			get_subtree.py -tree $tree -names <(grep ">" $fa | sed "s/>//") -out ${orf}.nwk
		fi


	done # end of for orf in !{orfs} ...
	'''

}




// process ALIGN_FOR_CODEML { // OLD WAY, no buffer

// 	label 'dnds'

// 	publishDir "${params.outdir}/alignments/", pattern: "*.aln"
// 	publishDir "${params.outdir}/trees/", pattern: "*.nwk"

// 	maxRetries 1
// 	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

// 	input:
// 		path A_fna,  stageAs : "A_fna/*"
// 		path B_fnas, stageAs : "B_fnas/*"
// 		path ortho_pairs_tsv,   stageAs : "ortho_pairs/*"
// 		path phylip_names
// 		path tree
// 		val orf

// 	output:
// 		tuple val(orf), path("${orf}.nwk"), path("${orf}.aln"), optional : true

// 	shell:
// 	'''
// 	# Turn groovy variables to bash variables
// 	focal=$(echo !{A_fna} | sed 's~.*/\\(.*\\)\\..*~\\1~')
// 	phylip_names=!{phylip_names}
// 	tree=!{tree}
// 	orf=!{orf}
// 	fa="${orf}_n_orthologs.fna"
// 	tab=$'\t'

// 	echo "focal = ${focal}"
// 	echo "phylip_names = ${phylip_names}"
// 	echo "tree = ${tree}"
// 	echo "orf = ${orf}"


// 	# First make a fasta with the ORF and its orthologs

// 	# Begin with the focal's ORF
// 	# As name of each sequence in the aligment has to be the genome's name (PHYLIP compliant), add the focal's phylip name as the header
// 	phylip_focal=$(awk -v name="${focal}" -F"," '$1==name {print $2}' $phylip_names)
// 	echo ">${phylip_focal}" > $fa
// 	# Add the ORF's sequence
// 	faOneRecord !{A_fna} $orf | tail -n +2 >> $fa

// 	# Then for each neighbor,
// 	for B_fna in !{B_fnas}
// 	do
// 		name=$(echo $B_fna | sed 's~.*/\\(.*\\)\\..*~\\1~')
// 		phylip_name=$(awk -v name="${name}" -F"," '$1==name {print $2}' $phylip_names)
// 		echo "current neighbor phylip name : ${phylip_name=}"

// 		# If there is an ortholog,
// 		if grep -q "^${orf}${tab}" ortho_pairs/${focal}_vs_${name}_orthologs.tsv
// 		then
// 			# Get it
// 			seq_to_add=$(grep "^${orf}${tab}" ortho_pairs/${focal}_vs_${name}_orthologs.tsv | awk -F"\t" '{print $2}')
// 			echo "ortholog found : ${seq_to_add}"

// 			# Add the neighbor phylip_name as header
// 			echo ">${phylip_name=}" >> $fa

// 			# Find its sequence and add it to the alignment fasta
// 			#faOneRecord ${B_fna} $seq_to_add | tail -n +2 >> $fa
// 			getonlyseq.sh ${B_fna} $seq_to_add >> $fa
// 		fi
// 	done

// 	seq_nb=$(grep -c ">" $fa)

// 	# Then, if there is more than one entry in the alignment fasta,
// 	if [ $seq_nb -gt 1 ]
// 	then
// 		# Remove terminal stop codons
// 		linearizefasta.sh $fa > ${fa}.tmp
// 		sed -E "/>/! s/...$//" ${fa}.tmp > $fa
// 		rm ${fa}.tmp

// 		# Align both nucleotides and AAs with MAFFT (see '-p F')
// 		translatorX.pl -i $fa -p F -o $orf

// 		# Get codons aligned with respect to amino acids 
// 		pal2nal.pl ${orf}.aa_ali.fasta ${orf}.nt_ali.fasta -output fasta > ${orf}_pal2nal.aln

// 		# If the below line is removed, you must linearize before format conversion (linearizefasta.sh ${orf}_pal2nal.aln > ${orf}_pal2nal.aln.oneline)
// 		remove_mostly_gapped_positions.py ${orf}_pal2nal.aln ${orf}_clean.aln

// 		# Turn the FASTA alignment to PHYLIP format (input_file, number of sequences, number of columns)
// 		FASTAtoPHYL.pl ${orf}_clean.aln $seq_nb $(awk 'FNR==2 {print length}' ${orf}_clean.aln)
// 		mv ${orf}_clean.aln.phy ${orf}.aln

// 		# Generate the phylogenetic tree with only the species that have a match
//    		get_subtree.py -tree $tree -names <(grep ">" $fa | sed "s/>//") -out ${orf}.nwk
// 	fi
// 	'''

// }






