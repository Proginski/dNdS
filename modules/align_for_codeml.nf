process ALIGN_FOR_CODEML {

	label 'dnds'

	publishDir "${params.outdir}/alignments/", pattern: "*_for_codeml.aln"

	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
	maxRetries 1

	input:
		path orthologs_fnas

	output:
		path('*_for_codeml.aln'), optional : true

	shell:
	'''
	for ortho_fna in !{orthologs_fnas.join(' ')}
	do

		orf=$(basename $ortho_fna _n_orthologs.fna)

		# Align both nucleotides and AAs with MAFFT (see '-p F')
		translatorX.pl -i $ortho_fna -p F -o $orf

		# Get codons aligned with respect to amino acids 
		pal2nal.pl ${orf}.aa_ali.fasta ${orf}.nt_ali.fasta -output fasta > ${orf}_pal2nal.aln

		# If the below line is removed, you must linearize before format conversion (linearizefasta.sh ${orf}_pal2nal.aln > ${orf}_pal2nal.aln.oneline)
		remove_mostly_gapped_positions.py ${orf}_pal2nal.aln ${orf}_clean.aln

		# Turn the FASTA alignment to PHYLIP format (input_file, number of sequences, number of columns)
		FASTAtoPHYL.pl ${orf}_clean.aln $(grep -c ">" ${orf}_clean.aln) $(awk 'FNR==2 {print length}' ${orf}_clean.aln)
		mv ${orf}_clean.aln.phy ${orf}_for_codeml.aln

	done
	'''
}