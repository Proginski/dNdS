process GET_SEQUENCES {

	input:
		path fasta
		path user_list

	output:
		path "sequences_list.txt"

	shell = ['/bin/bash','-u']
    """
    grep ">" $fasta | awk '{print \$1}' | sed 's/>//' > sequences_list.txt
	if [ -s $user_list ]
	then
		grep -vxf sequences_list.txt $user_list > unrecognized_sequences.txt
		if [ -s unrecognized_sequences.txt ]
		then
			echo "The following sequences are not present in the A FASTA file (see 'unrecognized_sequences.txt'):"
			cat unrecognized_sequences.txt
			exit 1
		else
			echo "All sequences in the list are present in the A FASTA file."
			mv $user_list sequences_list.txt
		fi
	fi
	"""
}