process NAMES_TO_PHYLIP {

	label 'local_job'

	input:
		path names

	output:
		path "${names}_phylip_names"

	"""
	awk '
        {
            short=substr(\$0,1,5)NR; lim_short=substr(short,1,10)
            final=gensub(/_/,"","g",lim_short)
            print \$0","final
        }
        ' $names > ${names}_phylip_names
	"""
}