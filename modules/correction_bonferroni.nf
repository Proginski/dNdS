process CORRECTION_BONFERRONI {

	publishDir "${params.outdir}/codeml_out/models_bonferroni_corrected", pattern: "*_models.tsv_bonferroni_corrected"
	publishDir "${params.outdir}", pattern: "under_selection_bonferroni_corrected.txt"

	input:
		path models

	output:
		path "*_models.tsv_bonferroni_corrected"
		path "under_selection_bonferroni_corrected.txt"

	script:
	"""
	# Count the models
	echo "${models.size()}"

	for models_tsv in ${models}
	do
		# Only process if there is at least one "true" in the file
		if grep -q "true\$" \${models_tsv}
		then
			# Multiply the p-value by the number of models
			awk -v total_models=${models.size()} '
				BEGIN{FS=OFS="\t"}
				FNR>1 {
					if(\$7==""){ \$7="NA";\$8="NA" }
					if(\$8!="NA"){ \$8=\$8*total_models }
				}
				{ print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8 }
			' \${models_tsv} > \${models_tsv}_bonferroni_corrected

			detect_selection_with_compared_models.py \${models_tsv}_bonferroni_corrected --pvalue 0.05
		else
			# If no "true" is found, just copy the file without changes
			cp \${models_tsv} \${models_tsv}_bonferroni_corrected
		fi
	done

	touch under_selection_bonferroni_corrected.tsv
	grep -c "true\$" *_models.tsv_bonferroni_corrected | awk -F":" '\$2>0 { print \$1 }' | sed 's/_models.tsv_bonferroni_corrected//g' > under_selection_bonferroni_corrected.txt
	"""
}