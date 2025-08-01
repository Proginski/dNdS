# dN/dS Analysis Pipeline

A Nextflow workflow for comprehensive evolutionary analysis using CodeML to detect selection pressure through dN/dS ratio calculations and likelihood ratio tests.

## Overview

This NEXTFLOW workflow performs custom dN/dS analysis by building codon-aware alignments from orthologous sequences and testing multiple evolutionary models using CodeML. It implements branch-specific and subtree-specific models inspired by Zhang et al. 2019 to detect positive selection across phylogenetic trees.

## Input Parameters

- **outdir** = name of the results directory
- **A** = focal nucleotide FASTA file with the query sequences
- **B** = subject nucleotides FASTA file with the orthologs
- **includeA** = boolean flag (default: false, unnecessary)
- **A_seqlist** = a subset list of the focal FASTA queries
- **ortho** = a directory with 2-column TSV files containing pairs of orthologs (should be named "{A}\_vs\_{anyB}_orthologs.tsv")
- **tree** = a NEWICK tree of the genomes
- **ctl_file** = should be "${projectDir}/ctl/codeml_one_omega_model.ctl" (others not tested)
- **style** = should be 'branch_models' (others not tested)
- **fasta_batch_size** = batch size for FASTA processing (default: 10000)
- **alignment_batch_size** = batch size for alignment processing (default: 100)
- **codeml_batch_size** = batch size for CodeML analysis (default: 10)

## Usage

```
nextflow run <repository>/dNdS -profile <SINGULARITY-APPTAINER-DOCKER> --outdir <OUTDIR> --A <FOCAL_FASTA> --B <SUBJECT_FASTA> --ortho <ORTHO_DIR> --tree <NEWICK_TREE>
```

## Container Requirements

The workflow is expected to be run with a container manager (Singularity, Apptainer, or Docker).

## Workflow Strategy

### Alignment Building

The pipeline builds codon alignments using the following strategy:

1. **Dual alignment with MAFFT**: Align both nucleotides and amino acids using translatorX.pl with '-p F' option
2. **Codon alignment**: Generate codon-aligned sequences using pal2nal.pl respecting amino acid alignment
3. **Gap removal**: Remove mostly gapped positions using remove_mostly_gapped_positions.py
4. **Format conversion**: Convert FASTA alignment to PHYLIP format for CodeML compatibility

```bash
# Align both nucleotides and AAs with MAFFT
translatorX.pl -i $ortho_fna -p F -o $orf

# Get codons aligned with respect to amino acids 
pal2nal.pl ${orf}.aa_ali.fasta ${orf}.nt_ali.fasta -output fasta > ${orf}_pal2nal.aln

# Remove mostly gapped positions
remove_mostly_gapped_positions.py ${orf}_pal2nal.aln ${orf}_clean.aln

# Convert to PHYLIP format for CodeML
FASTAtoPHYL.pl ${orf}_clean.aln $(grep -c ">" ${orf}_clean.aln) $(awk 'FNR==2 {print length}' ${orf}_clean.aln)
```

### Evolutionary Model Testing

For each alignment, the workflow tests multiple evolutionary models:

#### Basic Models
- **1omega**: Single omega (dN/dS) across entire tree
- **1omega_fixed1**: Single omega fixed to 1

#### Branch-Specific Models (Zhang et al. 2019 approach)
- **2omega**: Two-omega model with foreground branch having different omega
- **2omega_fixed1**: Two-omega model with foreground omega fixed to 1

#### Advanced Models
- **free_omega**: Free omega model allowing different omega for each branch

## Statistical Analysis

### Model Comparison
- **Likelihood ratio tests**: Compare nested models to detect significant differences
- **AIC/BIC criteria**: Evaluate model fit and complexity
- **Parameter estimation**: Extract dN/dS ratios and confidence intervals

### Selection Detection
- **P-value thresholds**: Default significance level of 0.05
- **Bonferroni correction**: Multiple testing correction available
- **Branch-specific selection**: Identify specific lineages under selection

## Output Files

### Analysis Results
- **\*_models.tsv_bonferroni_corrected**
- **under_selection_bonferroni_corrected.txt**

## Dependencies

- **Nextflow** >= 23.10.0
- **Singularity OR Apptainer OR Docker**

## Methodology Reference

The branch-specific and subtree-specific analysis approach is inspired by:
Zhang et al. 2019. "The evolution of gene expression in cis and trans." Nature Ecology & Evolution. https://doi.org/10.1038/s41559-019-0822-5

## Citation

If you use this pipeline in your research, please cite the associated publication and the Zhang et al. 2019 methodology reference.