# Rosa-hybrida-Samantha
Haplotype-resolved genome assembly and resequencing provide insights into the origin and domestication of modern rose.

In this study, the haplotype-resolved tetraploid genome of modern rose was assembled using the Pore-C data.
Due to the complex ancestral origins of modern roses, which result in a segmental allotetraploid genome characteristic among the chromosomes, correctly phasing between homologous chromosomes is challenging using only Hi-C technology. The Pore-C technology was utilized for precise orientation between homologous chromosomes to reduce assembly errors.

BamModifier.py: Used to split BAM files and generate bin.bam files.

calculate_sum.py: Used for gene retention analysis, calculating the proportion of gene retention within each window.

compare_genotypes.py: Used for analyzing VCF file of several samples, including calculating the proportion of each type of inconsistent sites among the samples.

extract_genes.py: Used for extracting genes within a specific interval.

vcf_genotype_stats.py: Used for counting the number of each type of site in VCF file for each sample.
