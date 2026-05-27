Scripts to functionally annotate SNPs obtained from 3RAD data from _Anolis cybotes_ and _Anolis shrevei_. 

# Input files:
1. Genome assembly .fasta file (genome.fa)
2. Genome annotation .gtf file (AnoSag2.1.gtf)
3. 3RAD assembly from ipyRAD in .fasta format (assembly.fa)
4. 3RAD variants file in .tsv format (variants_noprefix.tsv)

# The general pipeline goes as follows:
1. Annotate CDS in the genome (gtf_with_cds.py).\
   This is done because the original .gtf only includes exons, and start and stop codons. We add the CDS features to the annotation.\
   Output = AnoSag2.1_cds.gtf
2. Blast assembly against reference genome.\
   This is run in blast and we generate a .xml output
3. Get SNP coordinates in the reference genome (snp_genome_coords.py).\
   This script uses as input the blast.xml output from step 2 and the variants.tsv file to generate a file with the new coordinates for SNPs.\
   Output = snp_genome_coords.tsv
4. Annotate the SNPs (snps_in_cds.py, snps_in_3UTR.py, snps_in_5UTR.py, snps_in_intron.py).\
   Here we generate a file for each feature. We use the "AnoSag2.1_cds.gtf" and "snp_genome_coords.tsv" as inputs.
5. Get the effect of substitutions from SNPs in CDS (snp_cds_effects.py).\
   With this script we identify if the substitutions are synonymous or non-synonymous. The inputs are "AnoSag2.1_cds.gtf", "genome.fa", and "snp_genome_coords.tsv".
