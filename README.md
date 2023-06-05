Calculates the impact of a mutation on a Coding sequence using RNASeq data

Generates a .tsv file that separates the mutations into three types; Synonymous, Non-Synonymous and Non-Coding/Nonsense

Requires:
VCF for the mutations
GFF of the genome
Fasta for the genome

Basic Syntax:
mtclassifier.py --vcf <VCF file> --gff <GFF file> --Fasta <FASTA file> --output <Output file name> (optional -q <Quality score>)  
