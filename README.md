# rna_count.py
A Python tool that transforms RNA-seq data into an input file suitable for
Degust

# Inputs

## Positional arguments
<project_file>.<csv|tsv>:

    * A comma-delimited or tab-delimited file with the following columns:
      1. SampleID --- a string
      2. ReplicateID --- a string or an integer
      3. SeqType --- either SE or PE for single and paired-end, respectively
      4. Fastq1 --- /path/to/fastq1.fq.gz
      5. Fastq2 --- /path/to/fastq2.fq.gz (ignored if SeqType = SE)

<reference>.<gbk|gff>

    * A Genbank or GFF file with a FASTA file

## Optional arguments
--threads: number of threads to use when mapping (default: 16)
--outdir: /path/to/outdir (default: .)
--stranded: reverse, forward, both (default: reverse)
--feature: which feature to count mapping reads from (default: CDS)
            must be in the GBK or GFF file

# Outputs

counts.csv

    * a comma-delimited file with features along rows

stats.text

    * some QC data

# Todo

  * Separate genbank parsing from saving gff file in class GenBank
  * Create class Project 
