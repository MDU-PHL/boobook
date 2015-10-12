# `boobook`: a tool to count RNA-seq reads mapping to genomic features
A Python tool that transforms RNA-seq data into an input file suitable for
[Degust](http://vicbioinformatics.com/degust/index.html).

# Etymology

A southern boobook (*Ninox boobook*) is an owl from Australia, NZ, and New Guinea.
It eats loads of insects (more than any other owl in Australia), and will
regurgitate the remains as a pellet. These pellets are studied carefully by
biologists to figure out what the owls eat. In the same spirit, `boobook` will
eat up all your RNA-seq reads, and produce pellets (tables) of collated counts
of reads mapping to chosen genomic features. Biologists can then study the pellets
to understand if there are differences in gene expression across treatments.

# Depends
* `BWA` (>= 0.7.12-r1039)

* `samtools` (>= 1.2 (using htslib 1.2.1))

* `HTseq` (>=version 0.6.0)

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

    * a comma-delimited file with features along rows --- the file can be
     directly loaded into Digust

stats.text

    * some QC data

# TODO

 * Add an argument parser

 * Fix up the input tab file

    - Needs to include a column to let the program know to update count

 * Need to allow the count of more than one feature --- right now only CDS
