# `boobook`: a Python tool to count RNA-seq reads mapping to genomic features
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

<!--
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
-->

# Outputs

## Folder structure

When creating a new `boobook` analysis, `boobook` will create the following
folder structure:

        project_home/
        +----------- .reads_dump
        +----------- input.tab
        +----------- sample_counts/
        |             +------------ read_set1/
        |             |             +------- read_set1.fq.gz --> symlink
        |             |             +------- alignment.bam
        |             +------------- read_set2/
        |             |
        |             +------------- read_setn/
        +------------ reference/
        |             +-------- reference.fa
        |             +-------- reference.gff
        +------------ pellets/
                      +------- boobook_pellet.csv

* `.reads_dump`: a special file used to keep track of what has already been done
to avoid repeating analysis. **Don't touch!**.
* `input.tab`: the file specifying input reads. If you new read sets have become
available, just **add** new rows to this file, and run `boobook` again.
* `sample_counts`: a folder holding information about where the program can
find the reads, and where the `alignment.bam` is kept.
* `reference`: a folder that holds processed `reference.fa` and `reference.gff`
required by `boobook` to perform the alignment and counting.
* `pellets`: where `boobook` puts its pellets. The file `boobook_pellet.csv` is
described below.

## Output files

The outputs for `boobook` are kept in a folder called `pellets`. In the folder,
you will find a single file:

`boobook_pellet.csv`:

* Each row contains counts of reads mapping to a single feature across all
 treatments and replicates

* Along with the count, each row has metadata that can help identify the
feature. At the moment these include:
    - locus_tag
    - Product description (**this field is quoted**)
    - GeneID
    - GI
    - Gene (if present)


* The file contains a header that describes the content of the column
* a comma-delimited file with features along rows --- the file can be directly
loaded into Digust

<!--
stats.text

    * some QC data
-->
# TODO

 * Add an argument parser

 * Fix up the input tab file

    - Needs to include a column to let the program know to update count

 * Need to allow the count of more than one feature --- right now only CDS

 * Update documentation

# History

* 2015-10-14: Pellets can be saved, and rows contain metadata to identify the
feature.
