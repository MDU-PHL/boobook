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

# Inputs

## Positional arguments
* REF --- A Genbank file: the defaults assume that it contains a locus_tag

* INFILE:

    * A comma-delimited or tab-delimited file with the following columns
    (it assumes no header column):
      1. SampleID --- a string **Must uniquely identify each row**
      2. TreatmentID -- a string **Must not include spaces**
      3. ReplicateID --- a string or an integer
      4. Fastq1 --- /path/to/fastq1.fq.gz
      5. Fastq2 --- /path/to/fastq2.fq.gz (this column is OPTIONAL)

In this file, the `TreatmentID` + `ReplicateID` will form the column header. This
way, `Degust` can group samples by **treatment**. It is important that the
treatment IDs be **unique**. Thus, if you have different samples (e.g., samples
from different strains) that have undergone the same experimental treatment, then that
information should be included in the treatment ID. Here is an example file,
where two strains (**A**, and **B**), have been treated to the same two experimental treatments
(**T1** and **T2**), with two replicates of each experimental treatment per strain (this makes for a total of 2 x 2 x 2 = 8 rows):

|SampleID |TreatmentID|ReplicateID|Fastq1|
|---------|--------|----|------------------|
|Sample1  |A_T1    |1   |/path/to/A_T1_1.fq|
|Sample2  |A_T1    |2   |/path/to/A_T1_2.fq|
|Sample3  |A_T2    |1   |/path/to/A_T2_1.fq|
|Sample4  |A_T2    |2   |/path/to/A_T2_2.fq|
|Sample5  |B_T1    |1   |/path/to/B_T1_1.fq|
|Sample6  |B_T1    |2   |/path/to/B_T1_2.fq|
|Sample7  |B_T2    |1   |/path/to/B_T2_1.fq|
|Sample8  |B_T2    |2   |/path/to/B_T2_2.fq|


## Optional arguments
--work_dir: the folder where the boobook analysis (default: '.'). Will try to
            create the folder if it does not exist

--features: A comma separated list of features to count (e.g., CDS;snRNA). must
            match features in the Genbank file (default: CDS).

--qualifier: A unique identifier for individual elements of the genomic
             annotation --- we highly recommend using **locus_tag**
             (default: locus_tag)

--bwa_threads: [**BWA option**] Number of threads used by BWA during mapping
               (default: 16)

--sam_threads: [**SAMtools options**] Number of threads used by SAMtools when
               generating BAM file and sorting (default: 8)

--hts_stranded: [**HTSeq option**] Strandedness of RNAseq data. Possible options:
                'yes', 'no', and 'reversed' (default: 'reversed')

--hts_overlap: [**HTSeq option**] Tells **HTSeq** how to count overlapping features.
                The options are: 'union', 'intersection-strict', 'intersection-nonempty'.
                The user is directed to the HTSeq website to find out [more](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html).
                (default: 'reversed')
--add: Re-run analysis with additional samples in the infile or with additional
        features.

--re_align: Force re-alignment of all samples in the project folder. Will force
            re-count.

--re_count: Force re-count of all samples in the project folder.

--change_ref: Re-run analysis with a different reference. Will force re-alignment
              and re-count.

# Examples

Simple run with BWA alignment and HTSeq counts run in the current directory:

        boobook.py myref.gbk myinfile.tab

Outputting `boobook` analysis to `my_workdir`:

        boobook.py --work_dir my_workdir myref.gbk myinfile.tab

Assuming RNAseq data is not stranded:

        boobook.py --hts_stranded no myref.gbk myinfile.tab

Want to count both CDS and  snRNA:

        boobook.py --features CDS;snRNA myref.gbk myinfile.tab

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
* The file can be uploaded directly in to
[Degust](http://vicbioinformatics.com/degust/index.html) for analysis

<!--
stats.text

    * some QC data
-->
# TODO
    - Update the Counter class to make things a little more modular
    - Add support for sailfish and kallisto
    - Ensure column names for treatment are unique

# History
* 2015-11-10
    - Output contains a 'Length' column for each feature
    - The `--version` flag now works
    - Version name in GFF file is now consistent with the version option
    - Updated information on the input file
    - Added check to ensure that `TreatmentID` + `ReplicateID` leads to unique
        column name combinations to avoid an error in `Degust`. If an error is
        found, then `boobook` raises an error with information on why, and quits.
        The check is done early on, to avoid unnecessary waiting by the user.

* 2015-11-09:
    - Added argument parser
    - Updated documentation to include information on infile and arguments/options
    - Can now count any arbitrary number of feature types


* 2015-10-14: Pellets can be saved, and rows contain metadata to identify the
feature.
