#!/usr/bin/env python
'''
A Python tool to count reads mapped to genomic features
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.GenBank.Scanner import GenBankScanner
#from Bio.Sequencing.Applications import _samtools as samtools
from Bio.Application import _Option, AbstractCommandline, _Switch, _StaticArgument, _Argument
from operator import itemgetter
import HTSeq
import itertools
import os
import pickle
import subprocess
import re
import shutil
import sys
import time
import traceback
import warnings
import click

## default parameters
version = "boobook:0.1"

# to ensure click gives us the option for -h and --help
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

## modified count script from HTSeq

class UnknownChrom( Exception ):
   pass

def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def count_reads_in_features( sam_filename, gff_filename, samtype, order, stranded, \
      overlap_mode, feature_type, id_attribute, quiet, minaqual, samout ):

   def write_to_samout( r, assignment ):
      if samoutfile is None:
         return
      if not pe_mode:
         r = (r,)
      for read in r:
         if read is not None:
            samoutfile.write( read.original_sam_line.rstrip() +
               "\tXF:Z:" + assignment + "\n" )


   if samout != "":
      samoutfile = open( samout, "w" )
   else:
      samoutfile = None

   features = HTSeq.GenomicArrayOfSets( "auto", stranded != "no" )
   counts = {}

   # Try to open samfile to fail early in case it is not there
   if sam_filename != "-":
      open( sam_filename ).close()

   gff = HTSeq.GFF_Reader( gff_filename )
   i = 0
   try:
      for f in gff:
         if f.type == feature_type:
            try:
               feature_id = f.attr[ id_attribute ]
            except KeyError:
               raise ValueError, ( "Feature %s does not contain a '%s' attribute" %
                  ( f.name, id_attribute ) )
            if stranded != "no" and f.iv.strand == ".":
               raise ValueError, ( "Feature %s at %s does not have strand information but you are "
                  "running htseq-count in stranded mode. Use '--stranded=no'." %
                  ( f.name, f.iv ) )
            features[ f.iv ] += feature_id
            counts[ f.attr[ id_attribute ] ] = 0
         i += 1
         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d GFF lines processed.\n" % i )
   except:
      sys.stderr.write( "Error occured when processing GFF file (%s):\n" % gff.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d GFF lines processed.\n" % i )

   if len( counts ) == 0:
      sys.stderr.write( "Warning: No features of type '%s' found.\n" % feature_type )

   if samtype == "sam":
      SAM_or_BAM_Reader = HTSeq.SAM_Reader
   elif samtype == "bam":
      SAM_or_BAM_Reader = HTSeq.BAM_Reader
   else:
      raise ValueError, "Unknown input format %s specified." % samtype

   try:
      if sam_filename != "-":
         read_seq_file = SAM_or_BAM_Reader( sam_filename )
         read_seq = read_seq_file
         first_read = iter(read_seq).next()
      else:
         read_seq_file = SAM_or_BAM_Reader( sys.stdin )
         read_seq_iter = iter( read_seq_file )
         first_read = read_seq_iter.next()
         read_seq = itertools.chain( [ first_read ], read_seq_iter )
      pe_mode = first_read.paired_end
   except:
      sys.stderr.write( "Error occured when reading beginning of SAM/BAM file.\n" )
      raise

   try:
      if pe_mode:
         if order == "name":
            read_seq = HTSeq.pair_SAM_alignments( read_seq )
         elif order == "pos":
            read_seq = HTSeq.pair_SAM_alignments_with_buffer( read_seq )
         else:
            raise ValueError, "Illegal order specified."
      empty = 0
      ambiguous = 0
      notaligned = 0
      lowqual = 0
      nonunique = 0
      i = 0
      for r in read_seq:
         if i > 0 and i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d SAM alignment record%s processed.\n" % ( i, "s" if not pe_mode else " pairs" ) )

         i += 1
         if not pe_mode:
            if not r.aligned:
               notaligned += 1
               write_to_samout( r, "__not_aligned" )
               continue
            try:
               if r.optional_field( "NH" ) > 1:
                  nonunique += 1
                  write_to_samout( r, "__alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if r.aQual < minaqual:
               lowqual += 1
               write_to_samout( r, "__too_low_aQual" )
               continue
            if stranded != "reverse":
               iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
            else:
               iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" and co.size > 0 )
         else:
            if r[0] is not None and r[0].aligned:
               if stranded != "reverse":
                  iv_seq = ( co.ref_iv for co in r[0].cigar if co.type == "M" and co.size > 0 )
               else:
                  iv_seq = ( invert_strand( co.ref_iv ) for co in r[0].cigar if co.type == "M" and co.size > 0 )
            else:
               iv_seq = tuple()
            if r[1] is not None and r[1].aligned:
               if stranded != "reverse":
                  iv_seq = itertools.chain( iv_seq,
                     ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
               else:
                  iv_seq = itertools.chain( iv_seq,
                     ( co.ref_iv for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
            else:
               if ( r[0] is None ) or not ( r[0].aligned ):
                  write_to_samout( r, "__not_aligned" )
                  notaligned += 1
                  continue
            try:
               if ( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 ):
                  nonunique += 1
                  write_to_samout( r, "__alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               write_to_samout( r, "__too_low_aQual" )
               continue

         try:
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     fs = fs.union( fs2 )
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )
            if fs is None or len( fs ) == 0:
               write_to_samout( r, "__no_feature" )
               empty += 1
            elif len( fs ) > 1:
               write_to_samout( r, "__ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:
               write_to_samout( r, list(fs)[0] )
               counts[ list(fs)[0] ] += 1
         except UnknownChrom:
            write_to_samout( r, "__no_feature" )
            empty += 1

   except:
      sys.stderr.write( "Error occured when processing SAM input (%s):\n" % read_seq_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d SAM %s processed.\n" % ( i, "alignments " if not pe_mode else "alignment pairs" ) )

   if samoutfile is not None:
      samoutfile.close()

   # for fn in sorted( counts.keys() ):
   #    print "%s\t%d" % ( fn, counts[fn] )
   print "__no_feature\t%d" % empty
   print "__ambiguous\t%d" % ambiguous
   print "__too_low_aQual\t%d" % lowqual
   print "__not_aligned\t%d" % notaligned
   print "__alignment_not_unique\t%d" % nonunique
   return counts

class SamtoolsViewCommandline(AbstractCommandline):
    """Command line wrapper for samtools view.

    Extract/print all or sub alignments in SAM or BAM format, equivalent to::

        $ samtools view [-bchuHS] [-t in.refList] [-o output] [-f reqFlag]
                        [-F skipFlag] [-q minMapQ] [-l library] [-r readGroup]
                        [-R rgFile] <in.bam>|<in.sam> [region1 [...]]

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Example
    -------

    >>> from Bio.Sequencing.Applications import SamtoolsViewCommandline
    >>> input_file = "/path/to/sam_or_bam_file"
    >>> samtools_view_cmd = SamtoolsViewCommandline(input_file=input_file)
    >>> print(samtools_view_cmd)
    samtools view /path/to/sam_or_bam_file

    """
    def __init__(self, cmd="samtools", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("view"),
            _Option(["-@","threads"],
                "Number of BAM compression threads [0]",
                equate = False),
            _Option(["-T","ref"],
                "Reference sequence FASTA FILE [null]",
                equate = False),
            _Switch(["-b", "b"], "Output in the BAM format"),
            _Switch(["-c", "c"],
                    """Instead of printing the alignments, only count them and
                    print the total number.

                    All filter options, such as '-f', '-F' and '-q',
                    are taken into account"""),
            _Switch(["-h", "h"], "Include the header in the output"),
            _Switch(["-u", "u"],
                    """Output uncompressed BAM.

                    This option saves time spent on compression/decompression
                    and is thus preferred when the output is piped to
                    another samtools command"""),
            _Switch(["-H", "H"], "Output the header only"),
            _Switch(["-S", "S"],
                    """Input is in SAM.
                    If @SQ header lines are absent,
                    the '-t' option is required."""),
            _Option(["-t", "t"],
                    """This file is TAB-delimited.
                    Each line must contain the reference name and the
                    length of the reference, one line for each
                    distinct reference; additional fields are ignored.

                    This file also defines the order of the reference
                    sequences in sorting.
                    If you run   'samtools faidx <ref.fa>',
                    the resultant index file <ref.fa>.fai can be used
                    as this <in.ref_list> file.""",
                    filename=True, equate=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-o", "o"], "Output file",
                    filename=True, equate=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-f", "f"],
                    """Only output alignments with all bits in
                    INT present in the FLAG field""",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-F", "F"],
                    "Skip alignments with bits present in INT",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-q", "q"],
                    "Skip alignments with MAPQ smaller than INT",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-r", "r"],
                    "Only output reads in read group STR",
                    equate=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-R", "R"],
                    "Output reads in read groups listed in FILE",
                    filename=True, equate=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-l", "l"],
                    "Only output reads in library STR",
                    equate=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Switch(["-1", "fast_bam"],
                    "Use zlib compression level 1 to compress the output"),
            _Argument(["input", "input_file"],
                      "Input File Name", filename=True, is_required=True),
            _Argument(["region"], "Region", is_required=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class SamtoolsSortCommandline(AbstractCommandline):
    """Command line wrapper for samtools sort.

    Concatenate BAMs, equivalent to::

    $ samtools sort [-no] [-m maxMem] <in.bam> <out.prefix>

    See http://samtools.sourceforge.net/samtools.shtml for more details

    Example
    -------

    >>> from Bio.Sequencing.Applications import SamtoolsSortCommandline
    >>> input_bam = "/path/to/input_bam"
    >>> out_prefix = "/path/to/out_prefix"
    >>> samtools_sort_cmd = SamtoolsSortCommandline(\
                                                    input_bam=input_bam,\
                                                    out_prefix=out_prefix)
    >>> print(samtools_sort_cmd)
    samtools sort /path/to/input_bam /path/to/out_prefix

    """
    def __init__(self, cmd="samtools", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _StaticArgument("sort"),
            _Option(["-@", "threads"],
                "Set number of sorting and compression threads [1]",
                equate = False),
            _Switch(["-o", "o"], """Output the final alignment
                                    to the standard output"""),
            _Switch(["-n", "n"], """Sort by read names rather
                                    than by chromosomal coordinates"""),
            _Option(["-m", "m"], "Approximately the maximum required memory",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Argument(["input_bam"], "Input BAM file",
                      filename=True, is_required=True),
            _Argument(["out_prefix"], "Output prefix",
                      filename=True, is_required=True)
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class BWAindex(AbstractCommandline):
    def __init__(self, cmd = 'bwa', **kwargs):
        self.parameters = [
            _StaticArgument("index"),
            _Option(["-a", "algorithm"],
                "BWT construction algorithm: bwtsw or is [auto]",
                equate = False),
            _Option(["-p", "prefix"],
                "prefix of the index [same as fasta name]",
                equate = False),
            _Option(["-b", "block_size"],
                "block size for the bwtsw algorithm \
                (effective with -a bwtsw) [10000000]",
                equate = False),
            _Switch(['-6', "index_64"],
                "index files named as <in.fasta>.64.* instead of <in.fasta>.*"),
            _Argument(["in_fasta"],
                    "Input FASTA file",
                    filename = True,
                    is_required = True)
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        return

class BWAmem(AbstractCommandline):
    def __init__(self, cmd = 'bwa', **kwargs):
        self.parameters = [
            _StaticArgument("mem"),
        #Algorithm options
            _Option(["-t", "threads"],
                "number of threads [default: 1]",
                equate = False),
            _Option(["-k", "seed_len"],
                "minimum seed length [default: 19]",
                equate = False),
            _Option(["-w", "band_width"],
                "band width for banded alignment [default: 100]",
                equate = False),
            _Option(["-d", "off_diag"],
                "off-diagonal X-dropoff [100]",
                equate = False),
            _Option(["-r", "intern_seeds"],
                "look for internal seeds inside a seed longer \
                than {-k} * FLOAT [1.5]",
                equate = False),
            _Option(["-y", "seed_occur"],
                "seed occurrence for the 3rd round seeding [20]",
                equate = False),
            _Option(["-c", "skip_seeds"],
                "skip seeds with more than INT occurrences [500]",
                equate = False),
            _Option(["-D", "drop_chains"],
                "drop chains shorter than FLOAT fraction of the longest \
                overlapping chain [0.50]",
                equate = False),
            _Option(["-W", "discard_chain"],
                "discard a chain if seeded bases shorter than INT [0]",
                equate = False),
            _Option(["-m", "max_rounds"],
                "perform at most INT rounds of mate rescues for each read [50]",
                equate = False),
            _Switch(["-S", "skip_mate"],
                "skip mate rescue"),
            _Switch(['-P', "skip_pairing"],
                "skip pairing; mate rescue performed unless -S also in use"),
            _Switch(["-e", "discard_matches"],
                "discard full-length exact matches"),
        #Scoring options
            _Option(["-A", "match_score"],
                "score for a sequence match, which scales \
                options -TdBOELU unless overridden [1]",
                equate = False),
            _Option(["-B", "mis_penalty"],
                "penalty for a mismatch [4]",
                equate = False),
            _Option(["-O", "gap_open_penal"],
                "gap open penalties for deletions and insertions [6,6]",
                equate = False),
            _Option(["-E", "gap_exten_penal"],
                "gap extension penalty; a gap of size k cost \
                '{-O} + {-E}*k' [1,1]",
                equate = False),
            _Option(["-L", "clip_penal"],
                "penalty for 5'- and 3'-end clipping [5,5]",
                equate = False),
            _Option(["-U", "unpair_penal"],
                "penalty for an unpaired read pair [17]",
                equate = False),
            _Option(["-x", "read_type"],
                "read type. \
                Setting -x changes multiple parameters unless overriden [null] \
                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref) \
                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref) \
                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)",
                equate = False),
        #Input/output options
            _Switch(["-p", "smart_pair"],
                "smart pairing (ignoring in2.fq)"),
            _Option(["-R", "read_group"],
                "read group header line such as '@RG\tID:foo\tSM:bar' [null]",
                equate = False),
            _Option(["-H", "insert_header"],
                "insert STR to header if it starts with @; \
                or insert lines in FILE [null]",
                equate = False),
            _Switch(["-j", "treat_alt"],
                "treat ALT contigs as part of the primary assembly \
                (i.e. ignore <idxbase>.alt file)"),
            _Option(["-v", "verbose_level"],
                "verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]",
                equate = False),
            _Option(["-T","min_score"],
                "minimum score to output [30]",
                equate = False),
            _Option(["-h", "min_hit_score"],
                "if there are <INT hits with score >80% of the max score, \
                output all in XA [5,200]",
                equate = False),
            _Switch(["-a", "output_all"],
                "output all alignments for SE or unpaired PE"),
            _Switch(["-C", "append_comment"],
                "append FASTA/FASTQ comment to SAM output"),
            _Switch(["-V", "output_ref"],
                "output the reference FASTA header in the XR tag"),
            _Switch(["-Y", "soft_clip"],
                "use soft clipping for supplementary alignments"),
            _Switch(["-M", "mark_splits"],
                "mark shorter split hits as secondary"),
            _Option(["-I", "insert_size"],
                "FLOAT[,FLOAT[,INT[,INT]]] \
                specify the mean, standard deviation (10% of the mean if absent), \
                max (4 sigma from the mean if absent) and min of the insert \
                size distribution. FR orientation only. [inferred]",
                equate = False),
            _Argument(["ref"],
                "Input FASTA reference",
                filename = True,
                is_required = True),
            _Argument(["in_fq1"],
                    "Input FASTQ file 1",
                    filename = True,
                    is_required = True),
            _Argument(["in_fq2"],
                    "Input FASTQ file 2",
                    filename = True,
                    is_required = False)
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        return

class GBk:
    '''
    A class that reads and stores GenBank
    '''
    def __init__(self, path, force):
        self.reference = os.path.join(path, "reference")
        if os.path.isdir(self.reference):
            if not force:
                raise RuntimeError('''
                    Analysis has been run already. Use --change_ref flag to \n
                            re-run using the same directory with a different
                            reference'''
                    )
        else:
            os.makedirs(self.reference)
        return
    def __parse_location(self, coords):
        '''
        PRIVATE function to parse feature location
        '''
        start = int(coords.start) + 1 # added one because it outputs in Python
                                      # count (i.e., starting at zero)
        end = int(coords.end)
        strand = '+' if coords.strand else '-'
        return (start, end, strand)
    def read_gb(self, infile):
        '''
        A method to load a GenBank file into memory, and create a dictionary
        for future use in outputting GFF files, and a Reference FASTA file.

        >>> data = GBk()
        >>> data.read_gb(infile = "test/Sa6008.gbk")

        '''
        fn = os.path.basename(infile)
        base = os.path.splitext(fn)[0]
        self.fasta = os.path.join(self.reference, base+".fa")
        self.gff = os.path.join(self.reference, base+".gff")
        f = open(infile, 'r')
        self.gb = SeqIO.read(f, 'genbank')
        print "Read file " + infile + " and found:"
        print "Sequence %s with %i features" % (self.gb.description, len(self.gb.features))
        f.close()
        return

    def filter_features(self, features, qualifier):
        '''
        A method to filter out the features that one wants to count number of
        mapping reads

        >>> data = GBk()
        >>> data.read_gb(infile = "test/Sa6008.gbk")
        >>> data.filter_features(features = ['CDS', 'tRNA'])

        '''
        #transfomr qualifier into ascii
        # this allows the assert below to pass if the string is in unicode
        # format, as it comes out of click
        qualifier = qualifier.encode('ascii', 'ignore')
        if (not isinstance(features, list)):
            features = [features]
        assert type(qualifier) is str, "qualifier is not a string: %r" % qualifier
        self.features = {}
        tmp_count = 0
        tmp_dict = {}
        count_found_features = 0
        count_found_features_with_qualifiers = 0
        count_per_feature = dict((f, 0) for f in features)
        for (index, feat) in enumerate(self.gb.features):
            tmp_count += 1
            try:
                tmp_dict[feat.type] += 1
            except:
                tmp_dict[feat.type] = 1
            if feat.type in features:
                count_found_features += 1
                if qualifier in feat.qualifiers:
                    count_found_features_with_qualifiers += 1
                    count_per_feature[feat.type] += 1
                    qual = feat.qualifiers[qualifier]
                    if len(qual) is not 1:
                        raise("Expected single entry for qualifier, assuming the \
                        first element is correct.")
                    qual = qual[0]
                    try:
                        self.features[qual][feat.type] = feat
                    except:
                        self.features[qual] = {}
                        self.features[qual][feat.type] = feat
        print "Total features found = " + str(count_found_features)
        print "Total features found with qualifier = " + str(count_found_features_with_qualifiers)
        print "Counts per feature type:"
        for f in features:
            print "\t" + f + " = " + str(count_per_feature[f])
        print "#####"
        for f in tmp_dict:
            print "\t" + f + " = " + str(tmp_dict[f])
        return
    def write_gff(self):
        '''
        Write a GFF output file that can be read by htseq-count
        '''
        global version
        f = open(self.gff, 'w')
        x = len(self.features.keys())
        count1 = 0
        count2 = 0
        countPlus = 0
        refname = self.gb.id
        source = 'RefSeq'
        feature = version
        score = "."
        frame = 0
        gff_recs = []
        for k in self.features:
            feature_d = self.features[k]
            feature_k = feature_d.keys()
            if 'CDS' in feature_k:
                tmp_rec = feature_d['CDS']
            else:
                tmp_rec = feature_d[feature_k[0]]
            if feature_k[0] == 'tRNA':
                print feature_d[feature_k[0]]
            ltag = tmp_rec.qualifiers['locus_tag'][0]
            coords = tmp_rec.location
            (start, end, strand) = self.__parse_location(coords)
            attrib = "ID="+ltag+";Type="+feature_k[0]
            gff_recs.append((refname, source, feature, start, end, score, strand, frame, attrib))
            n_keys = len(feature_k)
            if n_keys == 1:
                count1 += 1
            elif n_keys == 2:
                count2 += 1
            else:
                countPlus += 1
        sorted_gff_recs = sorted(gff_recs, key=itemgetter(3))
        out_str = ""
        out_str += "##gff-version 3\n"
        out_str += "##source-version boobook.py 0.01\n"
        out_str += "##data This file was mode on " + time.strftime("%d/%m/%Y") + "\n"
        out_str += "##description " + self.gb.description + "\n"
        out_str += "##sequence-region " + refname + " " + "1 " + str(len(self.gb.seq)) + "\n"
        for r in sorted_gff_recs:
            out_str+="\t".join(map(str,r))+'\n'
        f.write(out_str)
        f.close()
        return
    def write_fasta(self):
        '''
        This function outputs the sequence data in the GenBank file to a
        FASTA file called ref.fa
        '''
        print self.gff
        print self.fasta
        f = open(self.fasta, "w")
        SeqIO.write(self.gb, f, 'fasta')
        f.close()
        print "Reference file was successfully written."
        return

class GFF:
    def __init__(self):
        return

class Counter:
    def __init__(self, type = 'b+h'):
        self.type = type
    def count(self,)

class ReadData:
    def __init__(self, path, force):
        self.workdir = path
        self.sample_path = os.path.join(self.workdir, "sample_counts")
        self.ref_path = os.path.join(self.workdir, "reference")
        self.results_path = os.path.join(self.workdir, "pellets")
        self.__dict__["reads"] = {}
        if os.path.isdir(self.sample_path):
            if not force:
                raise RuntimeError('''
                    Analysis has been run already. Use --add flag to re-run using the same directory'''
                    )
            else:
                print "## Loading pickle..."
                old_fn = os.path.join(self.workdir,".reads_dump")
                print old_fn, os.path.isfile(old_fn)
                if os.path.isfile(old_fn):
                    print "Found pickle..."
                    s = self.Load(old_fn)
                    self.__dict__ = s.__dict__
                return
        os.makedirs(self.sample_path)
        # os.makedirs(self.ref_path)
        # os.makedirs(self.results_path)
        return
    @classmethod
    def Load(cls, file):
        f = open(file, "rb")
        s = pickle.load(f)
        f.close()
        return s
    def __iter__(self):
        return iter(self.__dict__["reads"])
    def __getitem__(self, key):
        return self.__dict__["reads"][key]
    def __setitem__(self, key, value):
        self.__dict__["reads"][key] = value
        return
    def read_input(self, infile):
        '''
        Reads an input tab delimited file with one sample per row. Currently,
        assumes single read data set per sample.

        The file should have the following fields:
        ID \t Treatment ID \t Replicate Number \t /Path/To/Read.fq<.gz>

        Features specifies which features to count. It should be a colon separated
        values (e.g., CDS,tRNA). A single feature is also possible (e.g., CDS).

        The program will assume the same reference for all the samples. This
        seems reasonable given the analysis
        '''
        self.input_path = os.path.join(self.workdir, infile)
        fi = open(self.input_path)
        for line in fi:
            (sample_id, treat, rep, path) = tuple(line.strip().split("\t"))
            if sample_id not in self.__dict__["reads"]:
                self.__dict__["reads"][sample_id] = {"treat":treat, \
                                        "rep":rep, \
                                        "path":path}
            else:
                print "Sample %s found..." % sample_id
        fi.close()
        return
    def create_subfolders(self):
        for r in self.__dict__["reads"]:
            print "Creating necessary files for %s" % r
            tmp_dir = os.path.join(self.sample_path, r)
            if not os.path.isdir(tmp_dir):
                os.mkdir(tmp_dir)
            if re.search("\.gz$", self.__dict__["reads"][r]["path"]):
                tmp_reads = os.path.join(tmp_dir, r + "_R.fq.gz")
            else:
                tmp_reads = os.path.join(tmp_dir, r + "_R.fq")
            if not os.path.islink(tmp_reads):
                os.symlink(self.__dict__["reads"][r]["path"], tmp_reads)
            self.__dict__["reads"][r]["link"] = tmp_reads
            self.__dict__["reads"][r]["align"] = os.path.join(tmp_dir, "alignment")
        return
    def count_features(self, ref, stranded, overlap, force_align = False, force_count = False):
        global mapper_run
        global sam_view
        global sam_sort
        for r in self.__dict__["reads"]:
            import pdb; pdb.set_trace()
            outbam = self[r]['align'] + ".bam"
            mapper_run.in_fq1 = self[r]["link"]
            sam_sort.out_prefix = self[r]['align']
            print '#' * 80
            if not os.path.isfile(outbam) or force_align:
                print '#### RUNNING BWA MEM on %s' % r
                run_align1 = subprocess.Popen(str(mapper_run).split(), stdout = subprocess.PIPE, bufsize = -1)
                out,err = run_align1.communicate()
                print '#' * 80
                print '#### RUNNING SAMTOOLS VIEW on %s' % r
                run_align2 = subprocess.Popen(str(sam_view).split(), stdin = subprocess.PIPE, stdout = subprocess.PIPE, bufsize = -1)
                out,err = run_align2.communicate(out)
                print '#' * 80
                print '#### RUNNING SAMTOOLS SORT on %s' % r
                run_align3 = subprocess.Popen(str(sam_sort).split(), stdin = subprocess.PIPE, bufsize = -1)
                out,err = run_align3.communicate(out)
            else:
                print '### Found alignment file for %s' % r
            print '#' * 80
            if 'counts' not in self.__dict__["reads"][r] or force_count or force_align:
                print '#### COUNTING READS MAPPED TO FEATURES on %s' % r
                self[r]['counts'] = count_reads_in_features(sam_filename = outbam,\
                                        gff_filename = ref,\
                                        samtype = 'bam',\
                                        order = 'name',\
                                        stranded = stranded,\
                                        overlap_mode = overlap,\
                                        feature_type = version, \
                                        id_attribute = 'ID',\
                                        minaqual = 60, \
                                        quiet = False, \
                                        samout = "")
            else:
                print '#### FOUND COUNTS FOR %s' % r
        print "#" * 80
        print "### Done counting."
        pickle_dump = os.path.join(self.workdir,".reads_dump")
        pickle_obj = open(pickle_dump, "w")
        pickle.dump(self, pickle_obj)
        pickle_obj.close()
    def make_pellet(self, include_features):
        if not os.path.isdir(self.results_path):
            os.makedirs(self.results_path)
        pellet_fn = os.path.join(self.results_path,"boobook_pellet.csv")
        read_sets = sorted(self.__dict__["reads"].keys())
        feature_list = sorted(self.__dict__["reads"][read_sets[0]]["counts"].keys())
        collated_counts = []
        for f in feature_list:
            feat_counts = [f]
            product = ""
            gene_id = ""
            gi = ""
            gene = ""
            tmp_qual = self.features[f]
            for feat in include_features:
                # if feat == 'tRNA':
                #     print feat
                #     print tmp_qual
                try:
                    quals = tmp_qual[feat].qualifiers
                    feat_type = feat
                    try:
                        product = '"' + quals["product"][0] + '"'
                    except:
                        pass
                    try:
                        gene = quals["gene"][0]
                    except:
                        pass
                    try:
                        tmp_db = quals["db_xref"]
                        try:
                            gi = [ref[4:] for ref in tmp_db if re.match("GI", ref)][0]
                        except:
                            pass
                        try:
                            gene_id = [ref[7:] for ref in tmp_db if re.match("GeneID", ref)][0]
                        except:
                            pass
                    except:
                        pass
                except:
                    pass
            feat_counts.extend([product, gene_id, gi, gene, feat_type])
            for r in read_sets:
                feat_counts.append(str(self.__dict__["reads"][r]["counts"][f]))
            feat_counts = ",".join(feat_counts) + "\n"
            collated_counts.append(feat_counts)
        header = ["locus_tag", "Product", "GeneID", "GI", "Gene", "Feature"]
        for r in read_sets:
            tmp = self.__dict__["reads"][r]['treat'] + "_" + \
                self.__dict__["reads"][r]['rep']
            header.append(tmp)
        header = ",".join(header) + "\n"
        # print header
        # for i in range(5):
        #     print collated_counts[i]
        pellet_con = open(pellet_fn, "w")
        pellet_con.write(header)
        pellet_con.write("".join(collated_counts))
        pellet_con.close()
        print "Pellet successfully created!"


# implementing boobook options and arguments
@click.command(context_settings=CONTEXT_SETTINGS)
# Boobook options
@click.option("--work_dir", "-i", \
                help = '''
                The project folder for boobook''', \
                default = ".")
@click.option("--features", \
                help = '''
                A colon separated list of features to count (e.g. CDS;snRNA)
                (default: CDS)
                ''', \
                default = "CDS")
@click.option("--qualifier", \
                help = '''
                A unique ID to identify distinct elements of the annotation
                (e.g., locus_tag)
                ''', \
                default = 'locus_tag')
# BWA options
@click.option("--bwa_threads", \
                help = "[BWA option] Number of threads for BWA mapping (default: 16)", \
                default = 16)
# Samtools options
@click.option("--sam_threads", \
                help = "[SAMTOOLS opion] Number of threads for Samtools viewing and sorting (default: 8)", \
                default = 8)
#HTSeq options
@click.option("--hts_stranded", \
                help = '''[HTSeq option] Strandedness of the RNAseq data ('yes', 'no','reversed')
                (default: 'reversed')''', \
                default = 'reversed')
@click.option("--hts_overlap", \
                help = '''[HTSeq option] How to account for overlapping features when counting overlapping
                reads in HTSeq ('union', 'intersection-strict', 'intersection-nonempty')
                The recommended mode is 'union' (default: union)
                ''', \
                default = 'union')
#Forcing things to get redone options
@click.option("--add", \
                help = '''
                Re-run analysis with additional samples in the infile or with
                additional features''', \
                is_flag = True,
                default = False)
@click.option("--re_align", \
                help = '''
                Re-do alignment of all samples in the project folder''', \
                is_flag = True,
                default = False)
@click.option("--re_count", \
                help = '''
                Re-do count of all samples in the project folder''', \
                is_flag = True,
                default = False)
@click.option("--change_ref", \
                help = '''Re-run analysis with a different reference. Will force
                re-alignment and re-count.
                ''',
                is_flag = True,
                default = False)
#The default arguments
@click.argument('REF')
@click.argument('INFILE')
def boobook(infile, ref, \
                    work_dir, \
                    features, \
                    qualifier, \
                    bwa_threads, \
                    sam_threads, \
                    hts_stranded, \
                    hts_overlap, \
                    add, \
                    re_align, \
                    re_count,
                    change_ref):
    '''Boobook is a tool to transform your RNAseq data into a table of counts
    suitable for analysis using Degust (http://www.vicbioinformatics.com/degust/).


    It requires two inputs:

        * REF: A Genbank reference file\n
        * INFILE: A TAB-delimited file with information about the reads

    Please post any issues to: https://github.com/MDU-PHL/boobook/issues

    Boobook is distributed under GPL v3.

    Author: Anders Goncalves da Silva
    '''
    #make sure alignment and counting is re-done if changing the reference
    if(change_ref):
        print("The change_ref flag was triggered...")
        print("Forcing re-alignment and counting...")
        re_align = True
        re_count = True

    # loading the read data information
    reads = ReadData(work_dir, force = add)
    reads.read_input(infile)
    reads.create_subfolders()

    # loading reference information
    # parse features
    feat = features.split(";")
    reference = GBk(work_dir, force = change_ref)
    reference.read_gb(ref)
    reference.filter_features(features = feat, qualifier = qualifier)
    reference.write_gff() # create a standardized GFF to use in HTSeq
    reference.write_fasta() # create a FASTA reference for use in alignment
    reads.features = reference.features # transfer the features dict to the
                                        # reads object. This dict will then
                                        # be called by a reads method for counting

    # setup alignment
    #create the commandline call to index the reference
    mapper_index = BWAindex(in_fasta = reference.fasta)
    #create the commandline call to align the reads
    mapper_run = BWAmem(threads = bwa_threads, \
                        ref = reference.fasta, \
                        in_fq1 = "dummy.fq") # just a dummy this will be overwritten
                                             # when looping through the reads
    #create the commandline call to transform the SAM output into BAM
    sam_view = SamtoolsViewCommandline(threads = sam_threads, \
                                        # minimum mapping quality
                                        q = 60,  \
                                        # SAM input format
                                        S = True, \
                                        # create BAM output
                                        b = True, \
                                        #uncompressed BAM output
                                        u = True, \
                                        ref = reference.fasta, \
                                        input = "-") # Take the input from stdin
    sam_sort = SamtoolsSortCommandline(threads = sam_threads, \
                                        # take input from stdin
                                        input_bam = "-", \
                                        out_prefix = "dummy_prefix") # this will be
                                                                # overwritten when
                                                                # looping through the
                                                                # reads
    # run the actual commands.
    # first, index the reference
    # the str() prints the command, and the .split() splits the command string
    # by white-space, creating a list with the individual components as needed
    # by the subprocess.Popen command
    run_index = subprocess.Popen(str(mapper_index).split())
    # wait for the indexing to end before doing anything else
    run_index.wait()
    reads.count_features(reference.gff, \
                            stranded = hts_stranded, \
                            overlap = hts_overlap, \
                            force_align = re_align, \
                            force_count = re_count)
    reads.make_pellet(feat)
    print('#' * 40)
    print("Boobook is done!")
    print("Boobook pellet with all your counts can be found in {}/{}/{}.".format(work_dir, \
                                                            "pellets", \
                                                            "boobook_pellet.csv"))
    return

if __name__ == "__main__":
    boobook()
    # debug = False
    # if debug and os.path.isdir("test/sample_counts/"):
    #     shutil.rmtree("test/sample_counts/")
    # reads = ReadData("test/", force = True)
    # reads.read_input("infile.tab")
    # reads.create_subfolders()
    # reference = GBk("test/", force = True)
    # reference.read_gb(infile = "test/Saa6008_final.gbk")
    # reference.filter_features(features = ['CDS','scRNA'], qualifier = 'locus_tag')
    # reference.write_gff()
    # reference.write_fasta()
    # reads.features = reference.features
    # # for k in reads.features.keys():
    # #     print reads.features[k]['CDS']
    # mapper_index = BWAindex(in_fasta = reference.fasta)
    # mapper_run = BWAmem(threads = 72, ref = reference.fasta, in_fq1 = reads["Sa_JKD6009_L30T_1"]["link"])
    # sam_view = SamtoolsViewCommandline(threads = 8, q=60, S=True, b = True, u = True, ref = reference.fasta, input = "-")
    # sam_sort = SamtoolsSortCommandline(threads = 8, input_bam = "-", out_prefix = reads["Sa_JKD6009_L30T_1"]["align"])
    # run_index = subprocess.Popen(str(mapper_index).split())
    # run_index.wait()
    # reads.count_features(reference.gff, force_align = True)
    # reads.make_pellet(['CDS','scRNA'])
