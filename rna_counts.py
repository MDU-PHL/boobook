'''
A Python tool to count reads mapped to genomic features
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.GenBank.Scanner import GenBankScanner
from Bio.Sequencing.Applications import _samtools as samtools
from operator import itemgetter
import re
import sys
import time

class GBk:
    '''
    A class that reads and stores GenBank
    '''
    def __init__(self):
        pass
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
        if (not isinstance(features, list)):
            features = [features]
        assert type(qualifier) is str, "qualifier is not a string: %r" % qualifier
        self.features = {}
        count_found_features = 0
        count_found_features_with_qualifiers = 0
        count_per_feature = dict((f, 0) for f in features)
        for (index, feat) in enumerate(self.gb.features):
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
        return
    def write_gff(self, outfile):
        '''
        Write a GFF output file that can be read by htseq-count
        '''
        f = open(outfile, 'w')
        x = len(self.features.keys())
        count1 = 0
        count2 = 0
        countPlus = 0
        refname = self.gb.id
        source = 'RefSeq'
        feature = 'rna-count:0.1'
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
            ltag = tmp_rec.qualifiers['locus_tag'][0]
            coords = tmp_rec.location
            (start, end, strand) = self.__parse_location(coords)
            attrib = "ID="+ltag
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
        out_str += "##source-version rna_counts.py 0.01\n"
        out_str += "##data This file was mode on " + time.strftime("%d/%m/%Y") + "\n"
        out_str += "##description " + self.gb.description + "\n"
        out_str += "##sequence-region " + refname + " " + "1 " + str(len(self.gb.seq)) + "\n"
        for r in sorted_gff_recs:
            out_str+="\t".join(map(str,r))+'\n'
        f.write(out_str)
        f.close()
        return
    def write_fasta(self, filename):
        '''
        This function outputs the sequence data in the GenBank file to a
        FASTA file called ref.fa
        '''
        f = open(filename, "w")
        SeqIO.write(self.gb, f, 'fasta')
        f.close()
        print "Reference file was successfully written."
        return


class GFF:
    def __init__(self):
        return

class ReadData:
    def __init__(self):
        return

class BWA_mem:
    def __init__(self):
        return

if __name__ == "__main__":
    data = GBk()
    data.read_gb(infile = "test/Sa6008.gbk")
    data.filter_features(features = ['CDS', 'tRNA'], qualifier = 'locus_tag')
    data.write_gff("test/test.gff")
    data.write_fasta("test/ref.fa")
