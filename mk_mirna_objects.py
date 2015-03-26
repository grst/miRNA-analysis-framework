#!/usr/bin/env python
# -*- coding: utf-8 -*

__description__ = """Generates miRNA objects and pickles them for further usage"""

from read_results import *
from read_mirs import MatureReader
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from BCBio import GFF
import pickle
import argparse

class MirnaFactory(object):
    """reads the raw mirbase files and provides methods for
    fast access to the data, given a miRNA name"""

    VERBOSE = True
    ADD_EXT = 60

    def __init__(self, precursor_file, mature_file,
            coordinate_file, genome_path, constraint_file=None, **kwargs):
        """init the DataReader with filenames to the miRBase files.

        Args:
            precursor_file: fasta-file with the precursor sequence from mirbase
            mature_file: fasta-file with the mature sequences
            coordinate_file: gff3-file with genomic coordinates of the precursors
            genome_path: path to the directory with chromosomres.txt
            constraint_file: file in fasta format with secondary structure
                constraints in dotbracket notation (Optional)
            extension_length: Nucleotides to be added
                to the precursor sequence (Optional)
            annotation_files: Dictionary annotation:file_containing_ids. (Optional)
                Can be used to annotate the miRNAs based on a list of ids.
            remove_constraints: If set to True, the miRBase secondary structure
                constraints will be removed in the lower stem region. (Optinal,
                default False)
        """

        if constraint_file is not None:
            self.dotbracket = read_fasta_to_dict(constraint_file)
        else:
            self.dotbracket = None
        self.precursors = read_fasta_to_dict(precursor_file)
        self.mature_reader = MatureReader(mature_file, coordinate_file)
        self._mk_coordinate_dict(coordinate_file)
        self.genome_path = genome_path
        self.annotations = {}

        #additional nucleotides
        self.extension_length = 60
        if kwargs.get("extension_length") is not None:
            self.extension_length = kwargs["extension_length"]

        #remove constraints beyond cleavage site?
        self.remove_constraints = False
        if kwargs.get("remove_constraints") is not None:
            self.remove_constraints = kwargs["remove_constraints"]

        #process annotation files
        if kwargs.get("annotation_files") is not None:
            for keyword, annotation_file in kwargs["annotation_files"].iteritems():
               self.annotations[keyword] = read_to_dict(annotation_file)

    def _mk_coordinate_dict(self, gff_file):
        handle = gff_file
        handle.seek(0)
        limit_info = dict(gff_type = ["miRNA_primary_transcript"])
        self.coordinates = {}


        for rec in GFF.parse(handle, limit_info=limit_info):
            chromosome = rec.id
            for mirna in rec.features:
                name = mirna.qualifiers["Name"][0].lower()
                location = mirna.location

                self.coordinates[name] = {
                        "chromosome" : chromosome,
                        "location" : location
                }


    def produce_all(self, ids):
        """ create a hashmap with the correpsonding
        mirna objects for the names given in ids.

        Does perform the sequence extension based on the
        coordinates
        """
        all_mirnas = {}

        from hairpin import AnnotatedHairpin

        for name in ids:
            mirna = self.produce_mirna(name)
            all_mirnas[name] = mirna

        self._extend_mirnas(all_mirnas)

        for mirna in all_mirnas.itervalues():
            print "checking integrity %s" % mirna.name
            try:
                mirna.check_integrity()
            except AssertionError, msg:
                print >> sys.stderr, msg
                mirna.add_annotation("extension_fail")

        return all_mirnas


    def produce_mirna(self, name):
        """ create a mirna object by name, however does not perform
        the nucleotide-extension
        """
        name = name.lower()
        seq = self.precursors[name]

        if self.VERBOSE:
            print "processing: ", name
            print seq

        mature_5p, mature_3p = self.mature_reader.find_mature(seq, name)

        if self.VERBOSE:
            print (mature_5p, mature_3p)

        constraint = None if self.dotbracket is None else self.dotbracket.get(name)

        assert constraint is not None

        mirna = AnnotatedHairpin(name,
                seq,
                None if mature_5p[0] < 0 else mature_5p[0],
                None if mature_3p[1] < 0 else mature_3p[1] + 1,
                constraint=constraint)

        for keyword, id_list in self.annotations.iteritems():
            if id_list.get(name) is not None:
                mirna.add_annotation(keyword)

        return mirna


    def _extend_mirnas(self, mirnas):
        """fetch context information from the genome for
        calculation of features that are note within the
        actual pre-mirRNA hairpin sequence

        This method involves a hack:
            Extend ALL miRNAs (also those that have enough context)
            by a fixed (longer) length. After that trim all
            miRs to the desired length. (This fixes problems
            with e.g. 6068, where the cleavage site is not
            withing the miRBase annotated precrsor sequence)

        a lot of asserting here to avoid mistakes and to find
        wrongly annotated miRNAs
        """

        for name, mirna in mirnas.iteritems():
            chromosome = self.coordinates[name]["chromosome"]
            location = self.coordinates[name]["location"]
            with open(self.genome_path + chromosome + ".txt", "r") as f:
                tmp_ext = self.extension_length + self.ADD_EXT
                #correct position found?
                f.seek(location.start)
                ref = f.read(location.end - location.start).upper()
                if location.strand < 0:
                    ref = Seq(ref, generic_dna).reverse_complement()
                assert mirna.dna_seq == ref

                print name, " : extending"

                if self.VERBOSE:
                    print name
                    print "mirbase-sequence:"
                    print ref
                    print mirna.constraint
                    print mirna.cleavage_site

                #control: hairpin-sequence
                prec_seq_before = mirna.dna_seq[mirna.cleavage_site[0] :
                        mirna.cleavage_site[1]]

                #location.end refers to the position of the last character,
                #not the position after!
                cs_coords_5p = location.start
                cs_coords_3p = location.end

                #cleavage site coordinates correct
                assert (cs_coords_3p - cs_coords_5p  == mirna.length)

                cs_coords_5p -= tmp_ext
                #this refers to the position of the +60th base
                cs_coords_3p += tmp_ext

                f.seek(cs_coords_5p)
                #this returns exactly the sequence, including the +60th base
                extended_seq = f.read(cs_coords_3p - cs_coords_5p).upper()

                if location.strand < 0:
                    print name + " : reversed."
                    extended_seq = str(Seq(extended_seq,
                        generic_dna).reverse_complement())

                #length of extended sequence correct? (length nts from each CS)
                assert (len(extended_seq) == 2*tmp_ext
                        + mirna.length)

                #fill the dotbracket seqs beyond the cleavage site with dots
                extended_dotbracket = ("." * tmp_ext
                        + mirna.constraint
                        + "." * tmp_ext)

                if self.VERBOSE:
                    print "extended-sequence:"
                    print extended_seq
                    print extended_dotbracket

                #length of dotbracket equals length of sequence
                assert len(extended_dotbracket) == len(extended_seq)

                #trim length to 60nts beyond the cleavage site.
                tmp_co = (tmp_ext + mirna.cleavage_site[0]
                    - self.extension_length ,
                    tmp_ext + mirna.cleavage_site[1] + self.extension_length)
                extended_seq = extended_seq[tmp_co[0] : tmp_co[1]]
                extended_dotbracket = extended_dotbracket[tmp_co[0] : tmp_co[1]]
                if self.VERBOSE: print "trim at", tmp_co

                if self.VERBOSE:
                    print "trimmed-sequence:"
                    print extended_seq
                    print extended_dotbracket

                #length of trimmed sequence correct?
                assert len(extended_seq) == (2*self.extension_length
                    + mirna.cleavage_site[1] - mirna.cleavage_site[0])


                #update cleavage site, intentionally change private attribute
                mirna._cleavage_site = (self.extension_length,
                    self.extension_length
                    + mirna.cleavage_site[1] - mirna.cleavage_site[0])
                if self.VERBOSE: print "new CS", mirna.cleavage_site

                #remove constraints beyond the cleavage sites
                if self.remove_constraints:
                    extended_dotbracket = self._remove_constraints(
                            extended_dotbracket, mirna.cleavage_site)

                prec_seq_after = extended_seq[mirna.cleavage_site[0] :
                        mirna.cleavage_site[1]]

                if self.VERBOSE:
                    print "hairpin-sequence"
                    print prec_seq_before
                    print prec_seq_after

                #correctness of precursor sequences
                #equality is not required any more, since prec_before
                #can be shorter, when the annotated precursor is too short
                assert len(prec_seq_before) <= len(prec_seq_after)

                #update annotation of mature sequences
                mature_5p, mature_3p = self.mature_reader.find_mature(
                       extended_seq, mirna.name,
                       self.extension_length,
                       self.extension_length + len(prec_seq_after))

                new_mirna = AnnotatedHairpin(name, extended_seq,
                        mirna.cleavage_site[0], mirna.cleavage_site[1],
                        constraint=extended_dotbracket)
                new_mirna._annotation = mirna.annotation[:]

                assert new_mirna.cleavage_site[0] == self.extension_length
                assert (new_mirna.length - new_mirna.cleavage_site[1]
                        == self.extension_length)

                mirnas[name] = new_mirna

                if self.VERBOSE:
                    print "annotation", ",".join(mirna.annotation)
                    print "#" * 60
                    print ""

    @staticmethod
    def _remove_constraints(dotbracket, cs):
        """ remove constraints that are beyond the cleavage site
        in the lower stem as a matter of consistency
        """
        if MirnaFactory.VERBOSE:
            print "removing constraints"
            print dotbracket

        length = len(dotbracket)
        dotbracket = (
                "." * (cs[0]) +
                dotbracket[cs[0]:cs[1]] +
                "." * (len(dotbracket) - cs[1])
        )

        #balance dotbracket
        l = dotbracket.count("(")
        r = dotbracket.count(")")
        if l > r:
            """remove from left"""
            dotbracket = dotbracket.replace("(", ".", l -r)
        elif r > l:
            """remove from right"""
            rev = dotbracket[::-1]
            rev = rev.replace(")", ".", r -l)
            dotbracket = rev[::-1]

        if MirnaFactory.VERBOSE:
            print dotbracket

        assert(len(dotbracket) == length)

        return dotbracket

if __name__ == "__main__":

    argp = argparse.ArgumentParser(description=__description__)
    argp.add_argument('-p', '--precursor', type=file, action="store",
            dest="precursor_file", required=True,
            help="precursor file as downloaded from mirbase")
    argp.add_argument('-c', '--constraint', type=file, action="store",
            dest="dotbracket_file", required=False, default=None,
            help="additional folding constraints for miRNAs in dotbracket format")
    argp.add_argument('-m', '--mature', type=file, action="store",
            dest="mature_file", required=True,
            help="mature sequences from mirbase")
    argp.add_argument('-g', '--gff', type=file, action="store",
            dest="coordinate_file", required=True,
            help="gff file from mirbase with genome coordinates")
    argp.add_argument('-t', '--genome-path', type=str, action="store",
            dest="genome_path", required=True,
            help="path to folder with chromosome-text-files")
    argp.add_argument('-e', '--extension', type=int, action="store",
            dest="extension", required=False, default=30,
            help="the miRNAs will be extended by N nucleotides beyond the CS")
    argp.add_argument('-o', '--output_pickle', type=argparse.FileType('w'),
            action="store", dest="pickle_file", required=True,
            help="pickle mirna objects to this file")
    argp.add_argument('-f', '--output_fasta', type=argparse.FileType('w'),
            action="store", dest="output_file", required=True,
            help="write extended sequences to this fasta files")
    argp.add_argument('-a', '--annotation', type=file,
            action="store", dest="annotation_files", required=False,
            nargs='+',
            help="annotation files containing one mirna id per line")
    argp.add_argument('-l', '--annotation_labels', type=str,
            action="store", dest="annotation_labels", required=False,
            nargs='+',
            help="labels for annotation in the same order as annotation_files")

    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))

    assert len(args.annotation_files) == len(args.annotation_labels),\
            "specify a label for every annotation file"

    annotation_files = dict(zip(args.annotation_labels, args.annotation_files))

    factory = MirnaFactory(args.precursor_file,
           args.mature_file,
           args.coordinate_file,
           args.genome_path,
           extension_length=args.extension,
           annotation_files=annotation_files,
           constraint_file=args.dotbracket_file,
           remove_constraints=True)

    ids = factory.precursors.keys()

    all_mirs = factory.produce_all(ids)

    pickle.dump(all_mirs, args.pickle_file)

    for name, mirna in all_mirs.iteritems():
        print >>args.output_file, ">" + name
        print >>args.output_file, mirna.dna_seq

