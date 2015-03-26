####################################################################
# Functions for reading mirbase files
#
# author: Gregor Sturm (gregor.sturm@cs.tum.edu)
###################################################################

from read_results import *
from BCBio import GFF
import sys

class MatureReader:
    def __init__ (self, mirfile, gfffile):
        """read mature miRNA file from mirbase fasta-file and the
        genome.gff file."""
        self.mature = read_fasta_to_dict(mirfile)
        properties_by_name = {}
        properties_by_id = {}
        handle = gfffile

        """read ids"""
        limit_info = dict(gff_type = ["miRNA_primary_transcript"])
        for rec in GFF.parse(handle, limit_info=limit_info):
            for mirna in rec.features:
                name = mirna.qualifiers["Name"][0].lower()
                ID = mirna.qualifiers["ID"][0]
                properties = {}
                properties["name"] = name
                properties["ID"] = ID
                properties["location"] = mirna.location
                properties_by_name[name] = properties
                properties_by_id[ID] = properties

        """read matures annotations"""
        handle.seek(0)
        limit_info = dict(gff_type = ["miRNA"])
        for rec in GFF.parse(handle, limit_info=limit_info):
            for mirna in rec.features:
                matname = mirna.qualifiers["Name"][0].lower()
                derived_from = mirna.qualifiers["Derives_from"][0]
                properties = properties_by_id[derived_from]
                annotation_type = "unknown"
                if matname.endswith("-3p"):
                    annotation_type = "3p"
                elif matname.endswith("-5p"):
                    annotation_type = "5p"
                properties[annotation_type] = matname
                properties[annotation_type + "_location"] = mirna.location

        self.properties_by_name = properties_by_name

    def find_mature(self, seq, name, search_start=0, search_end=None):
        """get tuple ((5p-start, 5p-end), (3p-start, 3p-end)) if name
        is found in seq.

        Search is done with "intelligent matching". The annontated strand
        is rejected and the strand determined by starting position.

        Args:
            seq: precursor sequence (exactly the precursor sequence!)
            name: mirbase-name
            search_start: position from which on the prec seq should be searched
                for the mature sequence. (useful for large, repetitive
                flanking regions)

        Returns:
            pos5p, pos3p
            where each of those is a tuple (start, end). End refers to
            the last position of the mature sequence.

        """

        if search_end is None:
            search_end = len(seq)

        properties = self.properties_by_name.get(name)
        assert properties is not None, "name not in database"

        name5p = properties.get("5p")
        name3p = properties.get("3p")
        name_unknown = properties.get("unknown")

        assert not (name5p is None and name3p is None and name_unknown is None), """
            at least one sequence annotated"""

        assert not (name5p is not None and name3p is not None
                and name_unknown is not None), """not all three available"""

        assert not ((name5p is not None and name_unknown is not None) or
                (name3p is not None and name_unknown is not None)), """
            when one explicit given, there is no implicit entry"""


        def getPosition(properties, annotation_type, seq):
            mature_name = properties.get(annotation_type)
            #mature_loc = properties.get(annotation_type + "_location")
            #prec_loc  = properties.get("location")
            #if prec_loc.strand > 0:
            #    offset = mature_loc.start - prec_loc.start
            #else:
            #    """ reverse strand """
            #    offset = prec_loc.end - mature_loc.end

            #""" offset should be the value, but let's check """

            mature = self.mature.get(mature_name)
            assert mature is not None,mature_name
            uniq = seq.count(mature, search_start, search_end)
            if uniq != 1:
                print >>sys.stderr, "Read not matching uniquely %r" % name

            start = seq.find(mature, search_start, search_end)
#            assert start == offset, (start, offset)
            return (start, start + len(mature) - 1)

        if name_unknown is not None:
            pos_unknown = getPosition(properties, "unknown", seq)
            """ guess strand """
            if pos_unknown[0] < .5 * len(seq) - 10:
                return pos_unknown, (-1, -1)
            else:
                return (-1, -1), pos_unknown
        else:
            pos5p = pos3p = (-1, -1)
            if name5p is not None:
                pos5p = getPosition(properties, "5p",  seq)
            if name3p is not None:
                pos3p = getPosition(properties, "3p",  seq)
                assert pos5p[0] < pos3p[0]
            return pos5p, pos3p


