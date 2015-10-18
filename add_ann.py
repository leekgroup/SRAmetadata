#!/usr/bin/env python
"""
add_ann.py
Abhi Nellore
10/18/2015

Adds fields to all_SRA_introns.tsv.gz, which is read from stdin.
annotation files from command-line parameters and annotation from
GTF files specified as arguments of --annotations; writes to stdout.
all_SRA_introns.tsv.gz should have the following tab-separated fields on each
line:
1. chromosome
2. start position
3. end position
4. strand
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes in which junction was found
8. comma-separated list of numbers of reads in each corresponding sample from
    field 7

Fields added:
x,y,z ,where each is 1 or 0. x indicates whether the 5' splice site is in
annotation, y indicates whether the 3' splice site is in annotation, and
z indicates whether the junction is in annotation.

We executed:
gzip -cd all_SRA_introns.tsv.gz | pypy add_ann.py >[output_file]
"""
import os
import sys
import struct
import mmap
from operator import itemgetter
from bisect import bisect_right
from collections import defaultdict
import subprocess

class BowtieIndexReference(object):
    """
    Given prefix of a Bowtie index, parses the reference names, parses the
    extents of the unambiguous stretches, and memory-maps the file containing
    the unambiguous-stretch sequences.  get_stretch member function can
    retrieve stretches of characters from the reference, even if the stretch
    contains ambiguous characters.
    """

    def __init__(self, idx_prefix):

        # Open file handles
        if os.path.exists(idx_prefix + '.3.ebwt'):
            # Small index (32-bit offsets)
            fh1 = open(idx_prefix + '.1.ebwt', 'rb')  # for ref names
            fh3 = open(idx_prefix + '.3.ebwt', 'rb')  # for stretch extents
            fh4 = open(idx_prefix + '.4.ebwt', 'rb')  # for unambiguous seq
            sz, struct_unsigned = 4, struct.Struct('I')
        else:
            raise RuntimeError('No Bowtie index files with prefix "%s"'
                                    % idx_prefix)

        #
        # Parse .1.bt2 file
        #
        one = struct.unpack('<i', fh1.read(4))[0]
        assert one == 1

        ln = struct_unsigned.unpack(fh1.read(sz))[0]
        line_rate = struct.unpack('<i', fh1.read(4))[0]
        lines_per_side = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]
        ftab_chars = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]

        nref = struct_unsigned.unpack(fh1.read(sz))[0]
        # get ref lengths
        reference_length_list = []
        for i in xrange(nref):
            reference_length_list.append(struct.unpack('<i', fh1.read(sz))[0])

        nfrag = struct_unsigned.unpack(fh1.read(sz))[0]
        # skip rstarts
        fh1.seek(nfrag * sz * 3, 1)

        # skip ebwt
        bwt_sz = ln // 4 + 1
        line_sz = 1 << line_rate
        side_sz = line_sz * lines_per_side
        side_bwt_sz = side_sz - 8
        num_side_pairs = (bwt_sz + (2*side_bwt_sz) - 1) // (2*side_bwt_sz)
        ebwt_tot_len = num_side_pairs * 2 * side_sz
        fh1.seek(ebwt_tot_len, 1)

        # skip zOff
        fh1.seek(sz, 1)

        # skip fchr
        fh1.seek(5 * sz, 1)

        # skip ftab
        ftab_len = (1 << (ftab_chars * 2)) + 1
        fh1.seek(ftab_len * sz, 1)

        # skip eftab
        eftab_len = ftab_chars * 2
        fh1.seek(eftab_len * sz, 1)

        refnames = []
        while True:
            refname = fh1.readline()
            if len(refname) == 0 or ord(refname[0]) == 0:
                break
            refnames.append(refname.split()[0])
        assert len(refnames) == nref

        #
        # Parse .3.bt2 file
        #
        one = struct.unpack('<i', fh3.read(4))[0]
        assert one == 1

        nrecs = struct_unsigned.unpack(fh3.read(sz))[0]

        running_unambig, running_length = 0, 0
        self.recs = defaultdict(list)
        self.offset_in_ref = defaultdict(list)
        self.unambig_preceding = defaultdict(list)
        length = {}

        ref_id, ref_namenrecs_added = 0, None
        for i in xrange(nrecs):
            off = struct_unsigned.unpack(fh3.read(sz))[0]
            ln = struct_unsigned.unpack(fh3.read(sz))[0]
            first_of_chromosome = ord(fh3.read(1)) != 0
            if first_of_chromosome:
                if i > 0:
                    length[ref_name] = running_length
                ref_name = refnames[ref_id]
                ref_id += 1
                running_length = 0
            assert ref_name is not None
            self.recs[ref_name].append((off, ln, first_of_chromosome))
            self.offset_in_ref[ref_name].append(running_length)
            self.unambig_preceding[ref_name].append(running_unambig)
            running_length += (off + ln)
            running_unambig += ln

        length[ref_name] = running_length
        assert nrecs == sum(map(len, self.recs.itervalues()))

        #
        # Memory-map the .4.bt2 file
        #
        ln_bytes = (running_unambig + 3) // 4
        self.fh4mm = mmap.mmap(fh4.fileno(), ln_bytes,
                                flags=mmap.MAP_SHARED,
                                prot=mmap.PROT_READ)

        # These are per-reference
        self.length = length
        self.refnames = refnames

        # To facilitate sorting reference names in order of descending length
        sorted_rnames = sorted(self.length.items(),
                               key=lambda x: itemgetter(1)(x), reverse=True)
        self.rname_to_string = {}
        self.string_to_rname = {}
        for i, (rname, _) in enumerate(sorted_rnames):
            rname_string = ('%012d' % i)
            self.rname_to_string[rname] = rname_string
            self.string_to_rname[rname_string] = rname
        # Handle unmapped reads
        unmapped_string = ('%012d' % len(sorted_rnames))
        self.rname_to_string['*'] = unmapped_string
        self.string_to_rname[unmapped_string] = '*'

        # For compatibility
        self.rname_lengths = self.length

    def get_stretch(self, ref_id, ref_off, count):
        """
        Return a stretch of characters from the reference, retrieved
        from the Bowtie index.

        @param ref_id: name of ref seq, up to & excluding whitespace
        @param ref_off: offset into reference, 0-based
        @param count: # of characters
        @return: string extracted from reference
        """
        assert ref_id in self.recs
        stretch = []
        starting_rec = bisect_right(self.offset_in_ref[ref_id], ref_off) - 1
        assert starting_rec >= 0
        off = self.offset_in_ref[ref_id][starting_rec]
        buf_off = self.unambig_preceding[ref_id][starting_rec]
        '''Naive to scan these records linearly; obvious speedup is binary
        search'''
        for rec in self.recs[ref_id][starting_rec:]:
            off += rec[0]
            while ref_off < off and count > 0:
                stretch.append('N')
                count -= 1
                ref_off += 1
            if count == 0:
                break
            if ref_off < off + rec[1]:
                # stretch extends through part of the unambiguous stretch
                buf_off += (ref_off - off)
            else:
                buf_off += rec[1]
            off += rec[1]
            while ref_off < off and count > 0:
                buf_elt = buf_off >> 2
                shift_amt = (buf_off & 3) << 1
                stretch.append(
                    'ACGT'[(ord(self.fh4mm[buf_elt]) >> shift_amt) & 3]
                )
                buf_off += 1
                count -= 1
                ref_off += 1
            if count == 0:
                break
        # If the requested stretch went past the last unambiguous
        # character in the chromosome, pad with Ns
        while count > 0:
            count -= 1
            stretch.append('N')
        return ''.join(stretch)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--extract-splice-sites-path', type=str,
        default=os.path.join(os.path.dirname(__file__),
                                'extract_splice_sites.py'),
        help=('path to extract_splice_sites.py from HISAT v0.1.6-beta.'))
    parser.add_argument('--annotations', type=str, required=True, nargs='+',
        help='paths to GTF files encoding known junctions')
    parser.add_argument('--bowtie1-idx', type=str, required=True,
            help=('Path to basename of Bowtie 1 index with genome to which '
                  'FASTQs were aligned')
        )
    args = parser.parse_args()

    possible_combos = set([('GT', 'AG'), ('CT', 'AC'),
                            ('GC', 'AG'), ('CT', 'GC'),
                            ('AT', 'AC'), ('GT', 'AT')])
    plus_possibles = set([('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')])
    minus_possibles = set([('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')])

    annotated_junctions = set()
    annotated_5p = set()
    annotated_3p = set()
    refs = ['chr' + str(i) for i in xrange(1, 23)] + ['chrM', 'chrX', 'chrY']
    reference_index = BowtieIndexReference(args.bowtie1_idx)
    for annotation in args.annotations:
        extract_process = subprocess.Popen([sys.executable,
                                                args.extract_splice_sites_path,
                                                annotation],
                                                stdout=subprocess.PIPE)
        for line in extract_process.stdout:
            tokens = line.strip().split('\t')
            tokens[1] = int(tokens[1]) + 2
            tokens[2] = int(tokens[2])
            if not tokens[0].startswith('chr'):
                tokens[0] = 'chr' + tokens[0]
            if tokens[0] in refs:
                annotated_junctions.add(tuple(tokens[:-1]))
                left = reference_index.get_stretch(tokens[0], tokens[1] - 1, 2)
                right = reference_index.get_stretch(
                                    tokens[0], tokens[2] - 3, 2
                                )
                if (left, right) not in possible_combos: continue
                itsplus = True
                if (left, right) in plus_possibles:
                    annotated_5p.add(tuple(tokens[:-2]))
                    annotated_3p.add((tokens[0], tokens[2]))
                else:
                    assert (left, right) in minus_possibles
                    annotated_5p.add((tokens[0], tokens[2]))
                    annotated_3p.add(tuple(tokens[:-2]))
                    itsplus = False
                if itsplus != (tokens[3] == '+'):
                    print >>sys.stderr, ('extract_splice_sites sign disagrees '
                                         'with sign from reference')
        extract_process.stdout.close()
        exit_code = extract_process.wait()
        if exit_code != 0:
            raise RuntimeError(
                'extract_splice_sites.py had nonzero exit code {}.'.format(
                                                                    exit_code
                                                                )
            )

    for line in sys.stdin:
        tokens = line.strip().split('\t')
        left = int(tokens[1])
        right = int(tokens[2])
        if tokens[3] == '+':
            threep = right
            fivep = left
        else:
            fivep = right
            threep = left
        if (tokens[0], fivep) in annotated_5p:
            x = '1'
        else:
            x = '0'
        if (tokens[0], threep) in annotated_3p:
            y = '1'
        else:
            y = '0'
        if (tokens[0], left, right) in annotated_junctions:
            z = '1'
        else:
            z = '0'
        print '\t'.join([line.strip(), x, y, z])