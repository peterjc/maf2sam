#!/usr/bin/env python
"""Simple MIRA alignment format (MAF) to SAM format converter.

See: http://mira-assembler.sourceforge.net/docs/chap_maf_part.html
and: http://samtools.sourceforge.net/

The source code repository for this script is here:

http://github.com/peterjc/maf2sam

Copyright 2010-2011, Peter Cock, all rights reserved.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
"""
#
#v0.0.0 - Uses the gapped co-ordindates
#v0.0.1 - Use the ungapped co-ordindates
#v0.0.2 - Use object for each read
#v0.0.3 - Fill in pair partner info
#v0.0.4 - Simple command line interface
#       - Limited release 28 September 2010 by email
#v0.0.5 - Cope with BR, IC and IR lines in reads
#       - Limited release 30 September 2010 by email
#v0.0.6 - Refactor to cope with other MAF read lines
#       - Public release 30 September 2010 on MIRA mailing list
#v0.0.7 - Use stderr for missing Biopython error message
#       - Dated 3 November 2010
#v0.0.8 - Use handles with Bio.SeqIO to support pre-Biopython 1.54
#         (tested this on Biopython 1.47 and seems fine).
#v0.0.9 - Mate's RNAME and POS (known as MRNM and MPOS in SAM v1.2, renamed
#         as RNEXT and PNEXT in SAM v1.3) should default to * and 0.
#v0.0.10- Do not assume read names start with template name
#         (MIRA can be given this information explicitly in XML input)
#       - Ignores new BC line type in read blocks
#v0.0.11- Update CIGAR strings to use = and X (match and mismatch)
#         rather than just M (either), as per SAM v1.3 onwards.
#         WARNING - This is supported in samtools 0.1.18 onwards
#       - Converts sequences to upper case (since case is meaningless
#         in the SAM format, and not preserved in BAM format).
#v0.0.12- Set the SAM/BAM properly paired flag
#v0.1.00- Pre-parse the MAF file to get ST (sequence technology) and
#         SN (strain) lines in order to write @RG lines for the header
#         (using @RG PL platform and SM sample tags).
#       - Report @HD VN:1.4, i.e. we try to follow SAM spec v1.4
#v0.1.01- Record MD5 digest in @SG lines.
#
#PRERELEASE:
#v0.2.00- Produce either gapped or ungapped (padded or unpadded) SAM
#       - Internal option to produce CIGAR strings using M
#         (for testing with bits of samtools which don't like X/=)
#       - Record dummy reads for consensus/reference contig sequences
#         (as defined by Heng Li in pre-release SAM/BAM specification)
#       - Report file format version as 1.5
#       - Use P operators in CIGAR strings (unpadded SAM)
#       - Record MIRA's CT annotation using dummy reads with RT tags
#         (note CT lines before CS line so have to cache them).
#       - Record MIRA's RT annotation using PT tags
#
#
#TODO
# - Extend pre-parsing to record read offsets in file, so that we can
#   produce a sorted SAM file?
# - Could read contigs from MAF file itself? (On the other hand, the user
#   will need the unpadded reference FASTA to use the SAM output anyway)
# - Rewrite to avoid Biopython requirement?
# - insert size
# - Record original read name suffix in tags
# - Record any MIRA annotation in tags?
# - more testing

import sys
import re
import hashlib

CIGAR_M = True
RECORD_CT = True

if len(sys.argv)==3:
    ref = sys.argv[1]
    maf = sys.argv[2]
else:
    import os
    name = os.path.basename(sys.argv[0])
    print "Usage: Takes two command line arguments, (un)padded FASTA reference"
    print "file, and matching MAF file. Output is SAM format to stdout."
    print
    print "python %s EX_out.unpadded.fasta EX_out.maf > EX_out.unpadded.sam" % name
    print
    print "or if the script is marked as executable,"
    print
    print "./%s EX_out.unpadded.fasta EX_out.maf > EX_out.unpadded.sam" % name
    print
    print "Note that conventionally SAM/BAM have used an unpadded (ungapped)"
    print "reference sequence, but as of September 2011 this has been relaxed"
    print "to allow a padded (gapped) reference. This is intended to be used"
    print "for (de novo) assemblies, not mapping to a reference."
    print
    print "You can now produce either style SAM file, depending on if your"
    print "reference FASTA sequence is gapped/padded or not."
    print
    print "NOTE - This script does not accept ACE files as input."
    sys.exit(1)

try:
    from Bio.Seq import reverse_complement
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Requires Biopython\n")
    sys.exit(1)

def log(msg):
    sys.stderr.write("[maf2sam] %s\n" % msg.rstrip())

def seq_md5(seq):
    #TODO - Remove these asserts after testing
    assert " " not in seq
    assert seq == seq.upper()
    return hashlib.md5(seq).hexdigest()

class Read(object):
    def __init__(self, contig_name, read_name="", template_name="",
                 read_seq="", first_in_pair=True, ref_rc = False,
                 ref_pos=0, map_qual=255, insert_size = 0,
                 vect_left = 0, vect_right = 0,
                 qual_left = 0, qual_right = 0,
                 clip_left = 0, clip_right = 0,
                 seq_tech = "", strain = "",
                 tags=None, annotations=None):
        self.contig_name = contig_name
        self.read_name = read_name
        self.template_name = template_name
        self.read_seq = read_seq
        self.first_in_pair = first_in_pair
        self.ref_rc = ref_rc
        self.ref_pos = ref_pos
        self.map_qual = map_qual
        self.cigar = ""
        self.insert_size = insert_size
        self.vect_left = vect_left
        self.vect_right = vect_right
        self.qual_left = qual_left
        self.qual_right = qual_right
        self.clip_left = clip_left
        self.clip_right = clip_right
        self.seq_tech = seq_tech
        self.strain = strain
        if tags:
            self.tags = tags
        else:
            self.tags = []
        if annotations:
            self.annotations = annotations
        else:
            self.annotations = []
    
    def __repr__(self):
        return "Read(%r, %r, %r, %r, %r, %r, %r, %r, ...)" % (
            self.contig_name,
            self.read_name,
            self.template_name,
            self.read_seq,
            self.first_in_pair,
            self.ref_rc,
            self.ref_pos,
            self.map_qual)

    def is_paired(self):
        if not self.template_name:
            assert self.read_name
            return False
            self.template_name = self.read_name
        elif self.template_name != self.read_name:
            #Looks like a paired end read!
            return True
        else:
            return False
    
    def get_partner(self):
        if not self.is_paired():
            raise ValueError
        return cached_pairs[(self.template_name, not self.first_in_pair)]
    
    def need_partner(self):
        if not self.is_paired():
            return False
        return (self.template_name, not self.first_in_pair) not in cached_pairs
        
    def __str__(self):
        global cached_pairs, read_group_ids
        if self.ref_rc:
            flag = 0x10 #maps to reverse strand
            read_seq = reverse_complement(self.read_seq)
            read_qual = self.read_qual[::-1]
        else:
            flag = 0
            read_seq = self.read_seq
            read_qual = self.read_qual
        mate_ref_name = "*"
        mate_ref_pos = 0
        if not self.template_name:
            assert self.read_name
            self.template_name = self.read_name
        if self.is_paired():
            flag += 1 #paired
            if self.first_in_pair:
                flag += 0x40 #forward partner
            else:
                flag += 0x80 #reverse partner
            try:
                mate = self.get_partner()
            except KeyError:
                #Paired but no parter in ACE file
                flag += 0x08 #mate unmapped
            else:
                mate_ref_name = mate.contig_name
                mate_ref_pos = mate.ref_pos
                if mate_ref_name == self.contig_name:
                    #Since MIRA seems happy and both on same contig,
                    flag += 0x02 #properly aligned

        assert not self.tags
        read_seq_unpadded = read_seq.replace("*", "")
        read_qual_unpadded = "".join(q for (l,q) in zip(read_seq,read_qual) if l!="*")
        cigar = self.cigar
        #assert "M" not in cigar, cigar
        if "D" not in cigar and "P" not in cigar:
            #Sum of lengths of the M/I/S/=/X operations should match the sequence length
            #By construction there are no M entries in our CIGAR string.
            #TODO - Improve this check to consider D in CIGAR?
            if len(read_seq_unpadded) != sum(int(x) for x in cigar.replace("I","=").replace("S","=").replace("M","X").replace("X","=").split("=") if x):
                raise ValueError("%s vs %i for %s" % (cigar, len(read_seq_unpadded), read_seq))
        assert len(read_seq_unpadded) == len(read_qual_unpadded)
        line = "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\t%s\t%s" % \
            (self.template_name, flag, self.contig_name, self.ref_pos,
             self.map_qual, cigar,
             mate_ref_name, mate_ref_pos, self.insert_size,
             read_seq_unpadded, read_qual_unpadded)
        assert self.seq_tech
        line += "\tRG:Z:%s" % read_group_ids[(self.seq_tech, self.strain)]
        for tag in self.tags:
             assert not tag.startswith("RG:"), tag
             assert not tag.startswith("PT:"), tag
             line += "\t" + tag
        if self.annotations:
            annotations = []
            for start, end, tag, value in self.annotations:
                if start <= end:
                    strand = "+"
                else:
                    strand = "-"
                    start, end = end, start
                #These should already be 1-based padded reference coords
                assert 1 <= start <= end <= len(self.read_seq), \
                    "Problem with %s PT tag coordindates %i:%i (bounds 1:%s) for %s %s" \
                    % (self.read_name, start, end, len(self.read_seq), tag, value)
                annotations.append("%i|%i|%s|%s|%s" \
                                   % (start, end, strand, tag, value))
            line += "\tPT:Z:%s" % "|".join(annotations)
        return line

print "@HD\tVN:1.5\tSO:unsorted"
print "@CO\tConverted from a MIRA Alignment Format (MAF) file"

ref_lens = {}
ref_md5 = {}
handle = open(ref)
gapped_sam = False
for rec in SeqIO.parse(handle, "fasta"):
    #Note MIRA uses * rather than - in the output padded FASTA
    #However, for padded references SAM/BAM say use * for MD5
    seq = rec.seq.tostring().upper().replace("-","*")
    if not gapped_sam and "*" in seq:
        log("NOTE: Producing SAM using a gapped reference sequence.")
        gapped_sam = True
    md5 = seq_md5(seq)
    ref_md5[rec.id] = md5
    ref_lens[rec.id] = len(seq)
    print "@SQ\tSN:%s\tLN:%i\tM5:%s" % (rec.id, len(seq), md5)
handle.close()
if not ref_lens:
    log("No FASTA sequences found in reference %s" % ref)
    sys.exit(1)

#First pass though the MAF file to get info for read groups
seq_tech_strains = set() #will make into a list of 2-tuples
handle = open(maf)
tech = ""
strain = ""
for line in handle:
    if line.startswith("RD"):
        assert not tech and not strain
    elif line.startswith("ST\t"):
        tech = line[3:].strip()
    elif line.startswith("SN\t"):
        strain = line[3:].strip()
    elif line.startswith("ER"):
        seq_tech_strains.add((tech, strain))
        tech = ""
        strain = ""
handle.close()
seq_tech_strains = sorted(list(seq_tech_strains))
read_group_ids = dict()
for id, (tech, strain) in enumerate(seq_tech_strains):
    platform = tech.upper()
    if platform == "SANGER":
        platform = "CAPILLARY"
    elif platform == "SOLEXA":
        platform = "ILLUMINA"
    elif platform == "454":
        platform = "LS454"
    if platform not in ["CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "PACBIO"]:
        raise ValueError("Sequencing technology (ST line) %r not supported in SAM/BAM" % tech)
    assert len(strain.split())<=1, "Whitespace in strain %r (SN line)" % strain
    read_group_id = ("%s_%s" % (tech, strain)).strip("_")
    read_group_ids[(tech, strain)] = read_group_id
    print "@RG\tID:%s\tPL:%s\tSM:%s" % (read_group_id, platform, strain)
del strain, tech, seq_tech_strains
log("Identified %i read groups" % len(read_group_ids))

def decode_cigar(cigar):
    """Returns a list of 2-tuples, integer count and operator char."""
    count = ""
    answer = []
    for letter in cigar:
        if letter.isdigit():
            count += letter #string addition
        elif letter in "MIDNSHP=X":
            answer.append((int(count), letter))
            count = ""
        else:
            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
    return answer

assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [(14,'S'),(15,'M'),(1,'P'),(1,'D'),(3,'P'),(54,'M'),(1,'D'),(34,'M'),(5,'S')]

def count_cigar(cigar):
    answer = dict((operator,0) for operator in "MIDNSHP=X")
    for count, operator in decode_cigar(cigar):
        answer[operator] += count
    return answer

def cigar_seq_len(cigar):
    len = 0
    for count, operator in decode_cigar(cigar):
        if operator in "MIS=X":
            len += count
    return len

def make_ungapped_ref_cigar_m(contig, read):
    #For testing legacy code which expects CIGAR with M rather than X/=
    assert len(contig) == len(read)
    cigar = ""
    count = 0
    mode = ""
    for c,r in zip(contig, read):
        if c == "*" and r == "*":
            if mode!="P":
                if count: cigar += "%i%s" % (count, mode)
                mode = "P"
                count = 1
            else:
                count+=1
        elif c != "*" and r != "*":
            #alignment match/mismatch
            if mode!="M":
                if count: cigar += "%i%s" % (count, mode)
                mode = "M"
                count = 1
            else:
                count+=1
        elif c == "*":
            if mode!="I":
                if count: cigar += "%i%s" % (count, mode)
                mode = "I"
                count = 1
            else:
                count+=1
        elif r == "*":
            if mode!="D":
                if count: cigar += "%i%s" % (count, mode)
                mode = "D"
                count = 1
            else:
                count+=1
        else:
            assert False
    if count: cigar += "%i%s" % (count, mode)
    if len(read.replace("*", "")) != cigar_seq_len(cigar):
        raise ValueError("%s versus %i, %s" % (cigar, len(read.replace("*", "")), read))
    return cigar

def make_ungapped_ref_cigar(contig, read):
    #WARNING - This function expects contig and read to be in same case!
    assert len(contig) == len(read)
    cigar = ""
    count = 0
    mode = "" #Character codes in CIGAR string
    for c,r in zip(contig, read):
        if c == "*" and r == "*":
            if mode!="P":
                if count: cigar += "%i%s" % (count, mode)
                mode = "P"
                count = 1
            else:
                count+=1
        elif c != "*" and r != "*":
            #alignment match/mismatch
            #CIGAR in SAM v1.2 just had M for match/mismatch
            if c==r:
                if mode!="=":
                    if count: cigar += "%i%s" % (count, mode)
                    mode = "="
                    count = 1
                else:
                    count+=1
            else:
                if mode!="X":
                    if count: cigar += "%i%s" % (count, mode)
                    mode = "X"
                    count = 1
                else:
                    count+=1
        elif c == "*":
            if mode!="I":
                if count: cigar += "%i%s" % (count, mode)
                mode = "I"
                count = 1
            else:
                count+=1
        elif r == "*":
            if mode!="D":
                if count: cigar += "%i%s" % (count, mode)
                mode = "D"
                count = 1
            else:
                count+=1
        else:
            assert False
    if count: cigar += "%i%s" % (count, mode)
    if len(read.replace("*", "")) != cigar_seq_len(cigar):
        raise ValueError("%s versus %i, %s" % (cigar, len(read.replace("*", "")), read))
    return cigar

    
assert make_ungapped_ref_cigar("ACGTA" ,"ACGTA") == "5="
assert make_ungapped_ref_cigar("ACGTA" ,"CGTAT") == "5X"
assert make_ungapped_ref_cigar("ACGTA" ,"ACTTA") == "2=1X2="
assert make_ungapped_ref_cigar("ACG*A" ,"ACT*A") == "2=1X1P1="
assert make_ungapped_ref_cigar("ACGTA" ,"ACT*A") == "2=1X1D1="
assert make_ungapped_ref_cigar("ACG*A" ,"ACTTA") == "2=1X1I1="

def make_gapped_ref_cigar_m(contig, read):
    #For testing legacy code which expects CIGAR with M rather than X/= 
    #WARNING - This function expects contig and read to be in same case!
    assert len(contig) == len(read)
    cigar = ""
    count = 0
    d_count = 0
    mode = "" #Character codes in CIGAR string
    for c,r in zip(contig, read):
        if r == "*":
            if mode!="D":
                if count: cigar += "%i%s" % (count, mode)
                mode = "D"
                count = 1
            else:
                count+=1
        elif r != "*":
            #alignment match/mismatch
            if mode!="M":
                if count: cigar += "%i%s" % (count, mode)
                mode = "M"
                count = 1
            else:
                count+=1
    if count: cigar += "%i%s" % (count, mode)
    return cigar


def make_gapped_ref_cigar(contig, read):
    #WARNING - This function expects contig and read to be in same case!
    assert len(contig) == len(read)
    cigar = ""
    count = 0
    d_count = 0
    mode = "" #Character codes in CIGAR string
    for c,r in zip(contig, read):
        if r == "*":
            if mode!="D":
                if count: cigar += "%i%s" % (count, mode)
                mode = "D"
                count = 1
            else:
                count+=1
        elif r != "*":
            #alignment match/mismatch
            #CIGAR in SAM v1.2 just had M for match/mismatch
            if c==r:
                if mode!="=":
                    if count: cigar += "%i%s" % (count, mode)
                    mode = "="
                    count = 1
                else:
                    count+=1
            else:
                if mode!="X":
                    if count: cigar += "%i%s" % (count, mode)
                    mode = "X"
                    count = 1
                else:
                    count+=1
    if count: cigar += "%i%s" % (count, mode)
    return cigar

assert make_gapped_ref_cigar("ACGTA" ,"CGTAT") == "5X"
assert make_gapped_ref_cigar("ACGTA" ,"ACGTA") == "5="
assert make_gapped_ref_cigar("AC*TA" ,"ACGTA") == "2=1X2="
assert make_gapped_ref_cigar("AC*TA" ,"AC*TA") == "2=1D2="
assert make_gapped_ref_cigar("ACGTA" ,"AC*TA") == "2=1D2="
assert make_gapped_ref_cigar("ACGTA" ,"ACTTA") == "2=1X2="
assert make_gapped_ref_cigar("ACG*A" ,"ACT*A") == "2=1X1D1="
assert make_gapped_ref_cigar("ACGTA" ,"ACT*A") == "2=1X1D1="
assert make_gapped_ref_cigar("ACG*A" ,"ACTTA") == "2=2X1="

if gapped_sam:
    if CIGAR_M:
        make_cigar = make_gapped_ref_cigar_m
    else:
        make_cigar = make_gapped_ref_cigar
else:
    if CIGAR_M:
        make_cigar = make_ungapped_ref_cigar_m
    else:
        make_cigar = make_ungapped_ref_cigar


contig_lines_to_ignore = ['NR', #number of reads
                          'LC', #padded contig length
                          'CT', #consensus tag
                         ]
re_contig_lines_to_ignore = re.compile(r'^(%s)\t' % '|'.join(contig_lines_to_ignore))
read_lines_to_ignore = ['SV', #sequencing vector
                        'LR', #read length
                        'TF', #min estimated template length
                        'TT', #max estimated template length
                        'SF', #sequencing file
                        'AO', #align to original
                        'RT', #reads tag
                        'MT', #machine type 
                        'IB', #backbone
                        'IC', #coverage equivalent
                        'IR', #rail
                        'BC', #not sure what this is yet...
                        ]
re_read_lines_to_ignore = re.compile(r'^(%s)\t' % '|'.join(read_lines_to_ignore))
assert re_read_lines_to_ignore.match('LR\t2000\n')
assert re_read_lines_to_ignore.match('TF\t2000\n')
assert re_read_lines_to_ignore.match('TT\t5000\n')
assert not re_read_lines_to_ignore.match('LN\tFred\n')

cached_pairs = dict()
maf_handle = open(maf)
while True:
    line = maf_handle.readline()
    if not line: break
    assert line.startswith("CO\t"), line
    
    contig_name = line.rstrip().split("\t")[1]
    assert contig_name in ref_lens
    padded_con_seq = None
    padded_con_qual = None
    current_read = None
    mapping = None
    ct_tags = []

    while True:
        line = maf_handle.readline()
        if line == "EC\n":
            break
        elif line.startswith("CS\t"):
            padded_con_seq = line.rstrip().split("\t")[1].upper()
            if gapped_sam:
                mapping = None
                assert ref_lens[contig_name] == len(padded_con_seq), \
                    "Gapped reference length mismatch for %s" % contig_name
                assert ref_md5[contig_name] == seq_md5(padded_con_seq), \
                    "Gapped reference checksum mismatch for %s" % contig_name
                cigar = make_cigar(padded_con_seq, padded_con_seq)
            else:
                mapping = []
                index = 0
                for letter in padded_con_seq:
                    mapping.append(index)
                    if letter != "*":
                        index+=1
                assert ref_lens[contig_name] == len(padded_con_seq) - padded_con_seq.count("*"), \
                    "Ungapped reference length mismatch for %s" % contig_name
                assert ref_md5[contig_name] == seq_md5(padded_con_seq.replace("*","")), \
                    "Ungapped reference checksum mismatch for %s" % contig_name
                if CIGAR_M:
                    cigar = "%iM" % ref_lens[contig_name]
                else:
                    cigar = "%i=" % ref_lens[contig_name]
            #Record a dummy read for the contig sequence, FLAG = 516
            #print "%s\t516\t%s\t1\t0\t%s\t*\t0\t0\t%s\t*" \
            #      % (contig_name, contig_name, cigar, padded_con_seq.replace("*",""))
            for start, end, tag, text in ct_tags:
                #Had to wait for the CS line to see the padded reference
                flag = 768 #(filtered and secondary for annotation dummy reads)
                if start > end:
                    #Reverse strand
                    flag += 0x10
                    start, end = end, start 
                #Will potentially want P in the CIGAR string
                s = padded_con_seq[start-1:end]
                cigar = make_cigar(s, s)
                #At the time of writing, the samtools spec in SVN uses * for
                #the sequence of these dummy reads. However, that seems to
                #upset some viewers, e.g. StringIndexOutOfBoundsException from
                #Tablet - even with a proper sequence this is still not working
                #in IGV for me, it gives ArrayIndexOutOfBoundsException messages
                #in its log file.
                s = s.replace("*", "")
                if not gapped_sam:
                    assert mapping is not None and len(mapping) == len(padded_con_seq)
                    start = mapping[start-1] + 1 #SAM and MIRA one based
                print "*\t%i\t%s\t%i\t255\t%s\t*\t0\t0\t%s\t*\tRT:Z:%s|%s" \
                      % (flag, contig_name, start, cigar, s, tag, text)
                del s, cigar
                #Note where the CT tag described just an insert in the reference,
                #the dummy read "sequence" length is zero. The CIGAR string
                #will be just inserts (for a traditional unpadded reference,
                #which will be interesting for a viewer to display!) or just
                #just pads (for a padded reference, which Tablet shows fine
                #as a read made up of gaps only).
            ct_tags = []
        elif line.startswith("CQ\t"):
            assert len(padded_con_seq) == len(line.rstrip().split("\t")[1])
        elif line.startswith("CT\t") and RECORD_CT:
            try:
                #MIRA usually uses tabs, but I've seen spaces sometimes here:
                tag, start, end, text = line[3:].strip().split(None,3)
            except ValueError:
                try:
                    #See if there was no comment string
                    tag, start, end = line[3:].strip().split(None,2)
                    text = ""
                except ValueError:
                    raise ValueError("Problem with %r" % line)
            start = int(start)
            end = int(end)
            text = text.replace("\t", "%09").replace("|", "%A6").strip()
            ct_tags.append((start, end, tag, text))
        elif line == "\\\\\n":
            while line != "//\n":
                current_read = Read(contig_name)
                while True:
                    line = maf_handle.readline()
                    if line == "//\n":
                        break
                    elif line.startswith("RD\t"):
                        current_read.read_name = line.rstrip().split("\t")[1]
                        assert current_read.read_name
                    elif line.startswith("RS\t"):
                        current_read.read_seq = line.rstrip().split("\t")[1].upper()
                    elif line.startswith("RQ\t"):
                        current_read.read_qual = line.rstrip().split("\t")[1]
                        assert len(current_read.read_qual) == len(current_read.read_seq)
                    elif line.startswith("TN\t"):
                        current_read.template_name = line.rstrip().split("\t")[1]
                        #assert current_read.read_name.startswith(current_read.template_name)
                    elif line.startswith("DI\t"):
                        if line == "DI\tF\n":
                            current_read.first_in_pair = True
                        elif line == "DI\tR\n":
                            current_read.first_in_pair = False
                        else:
                            raise ValueError(line)
                    elif line.startswith("SL\t"):
                        current_read.vect_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("SR\t"):
                        current_read.vect_right = int(line.rstrip().split("\t")[1])
                    elif line.startswith("QL\t"):
                        current_read.qual_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("QR\t"):
                        current_read.qual_right = int(line.rstrip().split("\t")[1])
                    elif line.startswith("CL\t"):
                        current_read.clip_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("CR\t"):
                        current_read.clip_right = int(line.rstrip().split("\t")[1])
                    elif line.startswith("RT\t"):
                        #Read tag, will turn into part of a SAM PT tag
                        try:
                            tag, start, end, text = line[3:].strip().split(None,3)
                        except ValueError:
                            tag, start, end = line[3:].strip().split(None,2)
                            text = ""
                        start = int(start)
                        end = int(end)
                        #Must now convert from these padded reference coords
                        #to padded read coords... need the reads mapping
                        #position from the AT tags to do this.
                        text = text.replace("\t", "%09").replace("|", "%A6").strip()
                        current_read.annotations.append((start, end, tag, text))
                    elif line.startswith("ST\t"):
                        current_read.seq_tech = line.rstrip().split("\t")[1]
                    elif line.startswith("SN\t"):
                        current_read.strain =line.rstrip().split("\t")[1]
                    elif line == "ER\n":
                        #End of read - next line should be AT then //
                        pass
                    elif line.startswith("AT\t"):
                        #Assembles to
                        x1, y1, x2, y2 = [int(i) for i in line[3:-1].split()]
                        #AT Four integers: x1 y1 x2 y2
                        #The AT (Assemble_To) line defines the placement of the
                        #read in the contig and follows immediately the closing
                        #"ER" of a read so that parsers do not need to perform
                        #time consuming string lookups. Every read in a contig
                        #has exactly one AT line.
                        #The interval [x2 y2] of the read (i.e., the unclipped
                        #data, also called the 'clear range') aligns with the
                        #interval [x1 y1] of the contig. If x1 > y1 (the contig
                        #positions), then the reverse complement of the read is
                        #aligned to the contig. For the read positions, x2 is
                        #always < y2.
                        if x1 > y1:
                            current_read.ref_rc = True
                            #SAM stores these backwards:
                            cigar = make_cigar(padded_con_seq[y1-1:x1],
                                               reverse_complement(current_read.read_seq[x2-1:y2]))
                            #cigar = "%iM" % (x1-y1+1)
                            if x2 > 1:
                                cigar += "%iS" % (x2-1)
                            if y2 < len(current_read.read_seq):
                                cigar = "%iS%s" % (len(current_read.read_seq)-y2, cigar)
                        else:
                            cigar = make_cigar(padded_con_seq[x1-1:y1],
                                               current_read.read_seq[x2-1:y2])
                            #cigar = "%iM" % (y1-x1+1)
                            if x2 > 1:
                                cigar = "%iS%s" % (x2-1, cigar)
                            if y2 < len(current_read.read_seq):
                                cigar += "%iS" % (len(current_read.read_seq)-y2)
                        current_read.cigar = cigar
                        current_read.padded_pos = min(x1, y1)-1 #zero based
                        if gapped_sam:
                            current_read.ref_pos = current_read.padded_pos+1 #one based for SAM
                        else:
                            current_read.ref_pos = mapping[current_read.padded_pos]+1 #one based for SAM
                        break #End of this read
                    elif not line:
                        raise ValueError("EOF in read")
                    elif re_read_lines_to_ignore.match(line):
                        pass
                    else:
                        sys.stderr.write("Bad line in read: %s\n" % repr(line))
                        #Continue and hope we can just ignore it!
                if line == "//\n":
                    break
                if current_read.need_partner():
                    cached_pairs[(current_read.template_name, current_read.first_in_pair)] = current_read
                elif current_read.is_paired():
                    cached_pairs[(current_read.template_name, current_read.first_in_pair)] = current_read
                    print current_read.get_partner()
                    print current_read
                    #Clear from cache
                    del cached_pairs[(current_read.template_name, current_read.first_in_pair)]
                    del cached_pairs[(current_read.template_name, not current_read.first_in_pair)]
                else:
                    print current_read
        elif not line:
            raise ValueError("EOF in contig")
        elif re_contig_lines_to_ignore.match(line):
            pass
        else:
            sys.stderr.write("Bad line in contig: %s" % repr(line))
            #Continue and hope we can just ignore it!
    #print contig_name, ref_lens[contig_name]

if cached_pairs:
    log("Almost done, %i orphaned paired reads remain" % len(cached_pairs))
    #Special cases - paired reads where partner was not found in MAF file
    for current_read in cached_pairs.itervalues():
        print current_read
else:
    log("No orphaned paired reads")
log("Done")
