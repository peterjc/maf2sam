"""Simple MIRA alignment format (MAF) to SAM format converter.

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
#v000 - Uses the gapped co-ordindates
#v001 - Use the ungapped co-ordindates
#v002 - Use object for each read
#v003 - Fill in pair partner info
#
#TODO
# - insert size
# - properly paired flag?
# - Record origin read name suffix in tags
# - Record any MIRA annotation in tags?
# - testing!
from Bio.Seq import reverse_complement
from Bio import SeqIO

ref = "AA14X_out.unpadded.fasta"
maf = "AA14X_out.maf"

class Read(object):
    def __init__(self, contig_name, read_name="", template_name="",
                 read_seq="", first_in_pair=True, ref_rc = False,
                 ref_pos=0, map_qual=255, insert_size = 0,
                 vect_left = 0, vect_right = 0,
                 qual_left = 0, qual_right = 0,
                 clip_left = 0, clip_right = 0,
                 tags=""):
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
        self.tags = tags
    
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
        global cached_pairs
        if self.ref_rc:
            flag = 0x10 #maps to reverse strand
            read_seq = reverse_complement(self.read_seq)
            read_qual = self.read_qual[::-1]
        else:
            flag = 0
            read_seq = self.read_seq
            read_qual = self.read_qual
        mate_ref_name = -1
        mate_ref_pos = -1
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

        assert not self.tags
        read_seq_unpadded = read_seq.replace("*", "")
        read_qual_unpadded = "".join(q for (l,q) in zip(read_seq,read_qual) if l!="*")
        cigar = self.cigar
        if "D" not in cigar:
            #Need deletions offset
            if len(read_seq_unpadded) != sum(int(x) for x in cigar.replace("I","M").replace("S","M").split("M") if x):
                raise ValueError("%s vs %i for %s" % (cigar, len(read_seq_unpadded), read_seq))
        assert len(read_seq_unpadded) == len(read_qual_unpadded)
        return "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\t%s\t%s" % \
            (self.template_name, flag, self.contig_name, self.ref_pos,
             self.map_qual, cigar,
             mate_ref_name, mate_ref_pos, self.insert_size,
             read_seq_unpadded, read_qual_unpadded)

print "@HD\tVN:1.0\tSO:unsorted"
print "@CO\tConverted from a MIRA Alignment Format (MAF) file"

ref_lens = {}
for rec in SeqIO.parse(ref, "fasta"):
    assert "*" not in rec.seq
    ref_lens[rec.id] = len(rec)
    print "@SQ\tSN:%s\tLN:%i" % (rec.id, len(rec))

def make_cigar(contig, read):
    assert len(contig) == len(read)
    cigar = ""
    count = 0
    d_count = 0
    mode = ""
    for c,r in zip(contig, read):
        if c == "*" and r == "*":
            pass
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
            d_count+=1
        else:
            assert False
    if count: cigar += "%i%s" % (count, mode)
    if len(read.replace("*", "")) != sum(int(x) for x in cigar.replace("D","M").replace("I","M").split("M") if x) - d_count:
        raise ValueError("%s versus %i, %s" % (cigar, len(read.replace("*", "")), read))
    return cigar
    #return "%iM" % len(read)
    
assert make_cigar("ACGTA" ,"ACGTA") == "5M"
assert make_cigar("ACGTA" ,"ACTTA") == "5M"
assert make_cigar("ACG*A" ,"ACT*A") == "4M"
assert make_cigar("ACGTA" ,"ACT*A") == "3M1D1M", make_cigar("ACGTA" ,"ACT*A")
assert make_cigar("ACG*A" ,"ACTTA") == "3M1I1M", make_cigar("ACG*A" ,"ACTTA")

cached_pairs = dict()
maf_handle = open(maf)
while True:
    line = maf_handle.readline()
    if not line: break
    assert line.startswith("CO\t")
    
    contig_name = line.rstrip().split("\t")[1]
    assert contig_name in ref_lens
    padded_con_seq = None
    padded_con_qual = None
    current_read = None

    while True:
        line = maf_handle.readline()
        if line == "EC\n":
            break
        elif line.startswith("NR\t"):
            #number of reads
            pass
        elif line.startswith("LC\t"):
            #padded contig length
            pass
        elif line.startswith("CT\t"):
            #consensus tag
            pass
        elif line.startswith("CS\t"):
            padded_con_seq = line.rstrip().split("\t")[1]
            assert ref_lens[contig_name] == len(padded_con_seq) - padded_con_seq.count("*")
        elif line.startswith("CQ\t"):
            assert len(padded_con_seq) == len(line.rstrip().split("\t")[1])
        elif line == "\\\\\n":
            assert padded_con_seq
            mapping = []
            index = 0
            for letter in padded_con_seq:
                mapping.append(index)
                if letter != "*":
                    index+=1
                
            while line != "//\n":
                current_read = Read(contig_name)
                while True:
                    line = maf_handle.readline()
                    if line == "//\n":
                        break
                    elif line.startswith("RD\t"):
                        current_read.read_name = line.rstrip().split("\t")[1]
                        assert current_read.read_name
                    elif line.startswith("LR\t"):
                        #Optional length of read
                        pass
                    elif line.startswith("RS\t"):
                        current_read.read_seq = line.rstrip().split("\t")[1]
                    elif line.startswith("RQ\t"):
                        current_read.read_qual = line.rstrip().split("\t")[1]
                        assert len(current_read.read_qual) == len(current_read.read_seq)
                    elif line.startswith("SF\t"):
                        #name of seq file
                        pass
                    elif line.startswith("TN\t"):
                        current_read.template_name = line.rstrip().split("\t")[1]
                        assert current_read.read_name.startswith(current_read.template_name)
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
                    elif line.startswith("AO\t"):
                        #Align to Original... describes indels applied to the raw read
                        pass
                    elif line.startswith("ST\t"):
                        #sequencing technology: Sanger, 454, Solexa, SOLiD
                        pass
                    elif line.startswith("RT\t"):
                        #read tag
                        pass
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
                        current_read.ref_pos = mapping[current_read.padded_pos]+1 #one based for SAM
                        break #End of this read
                    elif not line:
                        raise ValueError("EOF in read")
                    else:
                        raise ValueError("Bad line in read: %s" % repr(line))
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
        else:
            raise ValueError("Bad line in contig: %s" % repr(line))
    #print contig_name, ref_lens[contig_name]

#Special cases - paired reads where partner was not found in ACE file
for current_read in cached_pairs.itervalues():
    print current_read
