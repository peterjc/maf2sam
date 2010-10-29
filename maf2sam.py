#v000 - Uses the gapped co-ordindates
#v001 - Use the ungapped co-ordindates
from Bio.Seq import reverse_complement
from Bio import SeqIO

ref = "AA14X_out.unpadded.fasta"
maf = "AA14X_out.maf"

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

maf_handle = open(maf)
while True:
    line = maf_handle.readline()
    if not line: break
    assert line.startswith("CO\t")
    
    contig_name = line.rstrip().split("\t")[1]
    assert contig_name in ref_lens
    padded_con_seq = None
    padded_con_qual = None

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
                read_name, template_name, read_seq, forward = "", "", "", True
                flag, ref_pos, map_qual = 0, 0, 255
                cigar = ""
                mate_ref_name = "="
                mate_ref_pos = 0
                insert_size = 0
                vect_left, vect_right = 0, 0
                qual_left, qual_right = 0, 0
                clip_left, clip_right = 0, 0
                tags=""
                while True:
                    line = maf_handle.readline()
                    if line == "//\n":
                        break
                    elif line.startswith("RD\t"):
                        read_name = line.rstrip().split("\t")[1]
                        assert read_name
                    elif line.startswith("LR\t"):
                        #Optional length of read
                        pass
                    elif line.startswith("RS\t"):
                        read_seq = line.rstrip().split("\t")[1]
                    elif line.startswith("RQ\t"):
                        read_qual = line.rstrip().split("\t")[1]
                        assert len(read_qual) == len(read_seq)
                    elif line.startswith("SF\t"):
                        #name of seq file
                        pass
                    elif line.startswith("TN\t"):
                        template_name = line.rstrip().split("\t")[1]
                        assert read_name.startswith(template_name)
                    elif line.startswith("DI\t"):
                        if line == "DI\tF\n":
                            forward = True
                        elif line == "DI\tR\n":
                            forward = False
                        else:
                            raise ValueError(line)
                    elif line.startswith("SL\t"):
                        vect_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("SR\t"):
                        vect_right = int(line.rstrip().split("\t")[1])
                    elif line.startswith("QL\t"):
                        qual_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("QR\t"):
                        qual_right = int(line.rstrip().split("\t")[1])
                    elif line.startswith("CL\t"):
                        clip_left = int(line.rstrip().split("\t")[1])
                    elif line.startswith("CR\t"):
                        clip_right = int(line.rstrip().split("\t")[1])
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
                            flag += 0x10 #maps to reverse strand
                            #SAM stores these backwards:
                            cigar = make_cigar(padded_con_seq[y1-1:x1],
                                               reverse_complement(read_seq[x2-1:y2]))
                            read_qual = read_qual[::-1]
                            read_seq = reverse_complement(read_seq)
                            #cigar = "%iM" % (x1-y1+1)
                            if x2 > 1:
                                cigar += "%iS" % (x2-1)
                            if y2 < len(read_seq):
                                cigar = "%iS%s" % (len(read_seq)-y2, cigar)
                        else:
                            cigar = make_cigar(padded_con_seq[x1-1:y1],
                                               read_seq[x2-1:y2])
                            #cigar = "%iM" % (y1-x1+1)
                            if x2 > 1:
                                cigar = "%iS%s" % (x2-1, cigar)
                            if y2 < len(read_seq):
                                cigar += "%iS" % (len(read_seq)-y2)
                        padded_pos = min(x1, y1)-1 #zero based
                        ref_pos = mapping[padded_pos]+1 #one based for SAM
                        break #End of this read
                    elif not line:
                        raise ValueError("EOF in read")
                    else:
                        raise ValueError("Bad line in read: %s" % repr(line))
                if line == "//\n":
                    break
                if not template_name:
                    assert read_name
                    template_name = read_name
                elif template_name != read_name:
                    #Looks like a paired end read!
                    flag += 1
                    if forward:
                        flag += 0x40 #forward partner
                    else:
                        flag += 0x80 #reverse partner
                assert not tags
                read_seq_unpadded = read_seq.replace("*", "")
                read_qual_unpadded = "".join(q for (l,q) in zip(read_seq,read_qual) if l!="*")
                if "D" not in cigar:
                    #Need deletions offset
                    if len(read_seq_unpadded) != sum(int(x) for x in cigar.replace("I","M").replace("S","M").split("M") if x):
                        raise ValueError("%s vs %i for %s" % (cigar, len(read_seq_unpadded), read_seq))
                assert len(read_seq_unpadded) == len(read_qual_unpadded)
                print "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\t%s\t%s" % \
                    (template_name, flag, contig_name, ref_pos, map_qual, cigar,
                     mate_ref_name, mate_ref_pos, insert_size,
                     read_seq_unpadded, read_qual_unpadded)
        elif not line:
            raise ValueError("EOF in contig")
        else:
            raise ValueError("Bad line in contig: %s" % repr(line))
    #print contig_name, ref_lens[contig_name]
 
 