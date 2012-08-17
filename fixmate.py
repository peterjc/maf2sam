#!/usr/bin/env python
#Quick hack to fix MIRA 3.9.3's SAM output while waiting for a fix.
#The specific error was setting MRNM/RNEXT to * when it should be =.
#
#This script will reset the MRNM and MPOS aka RNEXT and PNEXT fields
#with the information from the mate's POS and RNAME fields.
import sys
import os

if len(sys.argv) != 2:
   print "Expects one filename argument: Input SAM file"
   print
   print "Only supports pairs (not multiple fragments per read)."
   print
   print "Writes SAM to stdout."
   sys.exit(1)
in_file = sys.argv[1]

in_handle = open(in_file)

#Quick and dirty first implementation uses dictionaries in RAM!
mapping1 = dict()
mapping2 = dict()
sys.stderr.write("Starting first pass of %s\n" % in_file)
for line in in_handle:
    if line[0] == "@":
        continue
    qname, flag, rname, pos, rest = line.split("\t",4)
    flag = int(flag)
    if not flag & 0x1:
        #Not a pair
        continue
    elif flag & 0x100:
        #Secondary alignment
        continue
    elif flag & 0x40 and flag & 0x80:
        #Not first or last, but a middle fragment
        sys.stderr.write("Read %s is part of a mult-fragment read (not a simple pair)\n" % qname)
        sys.exit(1)
    elif flag & 0x40:
        #First
        if flag & 0x4:
            #Unmapped, will want to use partner's data for sorting
            pass
        else:
            mapping1[qname] = rname, int(pos)
    elif flag & 0x80:
        #Last
        if flag & 0x4:
            #Unmapped, will want to use partner's data for sorting
            pass
        else:
            mapping2[qname] = rname, int(pos)
    else:
        #A middle fragment, position unknown
        sys.stderr.write("Read %s is part of a mult-fragment read (not a simple pair)\n" % qname)
        sys.exit(1)
sys.stderr.write("Loaded positions of %i first and %i second reads\n" % (len(mapping1), len(mapping2)))
in_handle.seek(0)
sys.stderr.write("Starting second pass of %s\n" % in_file)
for line in in_handle:
    if line[0] == "@":
        sys.stdout.write(line)
        continue
    qname, flag, rname, pos, mapq, cigar, mate_ref, mate_pos, rest = line.split("\t",8)
    flag = int(flag)
    if not flag & 0x1:
        #Not a pair
        sys.stdout.write(line)
        continue
    elif flag & 0x100:
        #Secondary alignment - do anything?
        sys.stdout.write(line)
        continue
    elif flag & 0x40:
        #First
        if flag & 0x8:
            #Mate unmapped
            assert qname not in mapping2, line
            assert flag & 0x8, qname
            mate_ref, mate_pos = "*", 0
        else:
            mate_ref, mate_pos = mapping2[qname]
    else:
        assert flag & 0x80
        if flag & 0x8:
            #Mate unmapped
            mate_ref, mate_pos = mapping1[qname]
            assert qname not in mapping1, line
            assert flag & 0x8, qname
            mate_ref, mate_pos = "*", 0
        else:
            mate_ref, mate_pos = mapping1[qname]
    if flag & 0x4:
        #This read is unmapped - take partner's position for sorting
        ref, pos = mate_ref, mate_pos
    if rname == mate_ref and mate_ref != "*":
        mate_ref = "="
    #TODO - Check TLEN?
    sys.stdout.write("\t".join((qname, str(flag), rname, str(pos), mapq, cigar, mate_ref, str(mate_pos), rest)))

in_handle.close()
sys.stderr.write("Done\n")
