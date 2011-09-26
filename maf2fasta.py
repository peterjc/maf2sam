#!/usr/bin/env python
#Quick hack to turn MIRA's assembly format (MAF)
#into FASTA contig files. Use at your own risk.
import sys
import os

try:
   maf_file, padded_file, unpadded_file = sys.argv[1:]
except:
   print "Expects three filename arguments: MAF (input), padded FASTA (outout), unpadded FASTA (outut)"
   sys.exit(1)

handle = open(maf_file)
padded_handle = open(padded_file, "w")
unpadded_handle = open(unpadded_file, "w")

for line in handle:
    if line.startswith("CO\t"):
        contig = line[3:].strip()
        print contig
        padded_handle.write(">%s\n" % contig)
        unpadded_handle.write(">%s\n" % contig)
    elif line.startswith("CS\t"):
        seq = line[3:].strip()
        for i in xrange(0, len(seq), 80):
            padded_handle.write(seq[i:i+80] + "\n")
        seq = seq.replace("*", "").replace("-", "")
        for i in xrange(0, len(seq), 80):
           unpadded_handle.write(seq[i:i+80] + "\n")

handle.close()
padded_handle.close()
unpadded_handle.close()

print "Done"
