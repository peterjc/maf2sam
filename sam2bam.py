#!/usr/bin/env python
"""Quick script to convert SAM file(s) into BAM, sort and index it.

Just makes a series of calls to samtools (view, index, and idxstats).

This is not really intended for public use, but has been helpful
in my testing for MAF to SAM/BAM conversion.
"""
import os
import sys

def run(cmd):
    print cmd
    assert not os.system(cmd)

for sam in sys.argv[1:]:
    assert os.path.isfile(sam)
    assert sam.endswith(".sam")
    prefix = sam[:-4]
    bam = prefix + ".bam"
    if os.path.isfile(bam):
        os.remove(bam)
    if os.path.isfile(bam+".bai"):
        os.remove(bam+".bai")
    run("samtools view -b -S %s | samtools sort - %s" % (sam, prefix))
    assert os.path.isfile(bam)
    run("samtools index %s" % bam)
    assert os.path.isfile(bam+".bai")
    run("samtools idxstats %s" % bam)
