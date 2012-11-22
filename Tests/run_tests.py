#!/usr/bin/env python
import os
import sys
import itertools
import subprocess
from StringIO import StringIO

def update_maf2sam(maf, fasta, sam):
    assert os.path.isfile(maf)
    assert os.path.isfile(fasta)
    cmd1 = "../maf2sam.py %s %s > %s" % (f, m, s)
    cmd2 = "../sam2bam.py %s" % s
    for cmd in [cmd1, cmd2]:
        print cmd
        return_code = os.system(cmd)
        if return_code:
            sys.stderr.write("Return code %i" % return_code)
            sys.exit(return_code)

def test_depad_bam(fasta, padded_bam, unpadded_bam):
    #Compare the BAM files as both sorted
    cmd = ["/Users/pjcock/repositories/samtools/samtools", "depad",
           "-s", "-T", fasta, padded_bam]
    child = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    depad_out, stderr = child.communicate()
    return_code = child.returncode
    del child
    if return_code:
        print "FAILED, depad return code %i for:\n%s" % (return_code, " ".join(cmd))
        print stderr
        sys.exit(1)
    cmd = ["samtools", "view", "-h", unpadded_bam]
    child = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    maf2sam_out, stderr = child.communicate()
    return_code = child.returncode
    del child
    if return_code:
        print "FAILED, samtools view return code %i for:\n%s" % (return_code, " ".join(cmd))
        print stderr
        sys.exit(1)
    failed = False
    for old, new in itertools.izip(maf2sam_out.split("\n"), depad_out.split("\n")):
            if old.startswith("@") and new.startswith("@"):
                #Not checking header details
                continue
            old_id = old.split("\t",1)[0]
            new_id = new.split("\t",1)[0]
            if old_id != new_id:
                print "Sort problem? maf2sam %s vs %s from samtools depad" % (old_id, new_id)
                failed = True
            elif old != new:
                print
                print "maf2sam:"
                print repr(old)
                print
                print "samtools depad:"
                print repr(new)
                failed = True
                break
    #if failed:
    #    print
    #    print "FAILED, did not reproduce %s" % s
    #    sys.exit(1)


def test_maf2sam(maf, fasta, sam):
    assert os.path.isfile(maf)
    assert os.path.isfile(fasta)
    assert os.path.isfile(sam)
    print "../maf2sam.py %s %s > %s" % (f, m, s)

    child = subprocess.Popen(["../maf2sam.py", f, m],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    stdout, stderr = child.communicate()
    return_code = child.returncode
    del child
    if return_code:
        print "FAILED, return code %i" % return_code
        print stderr
        sys.exit(1)

    with open(sam, "rU") as wanted:
        for old, new in itertools.izip(wanted, stdout.split("\n")):
            if old != new+"\n":
                print
                print "Expected:"
                print repr(old)
                print
                print "Produced:"
                print repr(new)
                print
                print "FAILED, did not reproduce %s" % s
                sys.exit(1)

if "-g" in sys.argv:
    update = True
    print "Updating..."
else:
    print "Testing..."
    update = False

for d in os.listdir("."):
    if not os.path.isdir(d):
        continue
    for m in os.listdir(d):
        if not m.endswith(".maf"):
            continue
        m = os.path.join(d, m)
        for ref in ["padded", "unpadded"]:
            f = "%s.%s.fasta" % (m[:-4], ref)
            s = "%s.%s.sam" % (m[:-4], ref)
            if update:
                update_maf2sam(m, f, s)
                continue
            if not os.path.isfile(f):
                print "Missing %s" % f
                continue
            if not os.path.isfile(s):
                print "Missing %s" % s
                continue
            test_maf2sam(m, f, s)
        test_depad_bam("%s.padded.fasta" % m[:-4],
                       "%s.padded.bam" % m[:-4],
                       "%s.unpadded.bam" % m[:-4])
print "Done"
