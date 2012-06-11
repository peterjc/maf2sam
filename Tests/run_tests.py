#!/usr/bin/env python
import os
import sys
import itertools
import subprocess
from StringIO import StringIO

def update_maf2sam(maf, fasta, sam):
    assert os.path.isfile(maf)
    assert os.path.isfile(fasta)
    cmd = "../maf2sam.py %s %s > %s" % (f, m, s)
    print cmd
    return_code = os.system(cmd)
    if return_code:
        sys.stderr.write("Return code %i" % return_code)
        sys.exit(return_code)

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
print "Done"
