"""diff

Prints out where two hic files are different

First attempt; resolutions, chromosomes, norms hard-coded 
(these should instead be taken from respective hic files)
Only checks intra-chromosomal matrices
This is very slow particularly for FRAG since it looks at 
every single matrix entry
"""
import straw
import math
import argparse

bpresolutions=(2500000,1000000,500000,100000,50000,25000,10000,5000)
fragresolutions=(500,200,100,50,20,5,2,1) 
chrs=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')
norms=('NONE', 'KR', 'VC')

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    '''
    Python 2 implementation of Python 3.5 math.isclose()
    https://hg.python.org/cpython/file/tip/Modules/mathmodule.c#l1993
    '''
    # sanity check on the inputs
    if rel_tol < 0 or abs_tol < 0:
        raise ValueError("tolerances must be non-negative")

    # short circuit exact equality -- needed to catch two infinities of
    # the same sign. And perhaps speeds things up a bit sometimes.
    if a == b:
        return True

    # This catches the case of two infinities of opposite sign, or
    # one infinity and one finite number. Two infinities of opposite
    # sign would otherwise have an infinite relative tolerance.
    # Two infinities of the same sign are caught by the equality check
    # above.
    if math.isinf(a) or math.isinf(b):
        return False

    # equality check above would return false for nan, but we want
    # to return true
    if math.isnan(a) and math.isnan(b):
        return True

    # now do the regular computation
    # this is essentially the "weak" test from the Boost library
    diff = math.fabs(b - a)
    result = (((diff <= math.fabs(rel_tol * b)) or
               (diff <= math.fabs(rel_tol * a))) or
              (diff <= abs_tol))
    return result


parser = argparse.ArgumentParser(description='Print differences between hic files, or nothing if they are the same')
parser.add_argument('filenames', metavar='file', nargs=2, help='files to diff')
parser.add_argument('--frag_only', dest='frag_only', default=False, action='store_true', help='only examine BP resolutions')
parser.add_argument('--bp_only', dest='bp_only', default=False, action='store_true', help='only examine frag resolutions')
args = parser.parse_args()
file1=args.filenames[0]
file2=args.filenames[1]

if (not args.frag_only):
    for res in bpresolutions:
        for chr in chrs:
            for norm in norms:
                result=straw.straw(norm, file1, chr, chr, "BP", res)
                result2=straw.straw(norm, file2, chr, chr, "BP", res)
                for i in range(len(result[0])):
                    if (not isclose(result[0][i],result2[0][i]) or not isclose(result[1][i],result2[1][i]) or not isclose(result[2][i], result2[2][i])):
                        print("Difference at {0} {1} BP {2}".format(norm, chr, res))
                        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(result[0][i], result[1][i], result[2][i], result2[0][i], result2[1][i], result2[2][i]))

if (not args.bp_only):
    for res in fragresolutions:
        for chr in chrs:
            for norm in norms:
                result=straw.straw(norm, file1, chr, chr, "FRAG", res)
                result2=straw.straw(norm, file2, chr, chr, "FRAG", res)
                for i in range(len(result[0])):
                    if (not isclose(result[0][i],result2[0][i]) or not isclose(result[1][i],result2[1][i]) or not isclose(result[2][i], result2[2][i])):
                        print("Difference at {0} {1} FRAG {2}".format(norm, chr, res))
                        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(result[0][i], result[1][i], result[2][i], result2[0][i], result2[1][i], result2[2][i]))
