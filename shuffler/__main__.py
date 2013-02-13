"""\
A command-line interface to compare an observed relationship between 2 sets of
intervals to the distribution of that relationship derived from a number of
shufflings. This "relationship" can be anything that returns a single number,
e.g. total bases of overlap between the 2 sets, or number of overlaps.
"""
import argparse
import sys
import tempfile

shuffling_constraints_doc = """
    each of these specify some constraint on where to send the randomized
    intervals. any of them can be padded using a syntax like:
        -a query.bed:1000
    which will add 1000 bases to either end of each randomly placed interval.
        -b subject.bed:-1000:25
    will extend the interval by 1000 bases upstream and 25 bases downstream.
    if a strand (+/-) designation is present in the last column of the file,
    it will be used to determine up/dowstream.
"""

import random
from toolshed import reader, nopen
from .files import stream_file
from .shuffler import Shuffler, jaccard_value, jaccard_length \
        , num_intersections

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    g_inputs = p.add_argument_group("input intervals")
    g_inputs.add_argument("-a", metavar="SHUFFLED",
            help="bed file (the -a file is always shuffled)")
    g_inputs.add_argument("-b", metavar="UNSHUFFLED",
            help="bed file (the -b file is never shuffled)")

    g_constraints = p.add_argument_group("shuffling constraints",
            shuffling_constraints_doc)
    g_constraints.add_argument("-g", dest="genome", help=
        "genome version (to get chromosomes), e.g. mm8, dm3, hg19 or a file")
    # can have --domain or --domains not both
    g_doms = g_constraints.add_mutually_exclusive_group()
    g_doms.add_argument("--domain", "--include", dest="domain",
            help="(optional) shuffle -a intervals inside this domain. "
            "(may be specified multiple times)")
    g_doms.add_argument("--domains", help=
            "like domain, but the 4th column is used to split into separate "
            "domains a value is returned for each of the sub-domains")
    #
    g_constraints.add_argument("--exclude", help=
        "(optional) do not allow shuffled intervals to land in the "
        "intervals specified in this file (may be specified "
        "multiple times)")


    # misc
    p.add_argument("-n", dest="n", help="(optional) number of times to shuffle", type=int,
            default=1000)
    p.add_argument("-t", dest="threads", help="(optional) number of threads to use",
            type=int, default=1)
    p.add_argument("--seed", help="(optional) seed for random "
            "number generator", type=int, default=random.seed())
    p.add_argument("--metric", help="metric by which to evaluate overlap",
            default=jaccard_value)

    args = p.parse_args()
    if (args.a is None or args.b is None):
        sys.exit(not p.print_help())

    shuffle(args)

def tofile(fiter, fname):
    fh = nopen(fname, "w")
    for line in fiter:
        print >>fh, line
    fh.close()
    return fname

class _wrapper_fn(object):
    def __init__(self, command_string):
        self.command_string = command_string
        self.func_name = command_string
    def __call__(self, fh):
        out = open(tempfile.mktemp(dir="/tmp/"), "w")
        print out.name
        for row in fh:
            out.write(row)
        out.close()
        value = nopen("%s < %s" % (self.command_string, out.name)).next()
        return int(value)

def shuffle(args):

    a = tofile(stream_file(args.a), tempfile.mktemp(dir="/tmp"))
    b = tofile(stream_file(args.b), tempfile.mktemp(dir="/tmp"))

    value_fn = args.metric
    command = "bedtools jaccard -a %(query)s -b %(subject)s"
    if isinstance(value_fn, basestring):
        try:
            if value_fn == "num_intersections":
                command = "bedtools intersect -a %(query)s -b %(subject)s | wc -l"
            value_fn = globals()[value_fn]
        except KeyError:
            command_string = "|" + value_fn.lstrip("|")
            value_fn = _wrapper_fn(command_string)
            command = "bedtools intersect -a %(query)s -b %(subject)s"

    s = Shuffler(a, b, args.genome, value_fn, n=args.n,
                   seed=args.seed, map=args.threads if args.threads > None else map)
    print s.run(command=command, sims=True)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()

