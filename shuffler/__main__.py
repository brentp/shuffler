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
from .files import stream_file, parse_file_pad
from .shuffler import Shuffler, jaccard_values
import atexit
import os

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    g_inputs = p.add_argument_group("input intervals")
    g_inputs.add_argument("-a", metavar="SHUFFLED",
            help="bed file (the -a file is always shuffled)")

    bs = g_inputs.add_mutually_exclusive_group()
    bs.add_argument("-b", metavar="UNSHUFFLED",
            help="bed file (the -b file is never shuffled)")
    bs.add_argument("--bs",
            help="single bed file to be split on values in 4th column")

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
    g_constraints.add_argument("--exclude", "-excl", action="append", help=
        "(optional) do not allow shuffled intervals to land in the "
        "intervals specified in this file (may be specified "
        "multiple times)", default=[])

    # misc
    p.add_argument("-n", dest="n", help="(optional) number of times to shuffle", type=int,
            default=100)
    p.add_argument("-t", dest="threads", help="(optional) number of threads to use",
            type=int, default=1)
    p.add_argument("--seed", help="(optional) seed for random "
            "number generator", type=int, default=random.seed())
    p.add_argument("--metric", help="metric by which to evaluate overlap",
            default=jaccard_values)

    p.add_argument("--png", help="save a png of the distributions from the sims")

    args = p.parse_args()
    if (args.a is None or (args.b is None and args.bs == [])):
        sys.exit(not p.print_help())

    shuffle(args)

def tofile(fiter, fname):
    fh = nopen(fname, "w")
    for line in fiter:
        print >>fh, line.rstrip("\r\n")
    fh.close()
    atexit.register(os.unlink, fname)
    return fname


class _wrapper_fn(object):
    def __init__(self, command_string):
        self.command_string = command_string
        self.func_name = command_string
    def __call__(self, fh):
        out = tofile(fh, tempfile.mktemp(dir="/tmp/"))
        try:
            value = nopen("%s < %s" % (self.command_string, out)).next()
            return dict(value=float(value))
        except:
            print self.command_string
            raise

    def __str__(self):
        return "<user defined function: '%s'>" % \
                self.command_string.lstrip("|")

def plot(res, png):
    from matplotlib import pyplot as plt
    if isinstance(res, dict):
        res = [res]
    print res[0].keys()
    plot_keys = [k for k in res[0].keys() if hasattr(res[0][k], "__iter__") and "p_sims_gt" in res[0][k]]

    f, axarr = plt.subplots(len(res), len(plot_keys))
    f.set_size_inches((8, 12))

    for i, row in enumerate(res):
        for j, metric in enumerate(plot_keys):
            ax = axarr[i, j]
            ax = Shuffler.plot(row[metric], ax=ax)
            if ax.is_first_col():
                ax.set_ylabel(row['b'], fontsize=9, rotation="horizontal")
            ax.set_xticks([])
            ax.set_yticks([])
            if ax.is_last_row():
                ax.set_xlabel(metric, fontsize=9, rotation="horizontal")

    f.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, left=0.16, right=0.98)
    f.savefig(png)


def merge_excl(excl_list):
    if len(excl_list) == 1:
        excl = excl_list[0]
    else:
        excl = tempfile.mktemp(dir="/tmp")
        list(nopen("|cat %s | sort -k1,1 -k2,2n | bedtools merge -i - > %s" \
                % (" ".join(excl_list), excl)))
        atexit.register(os.unlink, excl)
    return "-excl %s" % excl

def gen_files(fname, col=3):
    files = {}
    for toks in reader(fname, header=False):
        key = toks[col]
        if not key in files:
            f = tempfile.mktemp(dir="/tmp")
            files[key] = open(f, "w")
        print >>files[key], "\t".join(toks)

    for key, fh in files.iteritems():
        fh.close()
        atexit.register(os.unlink, fh.name)
        yield key, fh.name


def shuffle(args):

    a = tofile(stream_file(args.a), tempfile.mktemp(dir="/tmp")) \
            if ":" in args.a else args.a
    value_fn = args.metric

    if args.b is None:
        binfo = parse_file_pad(args.bs)
        bs = []
        for f in gen_files(binfo['file']):
            binfo['file'] = f
            tmp = tempfile.mktemp(dir="/tmp")
            bs.append(to_file(stream_file(f, binfo), tmp) \
                    if abs(binfo['upstream']) + abs(binfo['downstream']) != 0
                    else f)
    else:

        bs = [(args.b, tofile(stream_file(args.b), tempfile.mktemp(dir="/tmp")) if ":" in
                args.b else args.b)]

    # print a  header
    prefix = "group\t" if args.bs else ""
    if value_fn == jaccard_values:
        print prefix + "\t".join(Shuffler.jaccard_metrics)
    else:
        print prefix + "\t%s" % value_fn


    command = "bedtools jaccard -a %(query)s -b %(subject)s"
    if isinstance(value_fn, basestring):
        if value_fn != jaccard_values:
            command = "bedtools intersect -a %(query)s -b %(subject)s -wo"
            command_string = "|" + value_fn.lstrip("|")
            value_fn = _wrapper_fn(command_string)

    shuffle_str = merge_excl(args.exclude) if args.exclude else ""

    results = []
    for bname, b in bs:
        s = Shuffler(a, b, args.genome, value_fn, n=args.n,
                shuffle_str=shuffle_str,
                       seed=args.seed, map=args.threads if args.threads > 1 else map)

        res = s.run(command=command, sims=True)
        line = ("%s\t" % ( bname)) if args.bs else ""
        if value_fn == jaccard_values:
            print line + "\t".join(("%.4g" % res[metric]['p_sims_gt']) for metric in \
                    Shuffler.jaccard_metrics)
        else:
            print line + "%s\n%.4g" % (value_fn, res['p_sims_gt'])
        res['b'] = bname
        results.append(res)
    if args.png:
        plot(results, args.png)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()

