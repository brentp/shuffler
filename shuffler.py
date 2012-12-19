from toolshed import reader, nopen
from tempfile import mktemp
from itertools import imap
import random
import os, glob
import sys
import signal

def _run(cmd):
    if cmd.startswith("|"): cmd = cmd[1:]
    list(nopen("|%s" % cmd))

class Shuffler(object):
    """
    given a query and a subject
    compare the original [value] of query and subject with the
    values after n shufflings of the query compared to the same subject.
    """
    def __init__(self, query, subject, genome, value_fn, n=10, shuffle_str="",
            seed=None, map=imap):
        # TODO -excl, -incl
        self.suffix = mktemp(dir='')

        if isinstance(map, (int, long)):
            import multiprocessing
            p = multiprocessing.Pool(map,  lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
            self.map = p.imap
        else:
            self.map = map

        self.query = mktemp(suffix=".sorted.%s" % self.suffix)
        self.value_fn = value_fn
        self.subject = mktemp(suffix=".sorted.%s" % self.suffix)
        _run("sort -k1,1 -k2,2n %s > %s" % (query, self.query))
        _run("sort -k1,1 -k2,2n %s > %s" % (subject, self.subject))
        self.n = n
        if not os.path.exists(genome):
            self.genome = mktemp(suffix="." + genome + ".%s" % self.suffix)
            _run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, size from %s.chromInfo" > %s' % (genome, self.genome))
        else:
            self.genome = genome

        self.shuffle_str = shuffle_str
        self.seed = 0 if seed is None else seed
    def __del__(self):

        for f in glob.glob("/tmp/*.sorted.%s" % self.suffix):
            os.unlink(f)
        if self.genome.endswith(self.suffix):
            os.unlink(self.genome)

    def run(self, command="bedtools jaccard -a %(query)s -b %(subject)s",
            sims=False):
        self.obs = _command(command, dict(query=self.query,
            subject=self.subject), self.value_fn)

        # call external functions self.map may be Pool.imap
        # must use list-comp here not generator to make sure they are same
        # order for same seed.
        sims = self.map(_shuffle_and_run_star, [(i + self.seed, self.shuffle_str, self.query,
            self.subject, self.genome, command, self.value_fn) for i in
            xrange(self.n)])

        #for i in xrange(self.n):
        #    sims.append(_shuffle_and_run(self.shuffle_str,
        #        self.query, self.subject, self.genome, command, self.value_fn))
        if not isinstance(sims, list):
            sims = list(sims)
        compare = self.compare(sims)
        if sims: compare['sims'] = sims
        return compare


    @classmethod
    def sim_compare(cls, obs, sims_output):
        n = len(sims_output)
        d = dict(observed=obs,
                n_sims=n,
                n_sims_gt=sum(1 for s in sims_output if s > obs),
                n_sims_eq=sum(1 for s in sims_output if s == obs),
                n_sims_lt=sum(1 for s in sims_output if s < obs),
            )
        d['p_sims_gt'] = d['n_sims_gt'] / float(n)
        d['p_sims_lt'] = d['n_sims_lt'] / float(n)
        return d

    def compare(self, sims_output):
        d = Shuffler.sim_compare(self.obs, sims_output)
        d['value_fuction'] = self.value_fn.func_name
        return d



def _shuffle_and_run_star(args):
    #print "seed:", args[0],
    #sys.stdout.flush()
    random.seed(args[0])
    return _shuffle_and_run(*args[1:])

def _shuffle_and_run(shuffle_str, query, subject, genome, command, value_fn):
    temp = _shuffle(shuffle_str, query, genome)
    args_dict = dict(query=temp, subject=subject)
    res = _command(command, args_dict, value_fn)
    os.unlink(temp)
    return res

def _command(command, args_dict, value_fn):
    res = nopen("|%s" % (command % args_dict))
    return value_fn(res)

def _shuffle(shuffle_str, query, genome):
    temp = mktemp(suffix=".shuffled")
    bed_seed = random.randint(0, 99999999999)
    _run('bedtools shuffle -seed %i %s -i %s -g %s | sort -k1,1 -k2,2n > %s' %
            (bed_seed, shuffle_str,
            query, genome, temp))
    return temp

# these are examples of value functions. they take an iterator
# that is the return of the command specified to Shuffler.run()
# so, e.g. when using | wc -l, num_intersections() simply takes
# the first value and returns it as an int.
def jaccard_length(res):
    for i, row in enumerate(res):
        row = row.split("\t")
        if row[0].isdigit(): return int(row[0])

def jaccard_value(res):
    for i, row in enumerate(res):
        row = row.split("\t")
        if row[0].isdigit(): return float(row[2])

def num_intersections(res):
    r = res.next()
    return int(r)

if __name__ == "__main__":

    SEED = 12221
    N_SIMS = 1500

    BASE = "/home/brentp/with_Brent/"

    import sys
    sys.path.insert(0, "%s/LOH_age/src/" % BASE)
    import lohcna

    imap = 18


    def shuff_compare(fname):
        early = Shuffler(fname,
            '%s/LOH_repli/data/features/data_c_constant_early.bed' % BASE, 'hg18',
            jaccard_length, n=N_SIMS, seed=SEED, map=imap).run(sims=True)

        late = Shuffler(fname,
            '%s/LOH_repli/data/features/data_e_constant_late.bed' % BASE, 'hg18',
            jaccard_length, n=N_SIMS, seed=SEED, map=imap).run(sims=True)

        return Shuffler.sim_compare(early['observed'] / float(late['observed']), [e
                        / float(l) for e, l in zip(early['sims'],
                            late['sims'])])['p_sims_gt']


    fnames = reader('%s/LOH_repli/data/filelist_GBM_f0_HAIB__HumanHap550.txt' % BASE,
            header=False)
    for i, (fname,) in enumerate(fnames):
        name = mktemp()
        lohcna.to_bed("%s/LOH_repli/data/%s" % (BASE, fname), name, lohcna.loh_fn)
        print shuff_compare(name)
        if i > 8: break

    #print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
    #    '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
    #    jaccard_value, n=100, shuffle_str="-excl ~/with_Brent/LOH_repli/data/features/bedtools.hg18.centromere").run()
    1/0


    print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        jaccard_length, n=10).run()

    print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        num_intersections, n=100).run("bedtools intersect -a %(query)s -b %(subject)s -u | wc -l")
