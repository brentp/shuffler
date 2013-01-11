from toolshed import reader, nopen
from tempfile import mktemp
from itertools import imap
import random
import os, glob
import sys
import signal

def _run(cmd):
    list(nopen("|%s" % cmd.lstrip("|")))

class Shuffler(object):
    """
    given a query and a subject
    compare the original [value] of query and subject with the
    values after n shufflings of the query compared to the same subject.
    """
    domain = None
    def __init__(self, query, subject, genome, value_fn, n=10, shuffle_str="",
            seed=None, map=imap, temp_dir="/tmp/"):
        # TODO -excl, -incl
        self.suffix = mktemp(dir='')

        self.temp_dir = temp_dir

        if isinstance(map, (int, long)):
            import multiprocessing
            p = multiprocessing.Pool(map,  lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
            self.map = p.imap
        else:
            self.map = map

        self.query = mktemp(suffix=".sorted.%s" % self.suffix, dir=self.temp_dir)
        self.value_fn = value_fn
        self.subject = mktemp(suffix=".sorted.%s" % self.suffix, dir=self.temp_dir)
        _run("sort -k1,1 -k2,2n %s > %s" % (query, self.query))
        _run("sort -k1,1 -k2,2n %s > %s" % (subject, self.subject))
        self.n = n
        if not os.path.exists(genome):
            self.genome_file = mktemp(suffix="." + genome + ".%s" % self.suffix,
                    dir=self.temp_dir)
            Shuffler.genome(genome, self.genome_file)
        else:
            self.genome_file = genome

        self.shuffle_str = shuffle_str
        self.seed = random.randint(0, sys.maxint) if seed is None else seed

    def set_domain(self, domain_file, pad=5000):
        """
        file(s) that indicate(s) which regions for which to allow shuffling.
        these are also applied to the query and the subject.
        """
        # TODO merge or union domain after multiple calls.
        self.domain = mktemp(suffix=".domain.%s" % self.suffix, dir=self.temp_dir)
        _run("bedtools slop -b %i -g %s -i %s > %s" % (pad, self.genome_file, domain_file, self.domain))
        for attr in ('query', 'subject'):
            new_file = mktemp(suffix=".sorted.%s" % self.suffix, dir=self.temp_dir)
            old_file = getattr(self, attr)
            _run("bedtools intersect -a %s -b %s -u | sort -k1,1 -k2,2n > %s" %
                                (old_file, self.domain, new_file))
            setattr(self, attr, new_file)
            os.unlink(old_file)
        self.shuffle_str = "-incl %s" % self.domain

    @classmethod
    def genome(cls, genome, outf):
        _run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, size from %s.chromInfo" > %s' % (genome, outf))
        return outf

    @classmethod
    def plot(self, sims, png):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("%i shufflings" % len(sims['sims']))
        ax.set_xlabel(sims['value_function'])
        ax.set_ylabel('count')
        n, bins, patches = ax.hist(sims['sims'], 20, normed=0, facecolor='green', alpha=0.75)
        red = ax.axvline(x=sims['observed'], color='r')
        ax.legend( (red, patches[0]), ('observed', 'simulated'))

        plt.savefig(png)

    def __del__(self):

        if getattr(self, "domain"):
            os.unlink(self.domain)

        for f in glob.glob("%s/*.sorted.%s" % (self.temp_dir, self.suffix)):
            os.unlink(f)
        if getattr(self, "genome_file", "").endswith(self.suffix):
            os.unlink(self.genome_file)

    def run(self, command="bedtools jaccard -a %(query)s -b %(subject)s",
            sims=False):
        args = dict(query=self.query, subject=self.subject)
        #print command % args
        self.obs = self.value_fn(nopen("|%s" % (command % args)))

        # call external functions self.map may be Pool.imap
        # must use list-comp here not generator to make sure they are same
        # order for same seed.
        sim_list = self.map(_shuffle_and_run_star, [(i + self.seed, self.shuffle_str, self.query,
            self.subject, self.genome_file, self.temp_dir, command, self.value_fn) for i in
            xrange(self.n)])

        if not isinstance(sim_list, list):
            sim_list = list(sim_list)
        compare = self.compare(sim_list)
        if sims: compare['sims'] = sim_list
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
        d['value_function'] = self.value_fn.func_name
        return d

def _shuffle_and_run_star(args):
    #print "seed:", args[0],
    #sys.stdout.flush()
    random.seed(args[0])
    return _shuffle_and_run(*args[1:])

def _shuffle_and_run(shuffle_str, query, subject, genome, temp_dir, command, value_fn):
    bed_seed = random.randint(0, sys.maxint)
    shuffle_cmd = 'bedtools shuffle -seed %i %s -i %s -g %s | sort -k1,1 -k2,2n ' \
            % (bed_seed, shuffle_str, query, genome)
    full_command = "%s |  %s" % (shuffle_cmd, command)
    args_dict = dict(query="-", subject=subject)
    res_iter = nopen("|%s" % full_command % args_dict)
    value = value_fn(res_iter)
    return value

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

    SEED = 122212
    N_SIMS = 1000

    BASE = "/home/brentp/with_Brent/"

    import sys
    sys.path.insert(0, "%s/LOH_age/src/" % BASE)
    import lohcna

    import multiprocessing
    p = multiprocessing.Pool(12,  lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
    imap = p.imap
    genome = '/dev/shm/hg18.genome'
    if not os.path.exists(genome):
        genome = Shuffler.genome('hg18', genome)

    def shuff_compare(fname, domain=None):
        shuff_str = "-maxTries 10"
        early = Shuffler(fname,
            '%s/LOH_repli/data/features/data_c_constant_early.bed' % BASE,
            genome,
            jaccard_length,
            #shuffle_str=shuff_str,
            n=N_SIMS, seed=SEED, map=imap, temp_dir="/dev/shm/")
        if domain:
            early.set_domain(domain)
        early = early.run(sims=True)

        late = Shuffler(fname,
            '%s/LOH_repli/data/features/data_e_constant_late.bed' % BASE,
            genome,
            jaccard_length,
            #shuffle_str=shuff_str,
            n=N_SIMS, seed=SEED, map=imap, temp_dir="/dev/shm/")
        if domain:
            late.set_domain(domain)

        late = late.run(sims=True)

        try:
            return Shuffler.sim_compare(early['observed'] / float(late['observed']), [e
                        / float(l) for e, l in zip(early['sims'],
                            late['sims'])])['p_sims_gt']
        except:
            return "NA"

    # res = open('OV.txt', 'w')
    # '%s/LOH_repli/data/filelist_OV_f0_HAIB__Human1MDuo.txt'
    fnames = [x[0] for x in reader('%s/LOH_repli/data/filelist_GBM_f0_HAIB__HumanHap550.txt'
            % BASE, header=False)]
    """
    res = open('%s/LOH_repli/data/GBMall.txt' % BASE, 'w')
    for i, f in enumerate(fnames):
        for j, line in enumerate(open('%s/LOH_repli/data/%s' % (BASE, f))):
            if j == 0 and i > 0: continue
            res.write(line)
    res.close()
    fnames = ['GBMall.txt']
    """

    for i, fname in enumerate(fnames):
        name = mktemp(dir="/dev/shm/")
        lohcna.to_bed("%s/LOH_repli/data/%s" % (BASE, fname), name, lohcna.loh_fn)
        domain = 'domain.bed'
        _run("sort -k1,1 -k2,2n %s | bedtools merge -d 5000 -i - > %s" %
                (name, domain))
        pair = (fname, shuff_compare(name, domain))
        print >>sys.stderr, pair
        #print >>res, "%s\t%s" % pair
        os.unlink(name)

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

