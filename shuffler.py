from toolshed import reader, nopen
from tempfile import mktemp
import random

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
            seed=None):
        # TODO -excl, -incl

        self.query = query + ".sorted"
        self.value_fn = value_fn
        self.subject = subject + ".sorted"
        _run("sort -k1,1 -k2,2n %s > %s" % (query, self.query))
        _run("sort -k1,1 -k2,2n %s > %s" % (subject, self.subject))
        self.n = n
        self.genome = mktemp(suffix="." + genome)
        _run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, size from %s.chromInfo" > %s' % (genome, self.genome))
        self.shuffle_str = shuffle_str
        self.seed = seed

    def run(self, command="bedtools jaccard -a %(query)s -b %(subject)s",
            sims=False):
        random.seed(self.seed)
        self.obs = self.command(command, dict(query=self.query,
            subject=self.subject))
        sims = []
        for i in xrange(self.n):
            temp = self.shuffle()
            sims.append(self.command(command, dict(query=temp, subject=self.subject)))
        compare = self.compare(sims)
        if sims: compare['sims'] = sims
        return compare

    def command(self, command, args_dict):
        res = nopen("|%s" % (command % args_dict))
        return self.value_fn(res)

    def shuffle(self):
        temp = mktemp(suffix=".shuffled")
        bed_seed = random.randint(0, 99999999999)
        _run('bedtools shuffle -seed %i %s -i %s -g %s | sort -k1,1 -k2,2n > %s' %
                (bed_seed, self.shuffle_str,
                self.query, self.genome, temp))
        return temp

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

    SEED = 1123
    N_SIMS = 200

    early = Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        jaccard_length, n=N_SIMS, seed=SEED).run(sims=True)
    late = Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_e_constant_late.bed', 'hg18',
        jaccard_length, n=N_SIMS, seed=SEED).run(sims=True)

    print Shuffler.sim_compare(early['observed'] / float(late['observed']), [e
                        / float(l) for e, l in zip(early['sims'], late['sims'])])


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

