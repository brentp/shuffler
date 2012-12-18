from toolshed import reader, nopen
from tempfile import mktemp

def _run(cmd):
    if cmd.startswith("|"): cmd = cmd[1:]
    list(nopen("|%s" % cmd))

class Shuffler(object):
    """
    given a query and a subject
    compare the original [value] of query and subject with the
    values after n shufflings of the query compared to the same subject.
    see below for examples of value_fn.
    """
    def __init__(self, query, subject, genome, value_fn, n=10, shuffle_str=""):
        # TODO -excl, -incl

        self.query = query + ".sorted"
        self.value_fn = value_fn
        self.subject = subject + ".sorted"
        _run("sort -k1,1 -k2,2n %s > %s" % (query, self.query))
        _run("|sort -k1,1 -k2,2n %s > %s" % (subject, self.subject))
        self.n = n
        self.genome = mktemp(suffix="." + genome)
        _run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, size from %s.chromInfo" > %s' % (genome, self.genome))
        self.shuffle_str = shuffle_str

    def run(self, command="bedtools jaccard -a %(query)s -b %(subject)s"):
        self.obs = self.command(command, dict(query=self.query,
            subject=self.subject))
        sims = []
        for i in xrange(self.n):
            temp = self.shuffle()
            sims.append(self.command(command, dict(query=temp, subject=self.subject)))
        return self.compare(sims)

    def command(self, command, args_dict):
        res = nopen("|%s" % (command % args_dict))
        return self.value_fn(res)

    def shuffle(self):
        temp = mktemp(suffix=".shuffled")
        _run('bedtools shuffle %s -i %s -g %s | sort -k1,1 -k2,2n > %s' %
                (self.shuffle_str,
                self.query, self.genome, temp))
        return temp

    def compare(self, sims_output):
        d = dict(observed=self.obs,
                n_sims=self.n,
                n_sims_gt=sum(1 for s in sims_output if s > self.obs),
                n_sims_eq=sum(1 for s in sims_output if s == self.obs),
                n_sims_lt=sum(1 for s in sims_output if s < self.obs),
            )
        d['p_sims_gt'] = d['n_sims_gt'] / float(self.n)
        d['p_sims_lt'] = d['n_sims_lt'] / float(self.n)
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

    print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        jaccard_value, n=100).run()
    print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        jaccard_length, n=100).run()

    print Shuffler('/tmp/hudsonalpha.org__HumanHap550__TCGA-02-0028-01A-01D-0184-06__snp_analysis.loh.txt.bed',
        '~/with_Brent/LOH_repli/data/features/data_c_constant_early.bed', 'hg18',
        num_intersections, n=100).run("bedtools intersect -a %(query)s -b %(subject)s | wc -l")
