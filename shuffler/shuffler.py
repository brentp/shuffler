from toolshed import reader, nopen
import atexit
import tempfile
from itertools import imap
import random
import os
os.environ['LC_ALL'] = 'C'
import sys
import signal

def _run(cmd):
    list(nopen("|%s" % cmd.lstrip("|")))
JACCARD_METRICS = "intersection jaccard n_intersections".split()

def rm(fname):
    try:
        os.unlink(fname)
    except OSError:
        pass

def mktemp(*args, **kwargs):
    f = tempfile.NamedTemporaryFile(*args, **kwargs)
    atexit.register(rm, f.name)
    try:
        return f.name
    finally:
        f.close()

def jaccard_values(res, keys=JACCARD_METRICS):
    for i, row in enumerate(res):
        row = row.split("\t")
        if row[0].isdigit():
            if len(row) > 3:
                return dict(zip(keys,
                    (int(row[0]), float(row[2]), int(row[3]))))
            else:
                return dict(zip(keys,
                    (int(row[0]), float(row[2]))))
    raise Exception("not found")

class Shuffler(object):
    """
    given a query and a subject
    compare the original [value] of query and subject with the
    values after n shufflings of the query compared to the same subject.
    """
    domain = None
    jaccard_metrics = JACCARD_METRICS
    def __init__(self, query, subject, genome, value_fn=jaccard_values, n=10,
            shuffle_str="", seed=None, map=imap,
            temp_dir=os.environ.get("TMPDIR", "/tmp/")):
        self.suffix = "shuffler"
        self.temp_dir = temp_dir

        if isinstance(map, (int, long)):
            import multiprocessing
            p = multiprocessing.Pool(map,
                    lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
            self.map = p.imap
        else:
            self.map = map
        self.value_fn = value_fn
        self.shuffle_str = shuffle_str
        self.seed = random.randint(0, sys.maxint) if seed is None else seed
        self.n = n

        self._prepare(query, subject, genome)

    def _prepare(self, query, subject, genome):
        self.query = mktemp(suffix=".sorted.%s" % self.suffix,
                dir=self.temp_dir)
        self.subject = mktemp(suffix=".sorted.%s" % self.suffix,
                dir=self.temp_dir)

        cut = "| cut -f 1-3 " if self.value_fn == jaccard_values else ""
        _run("awk 'NR > 1 || $1 != \"chrom\"' %s | sort -k1,1 -k2,2n %s > %s" \
                % (query, cut, self.query))
        _run("awk 'NR > 1 || $1 != \"chrom\"' %s | sort -k1,1 -k2,2n %s > %s" \
                % (subject, cut, self.subject))

        if not os.path.exists(genome):
            self.genome_file = mktemp(suffix="." + genome + ".%s" % self.suffix,
                    dir=self.temp_dir)
            Shuffler.genome(genome, self.genome_file)
        else:
            self.genome_file = genome

    def set_domain(self, domain_file, pad=5000):
        """
        file(s) that indicate(s) which regions for which to allow shuffling.
        these are also applied to the query and the subject.
        """
        # TODO merge or union domain after multiple calls.
        self.domain = mktemp(suffix=".domain.%s" % self.suffix,
                dir=self.temp_dir)
        _run("bedtools slop -b %i -g %s -i %s > %s" % (pad, self.genome_file,
                                                    domain_file, self.domain))
        for attr in ('query', 'subject'):
            new_file = mktemp(suffix=".sorted.%s" % self.suffix,
                              dir=self.temp_dir)
            old_file = getattr(self, attr)
            _run("bedtools intersect -a %s -b %s -u | sort -k1,1 -k2,2n > %s" %
                                (old_file, self.domain, new_file))
            setattr(self, attr, new_file)
        self.shuffle_str = "-incl %s" % self.domain

    @classmethod
    def genome(cls, genome, outf):
        _run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, size from %s.chromInfo" > %s' % (genome, outf))
        return outf

    @classmethod
    def plot(self, sims, png=None, ax=None):
        if not 'matplotlib.backends' in sys.modules:
            import matplotlib
            matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            ax.set_title("%i shufflings" % len(sims['sims']))
            #ax.set_xlabel(sims['value_function'])
            ax.set_ylabel('count')
        n, bins, patches = ax.hist(sims['sims'], 20, normed=0,
                facecolor='green', alpha=0.75, edgecolor='green',
                rwidth=1.0)
        red = ax.axvline(x=sims['observed'], color='blue', linewidth=2,
                alpha=0.70)
        if png is not None:
            ax.legend( (red, patches[0]), ('observed', 'simulated'))
            plt.savefig(png)
        return ax

    def run(self, cmd="bedtools jaccard -a %(query)s -b %(subject)s",
            sims=False):
        # jaccard not useful for most things.
        if self.value_fn != jaccard_values and "jaccard" in cmd:
            cmd = "bedtools intersect -sorted -a %(query)s -b %(subject)s -wo"
        args = dict(query=self.query, subject=self.subject)
        self.obs = self.value_fn(nopen("|%s" % (cmd % args)))

        # call external functions self.map may be Pool.imap
        # must use list-comp here not generator to make sure they are same
        # order for same seed.
        sim_list = self.map(_shuffle_and_run_star, [(i + self.seed, self.shuffle_str, self.query,
            self.subject, self.genome_file, self.temp_dir, cmd, self.value_fn) for i in
            xrange(self.n)])

        if not isinstance(sim_list, list):
            sim_list = list(sim_list)
        compare = self.compare(sim_list, sims=sims)
        if "p_sims_gt" in compare:
            compare['value_function'] = self.value_fn.func_name
        return compare

    def compare(self, sims_output, sims=False):
        d = Shuffler.sim_compare(self.obs, sims_output, sims)
        return d

    @classmethod
    def sim_compare(cls, obs, sims_output, sims=False):
        assert isinstance(obs, dict)
        assert len(sims_output[0].keys()) == len(obs.keys())
        n = len(sims_output)
        d = {'n_sims': n}
        for metric in obs.keys():
            d[metric] = dict(
                observed=obs[metric],
                # TODO: sims_output is list of dicts. less memory as dict of lists.
                n_sims_gt=sum(1 for s in sims_output if s[metric] > obs[metric]),
                n_sims_eq=sum(1 for s in sims_output if s[metric] == obs[metric]),
                n_sims_lt=sum(1 for s in sims_output if s[metric] < obs[metric]),
            )
            d[metric]['p_sims_gt'] = d[metric]['n_sims_gt'] / float(n)
            d[metric]['p_sims_lt'] = d[metric]['n_sims_lt'] / float(n)
            if sims:
                d[metric]['sims']=[s[metric] for s in sims_output]
        if len(obs) == 1:
            # if it's a single key we don't use the  metric sub-structure.
            d.update(d.pop(obs.keys()[0]))
        return d


def _shuffle_and_run_star(args):
    #sys.stdout.flush()
    random.seed(args[0])
    return _shuffle_and_run(*args[1:])

def _shuffle_and_run(shuffle_str, query, subject, genome, temp_dir, command,
                     value_fn):
    bed_seed = random.randint(0, sys.maxint)
    shuffle_cmd = 'bedtools shuffle -seed %i %s -i %s -g %s' \
            % (bed_seed, shuffle_str, query, genome)
    full_command = "%s |  %s" % (shuffle_cmd, command)
    args_dict = dict(query="<(sort -k1,1 -k2,2n -)", subject=subject)
    res_iter = nopen("|%s" % full_command % args_dict)
    value = value_fn(res_iter)
    assert isinstance(value, dict)
    return value

if __name__ == "__main__":
    pass
