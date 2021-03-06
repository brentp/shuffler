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

def count_length(bed):
    l = 0
    for toks in bed:
        l += int(toks[2]) - int(toks[1])
    return l


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

def merge_beds(excl_list, genome, prefix="ex"):
    if not os.path.exists(genome):
        fgen = mktemp()
        genome = Shuffler.genome(genome, fgen)

    if len(excl_list) == 1:
        excl = excl_list[0]
    else:
        excl = mktemp()
        _run("|cut -f 1-3 %s | sort -k1,1 -k2,2n | bedtools merge -i - > %s" \
                % (" ".join(excl_list), excl))

    bases = []
    for i, f in enumerate((genome, excl)):
        n_bases = 0
        for toks in reader(f, header=False):
            try:
                if i == 0:
                    n_bases += int(toks[1])
                else:
                    n_bases += (int(toks[2]) - int(toks[1]))
            except ValueError:
                pass
        bases.append(n_bases)

    #print >>sys.stderr, "# %scluding %5g out of %5g total bases (%.3g%%) in the genome" % \
    #        (prefix, bases[1] , bases[0], 100. * bases[1] / float(bases[0]))
    return excl


class Shuffler(object):
    """
    given a query and a subject
    compare the original [value] of query and subject with the
    values after n shufflings of the query compared to the same subject.
    """
    domain = None
    _exclude = _include = None
    jaccard_metrics = JACCARD_METRICS
    def __init__(self, query, subject, genome,
            value_fn=jaccard_values,
            excludes=None, includes=None,
            n=10,
            seed=None, map=imap,
            chrom=False,
            temp_dir=os.environ.get("TMPDIR", "/tmp/"),
            structure=None):
        self.suffix = "shuffler"
        self.temp_dir = temp_dir

        if isinstance(map, (int, long)):
            import multiprocessing
            p = multiprocessing.Pool(map,
                    lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
            self.map = p.imap
        else:
            self.map = map

        self.chrom = chrom
        self.value_fn = value_fn
        self.seed = random.randint(0, sys.maxint) if seed is None else seed
        self.n = n

        self._prepare(query, subject, excludes, includes, genome)
        self._set_structure(structure)

    def _set_structure(self, structure):
        """
        here, we want to intersect the query and subject bed files with the
        structure.bed file and give each set of intervals in query and bed
        that fall within (or have any overlap with) a unique, fake chromosome
        so that all shuffling is within that chromosome.
        in order to do this, we also have to create a fake genome file that
        contains the lengths of those chromosomes.
        """
        if structure in (None, ""): return
        self.chrom = True # has to be by chromosome.

        n_query_before = sum(1 for _ in nopen(self.query))
        n_subject_before = sum(1 for _ in nopen(self.subject))

        new_genome = open(mktemp(suffix='.fake_genome'), 'w')
        structure = "<(cut -f 1-3 %s)" % structure
        seen_segs = {}
        for bed in ('query', 'subject', 'exclude', 'include'):
            bed_path = getattr(self, "_" + bed, getattr(self, bed))
            if not bed_path: continue
            new_fh = open(mktemp(suffix='%s.fake' % bed), 'w')
            for toks in reader("|bedtools intersect -wo -a %s -b '%s' \
                    | sort -k4,4 -k5,5g" % (structure, bed_path), header=False):
                gtoks, btoks = toks[:3], toks[3:-1] # drop the bp overlap
                new_chrom = "_".join(gtoks)

                gtoks[1:] = map(int, gtoks[1:])
                btoks[1:3] = map(int, btoks[1:3])

                glen = gtoks[2] - gtoks[1] # fake chrom length.
                if new_chrom.startswith('chr'): new_chrom = new_chrom[3:]
                if not new_chrom in seen_segs:
                    # save it in the genome file.
                    print >>new_genome, "\t".join((new_chrom, str(glen)))
                seen_segs[new_chrom] = True

                # with partial overlap, we'll have a negative start or an
                # end outside the genome... for now, just truncate.

                # adjust the interval to its location the new chrom.
                btoks[0] = new_chrom
                btoks[1] = max(0, btoks[1] - gtoks[1]) # don't let it go below 0
                # chop to end of fake chrom.
                btoks[2] = min(btoks[2] - gtoks[1], glen - 1)
                assert 0 <= btoks[1] <= btoks[2] < glen
                btoks[1:3] = map(str, btoks[1:3])
                print >>new_fh, "\t".join(btoks)
            new_fh.close()
            setattr(self, bed, new_fh.name)
        new_genome.close()
        self.genome_file = new_genome.name
        #n_query_after = sum(1 for _ in nopen(self.query))
        #n_subject_after = sum(1 for _ in nopen(self.subject))
        #print >>sys.stderr, """applied structure:
        #    interval before and after:
        #    query %i %i
        #    subject %i %i"""\
        #            % (n_query_before, n_query_after,
        #                    n_subject_before, n_subject_after)

    def _prepare(self, query, subject, excludes, includes, genome):
        if not os.path.exists(genome):
            self.genome_file = mktemp(suffix="." + genome + ".%s" % self.suffix,
                    dir=self.temp_dir)
            Shuffler.genome(genome, self.genome_file)
        else:
            self.genome_file = genome


        self.query = mktemp(suffix=".sorted.%s" % self.suffix,
                dir=self.temp_dir)
        self.subject = mktemp(suffix=".sorted.%s" % self.suffix,
                dir=self.temp_dir)

        cut = "| cut -f 1-3 " if self.value_fn == jaccard_values else ""
        _run("sort -k1,1 -k2,2n %s %s > %s" % (query, cut, self.query))
        _run("sort -k1,1 -k2,2n %s %s > %s" % (subject, cut, self.subject))

        self.exclude = excludes
        self.include = includes

    def _set_exclude(self, excludes):
        # TODO: need to remove the q, s feats that overlap with exclue.
        if excludes:
            if isinstance(excludes, basestring): excludes = [excludes]
            #print >>sys.stderr, "merging excludes: ", excludes
            self._exclude = merge_beds(excludes, self.genome_file)

            for qs in ('query', 'subject'):
                bed = getattr(self, qs)
                n_orig = count_length(reader("|bedtools merge -i <(sort -k1,1 -k2,2n %s)" \
                         % bed, header=False))
                ex_bed = mktemp(suffix=qs)
                _run("|bedtools subtract -a %s -b %s | sort -k1,1 -k2,2n > %s" \
                                % (bed, self._exclude, ex_bed))

                n_after = count_length(reader(ex_bed, header=False))
                setattr(self, qs, ex_bed)

                if n_orig - n_after > 0:
                    print >>sys.stderr, "#removing %i of %i (%.3g%%) from %s" % \
                     (n_orig - n_after, n_orig,
                             100. * (n_orig - n_after) / n_orig, qs)
        else:
            self._exclude = ""

    def _get_exclude(self):
        return ("-excl %s" % self._exclude) if self._exclude else ""

    exclude = property(_get_exclude, _set_exclude)

    def _set_include(self, includes):
        if includes:
            if isinstance(includes, basestring): includes = [includes]
            self._include = merge_beds(includes, self.genome_file, prefix="in")
            self.set_domain(self._include, pad=0)
        else:
            self._include = ""

    def _get_include(self):
        return ("-incl %s" % self._include) if self._include else ""

    include = property(_get_include, _set_include)

    @property
    def shuffle_str(self):
        s = " ".join((self.exclude,
                      self.include,
                      "-chrom" if self.chrom else ""))
        return " ".join(s.split()) # remove empties

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

    @classmethod
    def genome(cls, genome, outf=None):
        if os.path.exists(genome): return genome
        if outf is None: outf = mktemp()
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
    shuffle_cmd = 'bedtools shuffle -maxTries 9999 -seed %i %s -i %s -g %s' \
            % (bed_seed, shuffle_str, query, genome)
    full_command = "%s |  %s" % (shuffle_cmd, command)
    args_dict = dict(query="<(sort -k1,1 -k2,2n -)", subject=subject)
    #print >>sys.stderr, full_command % args_dict
    res_iter = nopen("|%s" % full_command % args_dict)
    value = value_fn(res_iter)
    assert isinstance(value, dict)
    return value

if __name__ == "__main__":
    pass
