command-line interface to compare an observed relationship between 2 sets of
intervals to the distribution of that relationship derived from a number of
shufflings. This "relationship" can be anything that returns a single number,
e.g. total bases of overlap between the 2 sets, or number of overlaps.

parameters: 

    ------
    inputs
    ------

     -a: bed file (the -a file is always shuffled) 

     -b: bed file (the -b file is *NEVER* shuffled) 
     --bs: instead of b, specify a single file and output will be split into 
           groups by the final column (good for finding enrichment in ENCODE
           regions or many user-defined regions with a single call)

    -----------
    constraints 
    -----------
    
     --domain: (optional) shuffle -a intervals inside this domain.
     --include: [same as domain] may be specified multiple times.

     --exclude: (optional) do not allow shuffled intervals to land in the
               intervals specified in this file (may be specified multiple
               times)
    
     --metric: (optional, default=3 metrics of overlap, n_intersections and
               jaccard) may be any bash command, e.g.:
                'awk "BEGIN{s=0}{s+=$3-$2}END{print s}"'])
               or
                'perl my-evaluator.pl'
               or
                'wc -l'
              or anything that accepts intervals to stdin and prints a number
              to stdout: eg:
                'python interval-accepter.py'
              the metric is value to compare the between observed and shuffled.

     -g: genome version (to get chromosomes), e.g. mm8, dm3, hg19 or a file


    ----
    misc
    ----

     -n: (optional, default=1000) number of shufflings
             
     -t: (optional, default=1) number of threads to use.

     --seed: (optional) seed for the random number generator

     --png: (optional) a path to save a png of the simulated distribution
            vs the observed metric.

call like:

    $ python -m shuffler -n 1000 -a query.bed -b subject.bed \
        --domain domain.bed -t 8 -g hg19 --seed 42

output is a single integer indicating the number of simulations
with a metric < the observed

    $ python -m shuffler -a query.bed:500 -b subject.bed:-5000 --domain domain.bed    

this will pad each side of query.bed by 500 bases (this only affects the shuffled data
not the observed).
it will also extend intervals in subject.bed only upstream. if 4 columns are present,
it look for strand there if 6 columns are present it will check for strand there. if
strand is included, up-stream is based on that, otherwise, upstream is "less than".

    $ python -m shuffler -a query.bed -b subject.bed:-5000:100 --domain domain.bed    

will extend the intervals in subject.bed by 5kb upstream and 100 downstream.
for only downstream, use subject.bed:0:100

domains
-------

in the examples above, the --domain is a bed file with a single set of regions
the define the intervals where shuffled query intervals can be placed.

we extend this with --domains.

    $ shuffler -a query.bed -b subject.bed:-5000:100 --domains domains.bed    

here, the program will pull unique values from the 4th column of domains.bed
and create a domain.$value.bed file for each uniqe value. so, e.g. domains.bed
may contain 1000 TSS sites (labelled as "TSS" in the 4th column), 2100 exons,
2001 introns, and 407 UTRs. The program will create a file for each of those
feature types and report the number of shuffled values that are larger than
the observed for each of those types.

domains can also be extended:

    $ shuffler -a query.bed -b subject.bed:-5000:100 --domains domains.bed:10000

will extend each domain intervla by 10KB.


example
-------

We have the chromHMM segmentations in a BED file called: **h1hesc.ChromHMM.hg18.bed** the first few lines look like this:

    chr1    0   200 DnaseU
    chr1    200 400 CtcfO
    chr1    400 5863    Low
    chr1    5863    6263    CtcfO
    chr1    6263    6463    H4K20
    chr1    6463    81263   Quies
    chr1    81263   81463   Ctcf
    chr1    81463   123063  Quies
    chr1    123063  123463  DnaseU
    chr1    123463  128863  Low

where the last column is the name given to that particular chromHMM segmentation.
We also have another file with about 2000 top-secret intervals: **intervals.bed**
and a domain file that lists the possible locations where those intervals could be:
**charm.domain.bed**. We could use --include to specify that, but it's faster to
use --exclude. So we generate the complement with bedtools:

    bedtools complement -i charm.domain.bed -g hg18.genome > charm.domain.complement.bed

then our call to shuffler looks like this (annotated):

    python -m shuffler \
        -a intervals.bed \               # 2000 intervals of interest
        --bs h1hesc.ChromHMM.hg18.bed \  # last column indicates segmentation
        -g hg18.genome \                 # genome file
        -n 250 \                         # number of shufflings
        --seed 42 \                      # seed for reproducibility
        -t 14  \                         # use 14 threads
        --png shuffler.hmm.png \         # save image
        --exclude data/charm.domain.complement.bed  # don't shuffle inside these

and the output has one line for each unique value in the final column of
h1hesc.ChromHMM.hg18.bed:

    group	intersection	jaccard	n_intersections
    Art	0.3	0.3	0.156
    Ctcf	0.84	0.84	0.84
    CtcfO	0.716	0.716	0.36
    DnaseD	0.112	0.112	0.044
    DnaseU	0	0	0
    Elon	0.984	0.984	1
    ElonW	1	1	1
    Enh	0	0	0
    EnhF	0	0	0
    EnhW	0	0	0
    EnhWF	0	0	0
    FaireW	0.916	0.916	0.98
    Gen3'	1	1	1
    Gen5'	0	0	0
    H4K20	0	0	0
    Low	1	1	1
    Pol2	1	1	1
    PromF	0	0	0
    PromP	0.292	0.288	0
    Quies	0.884	0.884	0.688
    Repr	0.012	0.012	0
    ReprD	0	0	0
    ReprW	1	1	0.916
    Tss	1	1	1
    TssF	0	0	0.016

where each column is a "metric". For these, we are most interested
in the  first numeric column--which shows the proportion of times a shuffled
dataset had a higher overlap than the observed. We can see things like
Enhancer have a low probability of occurring by chance while TSS has
a p-value of 1. We can confer with a biologist to decide what this means.
In other cases, we may wish to look at the jaccard index which is the
ratio of intersection to union or the number of (unique) intersections.

Since we specified --png, we get this image:

![distribution](https://gist.github.com/brentp/bf7d3c3d3f23cc319ed8/raw/4d44893a42428172237e9cf84fc9528532be3749/ipf-hmm.png "Distribution of Simulated For intervals")

Where the blue line indicates the observed intersection, number of intersections or jaccard, respectively and
the green shows the distribution of simulations. Remember we got all this from a single command:

    python -m shuffler -a intervals.bed --bs h1hesc.ChromHMM.hg18.bed \
        -g hg18.genome -n 250 --seed 42 -t 14  --png shuffler.hmm.png \
        --exclude data/charm.domain.complement.bed

which ran 250 shufflings for each category from the 2.8 million chromHMM 
segementations on the 2.8 million intervals from chromHMM in about 5 minutes.

todo
----

see genome-wide plot (fig 2) and bi-plot (fig6) in zhang et al:
"Statistical analysis of the genomic distribution and correlation of regulatory elements in the ENCODE regions"
http://genome.cshlp.org/content/17/6/787.full.pdf


clustered-ness of regions within a single file

segway/chrom hmm segmentatnois:
http://www.nature.com/nmeth/journal/v9/n5/full/nmeth.1937.html
