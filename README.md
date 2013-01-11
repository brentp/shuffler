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

    -----------
    constraints 
    -----------
    
     --domain: (optional) shuffle -a intervals inside this domain.
     --include: [same as domain] may be specified multiple times.

     --exclude: (optional) do not allow shuffled intervals to land in the
               intervals specified in this file (may be specified multiple
               times)
    
     --metric: (optional, default=overlap [n_overlaps, etc. also accepted. can
              also be a bash command, e.g.:
                'awk "BEGIN{s=0}{s+=$3-$2}END{print s}"'])
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


call like:

    $ shuffler -n 1000 -a query.bed -b subject.bed \
        --domain domain.bed -t 8 -g hg19 --seed 42

output is a single integer indicating the number of simulations
with a metric < the observed

    $ shuffler -a query.bed:500 -b subject.bed:-5000 --domain domain.bed    

this will pad each side of query.bed by 500 bases (this only affects the shuffled data
not the observed).
it will also extend intervals in subject.bed only upstream. if 4 columns are present,
it look for strand there if 6 columns are present it will check for strand there. if
strand is included, up-stream is based on that, otherwise, upstream is "less than".

    $ shuffler -a query.bed -b subject.bed:-5000:100 --domain domain.bed    

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



todo
----

see genome-wide plot (fig 2) and bi-plot (fig6) in zhang et al:
"Statistical analysis of the genomic distribution and correlation of regulatory elements in the ENCODE regions"
http://genome.cshlp.org/content/17/6/787.full.pdf


clustered-ness of regions within a single file
