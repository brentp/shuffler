"""
chromHMM was run on each cell type separately
here, we get a "consensus region whenver more than 5 cell types (of 9)
had the same segmentation
"""
from toolshed import nopen, reader
import os
import os.path as op
from operator import itemgetter
from itertools import groupby
import sys

try:
    hg = sys.argv[1]
except:
    print >>sys.stderr, "send in hg18 or hg19 as first argument"
    sys.exit(1)
URL="http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/%s.combined.bb"

def run(cmd):
    print cmd
    list(nopen("|" + cmd.lstrip("|")))

run("mkdir -p tmp")


celltypes="Gm12878 H1hesc helas3 hepg2 huvec k562".lower().split()
os.chdir("tmp/")

fns = []
for c in celltypes:
    fn = "%s.%s.combined.bb" % (c, hg)
    fb = fn.replace(".bb", ".bed.gz")
    if not op.exists(fn):
        run("wget -O %s %s" % (fn, URL % c))
    if not op.exists(fb):
        run("bigBedToBed %s stdout | cut -f 1-4 | gzip -c > %s" % (fn, fb))
    print(fb)
    fns.append(fb)

segs=set(t[3].replace(" ", "_") for t in reader(fns[0], header=False))

print segs


base_cmd="bedtools multiinter -names %s -i" % " ".join(celltypes)

OUT="chromHMM-segway-combined.%s.bed" % hg
run("rm -f %s" % OUT)

seen = {}
for seg in segs:
    cmd = base_cmd
    # convert 1_Txn_Elongation to Txn_Elongation
    seg_nice = seg.split("_", 1)[-1]
    if seg_nice in seen: continue
    seen[seg_nice] = True
    for c in celltypes:
        fname = "%s.%s.combined.bed.gz" % (c, hg)
        run("zcat %s | grep -F \"%s\" \
                | cut -f 1-4 | sort -k1,1 -k2,2n > _tmpcomb%s.bed" % (fname, seg_nice, c))
        cmd += " _tmpcomb%s.bed" % c

    # keep places where this class had more than 5 cell types supporting an interval
    run("%s | awk -vseg=\"%s\" 'BEGIN{OFS=FS=\"\\t\"} $4 > 5 {print $1,$2,$3,$5,seg }' >> %s" \
                    % (cmd, seg_nice, OUT))


run("sort -k1,1 -k2,2n %s | gzip -c > ../%s.gz; rm %s" % (OUT, OUT, OUT))
os.chdir("..")
print OUT + ".gz"
# now merge if the last files are the same and the intervals are abutting
OUTM=open("chromHMM-segway-combined.%s.merged.bed.gz" % hg, 'w')

for toks in reader("|bedtools groupby -i %s.gz -g 1,5 -c 2,3 -o min,max" % OUT,
           header=False):
    print >>OUTM, "\t".join((toks[0], toks[2], toks[3], toks[1]))
OUTM.close()
print OUTM.name
