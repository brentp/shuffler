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
URL="http://hgdownload.cse.ucsc.edu/goldenPath/%s/encodeDCC/wgEncodeBroadHmm/" % hg

def run(cmd):
    list(nopen("|" + cmd.lstrip("|")))

run("mkdir -p tmp")


celltypes="Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf".split()
os.chdir("tmp/")

for c in celltypes:
    fn = "wgEncodeBroadHmm%sHMM.%s.bed.gz" % (c, hg)
    if not op.exists(fn):
        run("wget -O %s %s/wgEncodeBroadHmm%sHMM.bed.gz" % (fn, URL, c))

segs=set(t[3].replace(" ", "_") for t in reader('wgEncodeBroadHmmHmecHMM.%s.bed.gz' % hg, header=False))
print segs


base_cmd="bedtools multiinter -names %s -i" % " ".join(celltypes)

OUT="chromHMM.%s.bed" % hg
run("rm -f %s" % OUT)

seen = {}
for seg in segs:
    cmd = base_cmd
    # convert 1_Txn_Elongation to Txn_Elongation
    seg_nice = seg.split("_", 1)[1]
    if seg_nice in seen: continue
    seen[seg_nice] = True
    for c in celltypes:
        run("zcat wgEncodeBroadHmm%sHMM.%s.bed.gz | grep -F \"%s\" \
                | cut -f 1-4 | sort -k1,1 -k2,2n > _tmp%s.bed" % (c, hg, seg_nice, c))
        cmd += " _tmp%s.bed" % c

    # keep places where this class had more than 5 cell types supporting an interval
    run("%s | awk -vseg=\"%s\" 'BEGIN{OFS=FS=\"\\t\"} $4 > 5 {print $1,$2,$3,$5,seg }' >> %s" \
                    % (cmd, seg_nice, OUT))


run("sort -k1,1 -k2,2n %s | gzip -c > ../%s.gz; rm %s" % (OUT, OUT, OUT))
os.chdir("..")
print OUT + ".gz"
# now merge if the last files are the same and the intervals are abutting
OUTM=open("chromHMM.%s.merged.bed.gz" % hg, 'w')

for toks in reader("|bedtools groupby -i %s.gz -g 1,5 -c 2,3 -o min,max" % OUT,
           header=False):
    print >>OUTM, "\t".join((toks[0], toks[2], toks[3], toks[1]))
OUTM.close()
print OUTM.name
