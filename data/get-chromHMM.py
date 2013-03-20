"""
chromHMM was run on each cell type separately
here, we get a "consensus region whenver more than 5 cell types (of 9)
had the same segmentation
"""
from toolshed import nopen, reader
import os
import os.path as op
URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/"

def run(cmd):
    list(nopen("|" + cmd.lstrip("|")))

run("mkdir -p tmp")


celltypes="Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf".split()
os.chdir("tmp/")

for c in celltypes:
    if not op.exists("wgEncodeBroadHmm%sHMM.bed.gz" % c):
        run("wget %s/wgEncodeBroadHmm%sHMM.bed.gz" % (URL, c))

segs=set(t[3] for t in reader('wgEncodeBroadHmmHmecHMM.bed.gz', header=False))
print segs


base_cmd="bedtools multiinter -header -names %s -i" % " ".join(celltypes)

OUT="chromHMM.merged.bed"
run("rm -f %s" % OUT)

seen = {}
for seg in segs:
    cmd = base_cmd
    # convert 1_Txn_Elongation to Txn_Elongation
    seg_nice = seg.split("_", 1)[1]
    if seg_nice in seen: continue
    seen[seg_nice] = True
    for c in celltypes:
        run("zcat wgEncodeBroadHmm%sHMM.bed.gz | grep -F \"%s\" | cut -f 1-4 > _tmp%s.bed" % (c, seg_nice, c))
        cmd += " _tmp%s.bed" % c

    # keep places where this class had more than 5 cell types supporting an interval
    run("%s | awk -vseg=\"%s\" 'BEGIN{OFS=FS=\"\\t\"} $4 > 5 {print $1,$2,$3,$5,seg }' >> %s" \
                    % (cmd, seg_nice, OUT))


run("sort -k1,1 -k2,2n %s | gzip -c > ../%s.gz; rm %s" % (OUT, OUT, OUT))
print OUT + ".gz"
