# see: http://metatracks.encodenets.gersteinlab.org/
# cite: PMID 22950945
mkdir -p tmp/
cd tmp/
if [ ! -e ALL_All_merged.tar.gz ]; then
    wget http://metatracks.encodenets.gersteinlab.org/ALL_All_merged.tar.gz
fi
tar xzvf ALL_All_merged.tar.gz

cells="Gm12878 H1hesc Helas3 Hepg2 K562"
regions="HOT HOT_distal LOT LOT_distal BAR BIR DRM PRM"

OUT=gerstein-binding-sites.bed
rm -f $OUT
for region in $regions; do
    echo $region
    files=""
    for cell in $cells; do
        files="$files ${region}_${cell}_merged.bed "
    done
    # require at least 3 cell types in a given interval
    # extract the position and the cell-types list
    bedtools multiinter -i $files -names $cells \
        | awk -vregion=$region 'BEGIN{OFS=FS="\t"} $4 > 3 { print $1,$2,$3,$5,region }' \
    >> $OUT
        #| bedtools expand -c 4 \
done

sort -k1,1 -k2,2n $OUT > $OUT.tmp
mv $OUT.tmp $OUT
mv $OUT ../

