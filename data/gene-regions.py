import sys
from cruzdb import Genome

db = sys.argv[1]

db = Genome('sqlite:////usr/local/src/cruzdb/%s.db' % db)
refGene = db.refGene

for g in refGene:
    for feat in g.gene_features:
        print "\t".join(map(str, feat))
