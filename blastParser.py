from __future__ import print_function
from Bio import SearchIO
import sys
import re

fname = sys.argv[1]
ofname = fname + ".parseOut"
outhandle = open(ofname, "w")

#### filter param ####
cvrC = 0.80
idtC = 0.85
topN = 3

#### parse blastn output ####
rcds = SearchIO.parse(fname, 'blast-text')

outhandle.write(";".join(['id'] + ["hit"] * topN) + "\n")

for read in rcds:
    if len(read.hits) > topN:
        topC = [h.hit_description for h in read.hsps if h.ident_num /
                float(h.hit_span) > idtC and h.ident_num / float(h.query_span) > cvrC][0:topN]
    elif len(read.hits) <= topN and len(read.hits) > 0:
        topC = [h.hit_description for h in read.hsps if h.ident_num /
                float(h.hit_span) > idtC and h.ident_num / float(h.query_span) > cvrC]
        topC = topC + ["NA"] * (topN - len(topC))
    else:
        topC = ["NA"] * topN
    outhandle.write(";".join([read.id] + topC) + "\n")

outhandle.close()
