#!/usr/bin/env python3

import numpy as np
import pandas as pd
from biom import Table
from gemelli.rpca import rpca #, auto_rpca
import re
import sys

count_fname = sys.argv[1]
rank = sys.argv[2]

if(rank != "auto"):
  rank = int(rank)

# count_fname="data/sleep_apnea/data.11879.87_6.raw.tsv.gz"
# rank=2

dist_fname = re.sub("raw.tsv", f"gemelli_{rank}.euclidean.dist", count_fname)
axes_fname = re.sub("raw.tsv.gz", f"gemelli_{rank}.euclidean.pcoa.axes", count_fname)

counts = pd.read_csv(count_fname, index_col=0, sep = '\t').T

biom = Table(counts.to_numpy(),
             counts.index,
             counts.columns) 

ordination, distances = rpca(biom, n_components = rank, max_iterations = 5)

distances = distances.to_data_frame()
distances.insert(0, 'sample', distances.index)
distances.to_csv(dist_fname, sep='\t', index = False)


explained = ordination.proportion_explained

percent_explained=[]
pc_names=[]
new_names=[]
index = 1
for e in explained:
  if(e > 0):
    percent_explained = round(100 * e, ndigits = 2)
    pc_names.append(f"PC{index}")
    new_names.append(f"axis{index}_{percent_explained}")
    index = index+1


axes = ordination.samples[pc_names]
axes.columns = new_names
axes.insert(0, 'sample', axes.index)
axes.to_csv(axes_fname, sep='\t', index = False)