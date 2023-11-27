#!/usr/bin/env python3

from skbio.stats.composition import clr
import pandas as pd
import re
import sys


# this was lifted from deicode_benchmarking/simulations/negative_control.ipynb

count_fname = sys.argv[1]
counts = pd.read_csv(count_fname, index_col=0, sep = '\t')

output_fname = re.sub("raw", "oclr", count_fname)


clr = pd.DataFrame(clr(counts+1),
             columns=['OTU_'+str(x)
                      for x in range(counts.shape[1])],
            index=['sample_'+str(x) 
                      for x in range(counts.shape[0])])

clr.insert(0, 'sample', clr.index)
clr.to_csv(output_fname, sep = '\t', index = False)
