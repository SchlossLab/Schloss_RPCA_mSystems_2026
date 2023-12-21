#!/usr/bin/env python

# much of this adopted from deicode_benchmarking/simulations/scripts/cluster_results.py
# and deicode_benchmarking/simulations/negative_control.ipynb
#
# In deicode_benchmarking/simulations/scripts/cluster_results.py they used
# iteration=1000, which would seem to be prone to overfitting consdiering
# default is 5 and max recommended is 100. Using 1000 here since 5 give gave
# results that were larger than those reported in the paper (e.g. the Fst values
# were larger by OptSpace than by RPCA)

import pandas as pd
import numpy as np
from deicode.optspace import OptSpace
from scipy.spatial import distance
import sys
import re

count_file = sys.argv[1]
counts = pd.read_csv(count_file, index_col=0, sep = '\t')

base=re.sub("raw.tsv.gz", "", count_file)

subtmp=counts.copy()

X_spn = np.array(subtmp.copy()).astype(float)
X_spn[X_spn==0] = np.nan
opt_rnd=OptSpace(iteration=1000, rank = 2).fit(X_spn)

explained = opt_rnd.explained_variance_ratio
dist = pd.DataFrame(distance.cdist(opt_rnd.sample_weights, opt_rnd.sample_weights))
pcoa = pd.DataFrame(opt_rnd.sample_weights)


dist.columns = subtmp.index
dist.insert(0, 'sample', subtmp.index)
dist.to_csv(f'{base}optspace.euclidean.dist.gz', sep = '\t', index = False)

# deicode_benchmarking/simulations/negative_control.ipynb flips the axes
# relative to the percent explained, see the plot_pcoa function
pcoa.columns = [f"axis2_{round(100*explained[0], ndigits=2)}",
                f"axis1_{round(100*explained[1], ndigits=2)}"]
pcoa.insert(0, "sample", subtmp.index)
pcoa.to_csv(f'{base}optspace.euclidean.pcoa.axes', sep = '\t', index = False)
