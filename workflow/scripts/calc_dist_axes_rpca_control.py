#!/usr/bin/env python

# much of this adopted from deicode_benchmarking/simulations/scripts/cluster_results.py
# and deicode_benchmarking/simulations/negative_control.ipynb
#
# pegged iterations=5 and rand=2 to avoid overfitting. In deicode_benchmarking/simulations/scripts/cluster_results.py
# iteration=1000, which would seem to be prone to overfitting consdiering
# default is 5 and max recommended is 100

import pandas as pd
import numpy as np
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from scipy.spatial import distance
import sys
import re

count_file = sys.argv[1]
counts = pd.read_csv(count_file, index_col=0, sep = '\t')

base=re.sub("raw.tsv.gz", "", count_file)

subtmp=counts.copy()
subtmp_sub = subtmp.copy()

X_sparse=rclr().fit_transform(subtmp.copy())
opt_sig=OptSpace(iteration=5, rank = 2).fit(X_sparse)

explained = opt_sig.explained_variance_ratio
dist = pd.DataFrame(distance.cdist(opt_sig.sample_weights, opt_sig.sample_weights))
pcoa = pd.DataFrame(opt_sig.sample_weights)

dist.columns = subtmp.index
dist.insert(0, 'sample', subtmp.index)
dist.to_csv(f'{base}rpca.euclidean.dist.gz', sep = '\t', index = False)

# deicode_benchmarking/simulations/negative_control.ipynb flips the axes
# relative to the percent explained, see the plot_pcoa function
pcoa.columns = [f"axis2_{round(100*explained[0], ndigits=2)}",
                f"axis1_{round(100*explained[1], ndigits=2)}"]
pcoa.insert(0, "sample", subtmp.index)
pcoa.to_csv(f'{base}rpca.euclidean.pcoa.axes', sep = '\t', index = False)



X_spn = np.array(subtmp.copy()).astype(float)
X_spn[X_spn==0] = np.nan
opt_rnd=OptSpace(iteration=5, rank = 2).fit(X_spn)

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
dist.to_csv(f'{base}optspace.euclidean.pcoa.axes', sep = '\t', index = False)
