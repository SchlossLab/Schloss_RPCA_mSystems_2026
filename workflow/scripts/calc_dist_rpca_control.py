#!/usr/bin/env python

# much of this adopted from deicode_benchmarking/simulations/scripts/cluster_results.py
import pandas as pd
import numpy as np
import re
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from scipy.spatial import distance
import sys

count_file = sys.argv[1]
counts = pd.read_csv(count_file, index_col=0, sep = '\t')
control = re.sub("data/simulation/(.*)_raw.*", "\\1", count_file)

subtmp=counts.copy()
subtmp_sub = subtmp.copy()

X_sparse=rclr().fit_transform(subtmp_sub.copy())
U,s,V=OptSpace(iteration=1000).fit_transform(X_sparse)
distance.cdist(U,U)
df = pd.DataFrame(distance.cdist(U,U))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/{control}_rpca.4000.euclidean.dist.gz', sep = '\t', index = False)

X_spn = np.array(subtmp_sub.copy()).astype(float)
X_spn[X_spn==0] = np.nan
U_,s_,V_=OptSpace(iteration=1000).fit_transform(X_spn)
df = pd.DataFrame(distance.cdist(U_,U_))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/{control}_optspace.4000.euclidean.dist.gz', sep = '\t', index = False)
