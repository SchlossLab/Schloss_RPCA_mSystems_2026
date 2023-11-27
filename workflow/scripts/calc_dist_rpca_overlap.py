#!/usr/bin/env python3

# much of this adopted from deicode_benchmarking/simulations/scripts/cluster_results.py
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

depth_ = int(sys.argv[1])

subsampled = pd.read_csv('deicode_benchmarking/simulations/cluster_models/simulation_subsampled_noisy.csv.gz', index_col=[0,1,2,3])

depths_ = set(subsampled.index.get_level_values('sequence_depth'))

rank_ = 2
power_ = 20
# depth_ = 1000

subtmp=subsampled.loc[(rank_,power_,depth_,),:].copy().T 
subtmp_sub = subtmp.copy()

X_sparse=rclr().fit_transform(subtmp_sub.copy())
U,s,V=OptSpace(iteration=5, rank = 2).fit_transform(X_sparse)

distance.cdist(U,U)
df = pd.DataFrame(distance.cdist(U,U))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/overlap_rpca.{depth_}.euclidean.dist.gz', sep = '\t', index = False)

explained = sorted([s[0,0], s[1,1]], reverse=True)
a1 = round(100 * explained[0]/sum(explained), ndigits = 2)
a2 = round(100 * explained[1]/sum(explained), ndigits = 2)

pcoa = pd.DataFrame(U)
pcoa.columns = [f"axis1_{a1}", f"axis2_{a2}"]
pcoa.insert(0, 'sample', subtmp.index)
pcoa.to_csv(f'data/simulation/overlap_rpca.{depth_}.euclidean.pcoa.axes', sep = '\t', index = False)


X_spn = np.array(subtmp_sub.copy()).astype(float)
X_spn[X_spn==0] = np.nan
U_,s_,V_=OptSpace(iteration=5, rank = 2).fit_transform(X_spn)

df = pd.DataFrame(distance.cdist(U_,U_))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/overlap_optspace.{depth_}.euclidean.dist.gz', sep = '\t', index = False)

explained = sorted([s_[0,0], s_[1,1]], reverse=True)
a1 = round(100 * explained[0]/sum(explained), ndigits = 2)
a2 = round(100 * explained[1]/sum(explained), ndigits = 2)

pcoa = pd.DataFrame(U)
pcoa.columns = [f"axis1_{a1}", f"axis2_{a2}"]
pcoa.insert(0, 'sample', subtmp.index)
pcoa.to_csv(f'data/simulation/overlap_optspace.{depth_}.euclidean.pcoa.axes', sep = '\t', index = False)
