#!/usr/bin/env python

# much of this adopted from deicode_benchmarking/simulations/scripts/cluster_results.py
import pandas as pd
import numpy as np
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from scipy.spatial import distance
import sys

depth_ = int(sys.argv[1])

subsampled = pd.read_csv('deicode_benchmarking/simulations/cluster_models/simulation_subsampled_noisy.csv.gz', index_col=[0,1,2,3])

results = {}
depths_ = set(subsampled.index.get_level_values('sequence_depth'))

rank_ = 2
power_ = 20
# depth_ = 1000

subtmp=subsampled.loc[(rank_,power_,depth_,),:].copy().T 
subtmp_sub = subtmp.copy()
#meta on cluster
meta = np.array([1]*int(subtmp.shape[0]/2)+[2]*int(subtmp.shape[0]/2)).T
meta = pd.DataFrame(meta,index=subtmp.index,columns=['group'])

X_sparse=rclr().fit_transform(subtmp_sub.copy())
U,s,V=OptSpace(iteration=1000).fit_transform(X_sparse)
distance.cdist(U,U)
df = pd.DataFrame(distance.cdist(U,U))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/overlap_rpca.{depth_}.euclidean.dist.gz', sep = '\t', index = False)

X_spn = np.array(subtmp_sub.copy()).astype(float)
X_spn[X_spn==0] = np.nan
U_,s_,V_=OptSpace(iteration=1000).fit_transform(X_spn)
df = pd.DataFrame(distance.cdist(U_,U_))
df.columns = subtmp.index
df.insert(0, 'sample', subtmp.index)
df.to_csv(f'data/simulation/overlap_optspace.{depth_}.euclidean.dist.gz', sep = '\t', index = False)
