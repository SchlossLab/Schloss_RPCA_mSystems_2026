#!/usr/bin/env python3

# this was lifted from deicode_benchmarking/simulations/negative_control.ipynb
import pandas as pd
import numpy as np
from sklearn.utils import shuffle
from deicode_benchmarking.simulations.scripts.simulations import build_block_model
import sys

desired_depth = int(sys.argv[1])

seq_depth={500:3.05e2, 
            1000:6.1e2, 
            2000:1.25e3, 
            4000:2.5e3, 
            10000:6.05e3}

depth = seq_depth[desired_depth]
overlap_ = 0
rank_ = 2

#run model with fit variables and new variants
_,X_signal=build_block_model(rank_,  depth/60, depth/60, 
                             depth, depth
                             ,200,1000,overlap=overlap_
                             ,mapping_on=False)
X_signal=pd.DataFrame(X_signal.T,
                      columns=['OTU_'+str(x)
                                for x in range(X_signal.shape[0])],
                      index=['sample_'+str(x) 
                               for x in range(X_signal.shape[1])])
X_signal.insert(0, 'sample', X_signal.index)
X_signal.to_csv(f'data/simulation/positive.{desired_depth}.raw.tsv.gz', sep = '\t', index = False)


#run model with fit variables and new variants
X_random=np.random.randint(0,np.mean(X_signal.drop(columns = "sample").values)*2.3,(1000,200))
X_random=pd.DataFrame(X_random,
                      index=['OTU_'+str(x)
                                for x in range(X_random.shape[0])],
                      columns=['sample_'+str(x) 
                               for x in range(X_random.shape[1])])
X_random.index = shuffle(X_random).index
X_random.columns = shuffle(X_random.T).index
X_random=X_random.T
X_random.sort_index(inplace=True)
X_random=X_random.T
X_random.sort_index(inplace=True)
X_random=X_random.T
X_random.insert(0, 'sample', X_random.index)
X_random.to_csv(f'data/simulation/negative.{desired_depth}.raw.tsv.gz', sep = '\t', index = False)
