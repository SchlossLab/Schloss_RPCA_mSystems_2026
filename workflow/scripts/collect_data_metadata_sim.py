#!/usr/bin/env python3

# this was lifted from deicode_benchmarking/simulations/negative_control.ipynb
# and from deicode_benchmarking/simulations/scripts/cluster_build.py.
# note that the denomenator in the second and third arguments going to
# build_block_model was 15 in cluster_build.py and it was 60 in the
# negative_control.ipynb. I didn't see meaningful differences and went with 15
# since it was used in the overlap simulations, which got more exposure in the
# original paper

import pandas as pd
import numpy as np
from sklearn.utils import shuffle
from deicode_benchmarking.simulations.scripts.simulations import build_block_model
import sys
import random

desired_depth = int(sys.argv[1])
rep = int(sys.argv[2])

# desired_depth = 10000
# rep = 1

# i've been unable to get reproducible results when i peg the RNG seed.
np.random.seed(rep)
random.seed(rep)
np.random.RandomState(rep)


# the number to the right of the comment was the original version in the ipynb
# and script indicated above. I selected larger values that would reliably give
# a minimum sequencing depth per sample at or above the desired sequencing
# depth. I won't do anythign with 500 since it wasn't used in the paper and give
# a 0 negative control matrix (would need to add 1 to the mean abundance)

seq_depth={500:3.05e2, 
            1000:750, #6.1e2,
            2000:1500, #1.25e3,
            4000:2800, #2.5e3,
            10000:7250} #6.05e3

parameter = seq_depth[desired_depth]
rank = 2

nsamples = 200
nsamples_tag = int(nsamples / 2)
notus = 1000

################################################################################

def generate_counts_data(depth_, parameter_, scaling_factor_, overlap_,
                         rank_, nsamples_, notus_):
  
  X_ = []
  min_depth_ = 0
  while(min_depth_ < depth_):
    #run model with fit variables and new variants
    _, X_=build_block_model(rank_,
                            parameter_/scaling_factor_,
                            parameter_/scaling_factor_,
                            parameter_, parameter_,
                            nsamples_, notus_, overlap=overlap_,
                            mapping_on=False)
    X_=pd.DataFrame(X_.T,
                    columns=['OTU_'+str(x)
                              for x in range(notus_)],
                    index=['sample_'+str(x) 
                            for x in range(nsamples_)])
    min_depth_ = X_.sum(axis = "columns").min()
   
  return X_

################################################################################
# now we want to subsample each sample to desired_depth number of sequences and
# make the data frame wide again with samples in rows and otus in columns
################################################################################

def subsample_counts_data(X_, depth_):
  
  X_.insert(0, 'sample', X_.index)
  X_ = X_.melt(id_vars = ["sample"], var_name = "otu")
  X_ = X_.take(X_.index.repeat(X_.value)).drop(columns="value")
  X_ = X_.groupby('sample').sample(n=depth_)
  X_ = X_.value_counts().reset_index(name='count')
  X_ = X_.pivot_table(index = "sample", columns = "otu",
                      values = "count", fill_value = 0).astype(int)
  
  return(X_)

################################################################################
#
# Positive control
#
################################################################################

overlap = 0

X_pos = generate_counts_data(desired_depth, parameter, 60, overlap,
                             rank, nsamples, notus)
X_pos = subsample_counts_data(X_pos, desired_depth)
X_pos.insert(0, 'sample', X_pos.index)

X_pos.to_csv(f'data/positive/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz',
                sep = '\t', index = False)


################################################################################
#
# Negative control - modified code to make it more reliable at finding a
# solution with at least desired_depth in every sample
#
################################################################################

X_neg = []
min_depth = 0

while min_depth < desired_depth :
  X_neg = np.random.randint(0,
                            2 * desired_depth / notus + 2,
                            # np.mean(X_pos.drop(columns = "sample").values)*2.3 + 3,
                            (notus, nsamples))
  X_neg = pd.DataFrame(X_neg,
                      index=['OTU_'+str(x)
                                for x in range(notus)],
                      columns=['sample_'+str(x) 
                                for x in range(nsamples)])
  min_depth = X_neg.sum(axis = "index").min()

# I don't know that the following shuffles are necessary, but they were in the
# original code...
X_neg.index = shuffle(X_neg).index
X_neg.columns = shuffle(X_neg.T).index
X_neg=X_neg.T
X_neg.sort_index(inplace=True)
X_neg=X_neg.T
X_neg.sort_index(inplace=True)
X_neg=X_neg.T

X_neg = subsample_counts_data(X_neg, desired_depth)
X_neg.insert(0, 'sample', X_neg.index)

X_neg.to_csv(f'data/negative/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz',
                sep = '\t', index = False)

################################################################################
#
# Overlap control
#
################################################################################

overlap = 20

X_overlap = generate_counts_data(desired_depth, parameter, 15, overlap,
                                 rank, nsamples, notus)
X_overlap = subsample_counts_data(X_overlap, desired_depth)
X_overlap.insert(0, 'sample', X_overlap.index)

X_overlap.to_csv(f'data/overlap/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz',
                sep = '\t', index = False)

################################################################################
#
# Design files
#
################################################################################

design = pd.DataFrame(["A"] * int(nsamples/2) + ["B"] * int(nsamples/2),
                      columns=["treatment"],
                      index=['sample_'+str(x) 
                              for x in range(nsamples)])
design.insert(0, 'sample', design.index)

design.to_csv(f'data/positive/data.{desired_depth}.{nsamples_tag}_{rep}.design',
                sep = '\t', index = False)
design.to_csv(f'data/negative/data.{desired_depth}.{nsamples_tag}_{rep}.design',
                sep = '\t', index = False)
design.to_csv(f'data/overlap/data.{desired_depth}.{nsamples_tag}_{rep}.design',
                sep = '\t', index = False)
