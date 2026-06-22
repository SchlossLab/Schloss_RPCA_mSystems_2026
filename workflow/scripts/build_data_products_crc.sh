#!/usr/bin/env bash

DATA=$1
PROCESSORSSPLIT=$2
PROCESSORSCLUSTER=1

gunzip -kf data/crc/data.fasta.gz data/crc/data.count_table.gz data/crc/data.taxonomy.gz

mothur "#set.seed(seed=19760620);
  set.dir(input=$DATA);
  cluster.split(fasta=data.fasta, count=data.count_table, taxonomy=data.taxonomy, taxlevel=4, 
                processors=$PROCESSORSSPLIT, runsensspec=FALSE, cluster=FALSE);
  cluster.split(file=current, count=current, cutoff = 0.03,  processors=$PROCESSORSCLUSTER, runsensspec=FALSE);
  make.shared();
  get.oturep(list = data.opti_mcc.list, count = data.count_table, fasta = data.fasta, method = abundance);
  dist.seqs(fasta = data.opti_mcc.0.03.rep.fasta, processors = $PROCESSORSSPLIT, output = lt);
  clearcut(phylip = data.opti_mcc.0.03.rep.phylip.dist)"


mv $DATA/data.opti_mcc.shared $DATA/data.otu.shared
mv $DATA/data.opti_mcc.0.03.rep.phylip.tre $DATA/data.otu.tre
mv $DATA/data.opti_mcc.0.03.rep.count_table $DATA/data.otu.count_table
mv $DATA/data.opti_mcc.0.03.rep.fasta $DATA/data.otu.fasta

gzip -f $DATA/data.otu.shared

rm -f $DATA/data.fasta $DATA/data.count_table $DATA/data.taxonomy
rm -f $DATA/data.opti_mcc.* $DATA/data.file
