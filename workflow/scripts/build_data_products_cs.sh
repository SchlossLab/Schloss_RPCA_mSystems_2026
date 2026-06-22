#!/usr/bin/env bash

# note that there are a handful of unique sequences in both datasets that did 
# not align very well. they did not represent a large number of sequences. i
# just left them in 

FASTA=$1
COUNT_TABLE=$2

# FASTA=data/sleep_apnea/data_treatment.fasta
# COUNT_TABLE=data/sleep_apnea/data_treatment.tsv.gz

gunzip -k $COUNT_TABLE

mothur "#align.seqs(fasta = $FASTA, reference = data/references/silva.v4.align, flip = TRUE, processors = 1);
  filter.seqs(vertical = TRUE);
  dist.seqs(output = phylip);
  clearcut()"


PRETTY_TREE=`echo $FASTA | sed "s/fasta/tre/"`
MOTHUR_TREE=`echo $FASTA | sed "s/fasta/filter.phylip.tre/"`
mv $MOTHUR_TREE $PRETTY_TREE

rm -f `echo $FASTA | sed "s/fasta/filter*/"`
rm -f `echo $FASTA | sed "s/fasta/flip.accnos/"`
rm -f `echo $FASTA | sed "s/fasta/align*/"`
rm -f `echo $COUNT_TABLE | sed "s/.gz//"`
