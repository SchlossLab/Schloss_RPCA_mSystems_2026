#!/usr/bin/env bash


wget -N -P data/references/silva https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_2.tgz
tar xvzf data/references/silva/silva.nr_v138_2.tgz -C data/references/silva

#these coordinates remove the primers
mothur "#pcr.seqs(fasta=data/references/silva/silva.nr_v138_2.align, start=13862, end=23444, keepdots=F, processors=1);
    degap.seqs();
    unique.seqs()"

grep ">" data/references/silva/silva.nr_v138_2.pcr.ng.unique.fasta | cut -c 2-  > data/references/silva/silva.nr_v138_2.pcr.ng.accnos

mothur "#get.seqs(fasta=data/references/silva/silva.nr_v138_2.pcr.align, accnos=data/references/silva/silva.nr_v138_2.pcr.ng.accnos);
    screen.seqs(minlength=240, maxlength=275, maxambig=0, maxhomop=8, processors=1)"

mv data/references/silva/silva.nr_v138_2.pcr.pick.good.align data/references/silva.v4.align

rm -rf data/references/silva mothur*logfile #garbage collection
