#!/usr/bin/env bash

TREE=$1
COUNT_TABLE=$2
PROCESSORS=$3

# TREE=data/crc/data.10439.all_1.raw.tre
# COUNT_TABLE=data/crc/data.10439.all_1.raw.count_table
# PROCESSORS=16

mothur "#unifrac.weighted(tree = $TREE, count = $COUNT_TABLE, processors = $PROCESSORS, distance = square)"

MOTHUR_FILE=$TREE"1.weighted.phylip.dist"
DISTANCE=`echo $TREE | sed -e "s/tre/wunifrac.dist/"`

HEADER=`cut -f 1 $MOTHUR_FILE | tail -n +2 | paste -s`
echo sample$'\t'$HEADER | sed -E 's/[[:space:]]+/\t/g' > $DISTANCE
tail -n +2 $MOTHUR_FILE | sed -E 's/[[:space:]]+/\t/g' >> $DISTANCE
gzip $DISTANCE

rm $MOTHUR_FILE
rm $TREE"1.wsummary"
