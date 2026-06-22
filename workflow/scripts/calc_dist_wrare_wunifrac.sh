#!/usr/bin/env bash

TREE=$1
COUNT_TABLE=$2
ITERATIONS=$3
PROCESSORS=$4

# TREE=data/crc/data.1000.35_1.raw.tre
# COUNT_TABLE=data/crc/data.1000.35_1.raw.count_table
# ITERATIONS=1
# ITERATIONS=5
# PROCESSORS=16

TEMP_TREE=${TREE}${ITERATIONS}
cp $TREE $TEMP_TREE

mothur "#unifrac.weighted(tree = $TEMP_TREE, count = $COUNT_TABLE, processors = $PROCESSORS, iters = $ITERATIONS, subsample = TRUE, distance = square)"

MOTHUR_FILE=${TREE}${ITERATIONS}1.weighted.ave.dist

if [ $ITERATIONS -eq 1 ]; then
  DISTANCE=`echo "$TREE" | sed -e "s/raw.tre/sub.wunifrac.dist/"`
elif [ $ITERATIONS -gt 1 ]; then
  DISTANCE=`echo "$TREE" | sed -e "s/raw.tre/rare.wunifrac.dist/"`
fi

HEADER=`cut -f 1 $MOTHUR_FILE | tail -n +2 | paste -s`
echo sample$'\t'$HEADER | sed -E 's/[[:space:]]+/\t/g' > $DISTANCE
tail -n +2 $MOTHUR_FILE | sed -E 's/[[:space:]]+/\t/g' >> $DISTANCE
gzip -f $DISTANCE

rm ${TREE}${ITERATIONS}1.weighted.ave.dist
rm ${TREE}${ITERATIONS}1.weighted.phylip.dist
rm ${TREE}${ITERATIONS}1.weighted.std.dist
rm ${TREE}${ITERATIONS}1.wsummary
rm $TEMP_TREE
