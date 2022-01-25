#!/bin/bash

# Use this pipeline to merge replicated encode files.
# acc_id.txt should be a tab separated file with the first column containing accession numbers
    # The second column should be the relative concatenated ID
    # i.e. <accession1>   read1
    #      <accession2>   read1
  # This file will be used in filling out the pipeline inputs

for i in `awk '{print $2}' acc_id.txt | sort | uniq`
do

echo $i

CMD=`awk -v i="$i" '$2==i {aggr=aggr " " $1".fastq.gz"} END {print aggr}' acc_id.txt`

echo $CMD

cat $CMD > ../fastq/$i.fastq.gz


done
