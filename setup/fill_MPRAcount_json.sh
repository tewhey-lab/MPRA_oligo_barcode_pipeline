#!/bin/bash

#Assign variables to be built into the inputs json

PROJ="<project_name>"
FLOC="<path/to/general/project/folder/${PROJ}"
OUT="${PROJ}_<suffix>" # Suffix should match the one from fill_MPRAmatch.json
ACC="${FLOC}/setup/<path/to/replicate_index>"
ID=`awk '{print $2}' ${ACC} | uniq | awk '{aggr=aggr",\""$1"\""} END {print aggr}' | sed 's/,//'`
REPS=`awk '{print $2}' ${ACC} | uniq | awk '{aggr=aggr",\"'${FLOC}'/fastq/"$1".fastq.gz\""} END {print aggr}' | sed 's/,//'`
PARS="${FLOC}/<MPRAmatch_output_dir>/${OUT}.merged.match.enh.mapped.barcode.ct.parsed" #adjust directory to match directory organization

cp <path/to/cloned/repository>/inputs/MPRAcount_sample_inputs.json ../submission/MPRAcount_${PROJ}_inputs.json #If using different organization change from submission to desired location

jq --arg PROJ ${PROJ} --arg FLOC ${FLOC} --arg ACC ${ACC} --arg OUT ${OUT} --arg ID ${ID} --arg REPS ${REPS} --arg PARS ${PARS} -M '. + {"MPRAcount.out_directory":'\"${FLOC}/\"',
    "MPRAcount.id_out":'\"$OUT\"', "MPRAcount.parsed":'\"$PARS\"', "MPRAcount.acc_id":'\"$ACC\"', 
    "MPRAcount.replicate_fastq":'\[$REPS\]', "MPRAcount.replicate_id":'\[$ID\]'}' ../submission/MPRAcount_${PROJ}_inputs.json > <path/to/cloned/repository>/inputs/MPRAcount/MPRAcount_${PROJ}_inputs.json #change submission to match above
