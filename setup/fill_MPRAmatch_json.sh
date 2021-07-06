#!/bin/bash

#Assign variables to be built into the inputs json

PROJ="<project_name>"
FLOC="<path/to/general/project/folder/${PROJ}"
OUT="${PROJ}_<suffix>" 
READ1="${FLOC}/fastq/${PROJ}_read1.fastq.gz"
READ2="${FLOC}/fastq/${PROJ}_read2.fastq.gz"
FASTA="${FLOC}/setup/${PROJ}_reference.fasta.gz"

cp <path/to/cloned/repository>/inputs/MPRAmatch_sample_inputs.json ../submission/MPRAmatch_${PROJ}_inputs.json #If using different organization change from submission to desired location

jq --arg PROJ ${PROJ} --arg FLOC ${FLOC} --arg OUT ${OUT} --arg READ1 ${READ1} --arg READ2 ${READ2} --arg FASTA ${FASTA} -M '. + {"MPRAmatch.read_a":'\"${READ1}/\"',
    "MPRAmatch.read_b":'\"${READ2}/\"', "MPRAmatch.reference_fasta":'\"${FASTA}/\"',
    "MPRAmatch.id_out":'\"${OUT}/\"',"MPRAmatch.out_directory":'\"${FLOC}/\"'}' ../submission/MPRAmatch_${PROJ}_inputs.json > <path/to/cloned/repository>/inputs/MPRAmatch/MPRAmatch_${PROJ}_inputs.json #change submission to match above
