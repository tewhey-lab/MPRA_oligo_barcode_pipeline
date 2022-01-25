## Input setup for MPRAmatch and MPRAcount

### `fill_MPRAmatch_json.sh`

This script provides an easy method to fill out the input JSON needed for `MPRAmatch` to run. Below is an explanation of the variables to fill out to make it work.

  - `PROJ`  : "Your project name"
  - `FLOC`  : "<path/to/general/project/folder/${PROJ}>" (Should not contain a "/" at the end)
  - `OUT`   : "${PROJ}_\<suffix>" (The prefix for all files that will be written)
  - `READ1` : "${FLOC}/fastq/${PROJ}_read1.fastq.gz" (Pointer to the read1 fastq file)
  - `READ2` : "${FLOC}/fastq/${PROJ}_read2.fastq.gz" (Pointer to the read2 fastq file)
  - `FASTA` : "${FLOC}/setup/${PROJ}_reference.fasta.gz" (Pointer the the reference fasta file)

**NB** You will also need to input the location of the sample file and the destination file within the script

### `fill_MPRAcount_json.sh`

This script provides an easy method to fill out the input JSON needed for `MPRAcount` to run. Below is an explanation of the variables to fill out to make it work.

  - `PROJ`  : "Your project name"
  - `FLOC`  : "<path/to/general/project/folder/${PROJ}>" (Should not contain a "/" at the end)
  - `OUT`   : "${PROJ}_\<suffix>" (This should match what was used for `fill_MPRAmatch_json.sh`)
  - `ACC`   : "${FLOC}/setup/<path/to/replicate_index>" (Tab separated file with 4 columns, example [here](https://github.com/tewhey-lab/MPRA_oligo_barcode_pipeline/blob/master/setup/acc_id_example.txt))
  - `PARS`  : "${FLOC}/<MPRAmatch_output_dir>/${OUT}.merged.match.enh.mapped.barcode.ct.parsed" (Pointer to the parsed file output from MPRAmatch)

**NB** You will also need to input the location of the sample file and the destination file within the script

### Set up the `acc_id` file

This file consists of 4 tab separated columns further described (in order) below.

  1. **Accession** This is the accession number if pulling from ENCODE, but it can also be an internal accession number, a dummy variable, or the sequencing ID number if you don't rename the files earlier on. These should all be unique.
  2. **ID** This is an identifier as to what the file contains. i.e. If the file represents part or all of K562 Replicate 1 you might use "k562_r1", this can be repeated for as many files fit the description
  3. **Cell-type** What cell-type goes with this file, this can be repeated for as many files as necessary.
  4. **DNA/RNA** Is the sample DNA or RNA. Please use "DNA" and "RNA" for this column and not a flag.

This file is used to build the condition table that is input into the [`MPRAmodel`](https://github.com/tewhey-lab/MPRAmodel) pipeline, so please bear that in mind when creating your `acc_id` file.
