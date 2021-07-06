# MPRA oligo/barcode reconstruction and barcode counting pipeline (two parts)
### Written in the Written Description Language (WDL) version 1.0 more info [here](https://github.com/openwdl/wdl)

## Before running the pipeline
* Have the latest version of Cromwell and Womtool in your workspace
  * `conda install -c bioconda cromwell`
  * `conda install -c bioconda womtool`

* Have modules for FLASH2, minimap2 (version 2.17), preseq, pandas, reshape2, ggplot2, gridextra, and Biopython available
  * `conda install -c bioconda  flash2 minimap2=2.17 preseq pandas biopython r-ggplot2`
  * `conda install -c conda-forge r-reshape2`
  * `conda install -c r r-gridextra`

* Make sure all the available scripts (except for the WDL itself) are in a known directory (you will need to provide the path to this directory)

## Running the WDL
* Validate the file
  `womtool validate <pipeline_name>.wdl`

  **NB: use the version number for your version of Womtool downloaded above**

* Generate inputs file
  `womtool inputs <pipeline_name>.wdl > <your_projects_name>_inputs.json`

  **NB: see the "Filling in the json" section below for detailed description of input needed**

* Run the pipeline with your inputs
  `cromwell run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json`

**To submit to slurm** Make sure that you give the pipeline enough memory to run, if the pipeline fails the first time you run it, look at the end of the slurm output file to determine whether you need to give it more time or more memory
  * `sbatch -p compute -q batch -t 24:00:00 --mem=45GB -c 8 --wrap "cromwell run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json"`

## Filling in the json
A generalized filled in example of each .json is below

_MPRAmatch.wdl_
 ```
     {
       "MPRAmatch.read_a": "full/path/to/read/1.fastq.gz",
       "MPRAmatch.read_b": "full/path/to/read/2.fastq.gz",
       "MPRAmatch.reference_fasta": "/full/path/to/reference/fasta.fa",
       "MPRAmatch.read_b_number": "2",
       "MPRAmatch.seq_min": "100",
       "MPRAmatch.working_directory": "full/path/to/script/location",
       "MPRAmatch.out_directory": "full/path/to/output/directory/"
       "MPRAmatch.id_out": "Your_Project_ID",
       "MPRAmatch.barcode_link": "6 bases at the barcode end of the sequence linking the barcode and oligo",
       "MPRAmatch.oligo_link": "4 bases at the oligo end of the sequence linking the barcode and oligo",
       "MPRAmatch.end_oligo_link": "4 bases indicating the oligo is no longer being sequenced"

     }
 ```
_MPRAcount.wdl_
 ```
     {
       "MPRAcount.parsed": "full/path/to/MPRAmatch/output.merged.match.enh.mapped.barcode.ct.parsed",
       "MPRAcount.read_b_number": "2 (same as used for MPRAmatch)",
       "MPRAcount.working_directory": "full/path/to/script/location",
       "MPRAcount.out_directory": "full/path/to/output/directory/",
       "MPRAcount.id_out": "Your_Project_ID",
       "MPRAcount.flags": "-ECSM -A 0.05 (suggested)",
       "MPRAcount.replicate_fastq": ["full/path/to/celltype1/rep1.fastq.gz", "full/path/to/celltype1/rep2.fastq.gz", "full/path/to/celltype1/rep3.fastq.gz", "full/path/to/celltype2/rep1.fastq.gz",...],
       "MPRAcount.replicate_id": ["Celltype1_r1", "Celltype1_r2", "Celltype1_r3", "Celltype2_r1", ...]

     }
 ```

**NB: For the replicate fastq files, a single fastq.gz file is required for each replicate, if there are more than one they should be concatenated together**

## Outputs Needed at later steps
It is suggested that you note the job id generated within cromwell for assistance finding these files at a later date.

The output file from MPRAmatch needed as input for the MPRAcount pipeline can be found at:
  * Parsed File      : `/specified/output/directory/<id_out>.merged.match.enh.mapped.barcode.ct.parsed`

The following files can then be input into the R pipeline [here](https://github.com/tewhey-lab/MPRA_tag_analysis) for analysis:
  * Count File       : `cromwell-executions/MPRAcount/<job_id>/call-make_count_table/execution/<id_out>.count`
  * Attributes Table : The output of `make_project_list.pl` and `make_attributes_oligo.pl`. If you want to bypass the use of `make_project_list.pl`, you can pass `make_attributes_oligo.pl` a tab delimited file with two columns, the first column should be the oligo names and the second column should be the project(s) that the associated oligo belongs to. If an oligo belongs to multiple projects they should be separated by commas.

## How the Pipelines work

_MPRAmatch_

![Graphical Pipeline](graphics/MPRAmatch_pipeline.svg?raw=true "MPRAmatch Graphical Pipeline")

The above image is a graphical representation of the MPRAmatch pipeline. Green objects represent files and information provided to the pipeline which are passed directly to a script or program, blue objects are the calls to the modules called for above, yellow objects refer to scripts written for the pipeline, and the barcode-oligo dictionary is in red.

The two fastq files from the initial barcode-oligo sequencing are fed into FLASH2 in order to merge them into a single fastq. The merged fastq is then passed to a script which pulls the barcode and oligo sequences for each record in the fastq based on the linker sequences between the barcode and oligo, and at the end of the oligo. The barcode/oligo pair information is rearranged into a FASTA format and passed to MiniMap2 along with the reference fasta. The resulting SAM file is parsed for the Oligo name, barcode sequence, CIGAR, and error information for each mapped record. The number of times each barcode appears for each oligo is then counted; the output is passed to preseq to determine sequencing depth, and parsed to resolve barcodes which map to multiple oligos.

_MPRAcount_

![Graphical Pipeline](graphics/MPRAcount_pipeline.svg?raw=true "MPRAcount Graphical Pipeline")

The above image is a graphical representation of the MPRAcount pipeline. The tiled green object represents arrays of files and information passed to the pipeline, the turquoise object represents the output of the MPRAmatch pipeline, yellow objects refer to scripts written for the pipeline, and the final count table is in red.

The array of fastq files are processed in parallel and the barcodes, along with the barcode-oligo dictionary, are passed to a script which associates the barcodes with the matching oligo. These files are then passed to a script which organises the matched barcodes and oligos into a count table.
