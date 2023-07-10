# MPRA oligo/barcode reconstruction and barcode counting pipeline (two parts)
### Written in the Written Description Language (WDL) version 1.0 more info [here](https://github.com/openwdl/wdl)

## Before running the pipeline
We have the whole environment containerized, the definition file is available in the `environment/` folder of this repository.

If you are unable to run the pipeline via a container, then set up your environment as described below:

* Have modules for Cromwell, Womtool, FLASH2, minimap2 (version 2.17), preseq, pandas, reshape2, ggplot2, gridextra, and Biopython available (`.yml` of this conda enviornment can be found in the environment tab)
  * `conda install -c bioconda cromwell womtool flash2 minimap2=2.17 preseq pandas samtools`
  * `conda install -c conda-forge r-reshape2 r-ggplot2 r-gridextra biopython`

* Make sure all the available scripts (except for the WDL itself) are in a known directory (you will need to provide the path to this directory)

* WDL does not get rid of intermediate files. The  pipelines are set up to relocate files that are important for later use from where the pipeline is run to a more permanent location. Consider running the pipeline in a scratch area so you don't have to go and delete other intermediate files after the pipeline completes itself. If you do opt to delete the files manually, please check that the relocation at the end of the pipeline has completed.

## Suggested directory organization chart:
```
	- overall_project/
		- specific_project/
			- MPRAmatch_output/
			- MPRAcount_output/
			- MPRAmodel_output/
			- setup/ 	(copy from cloned repository)
			- fastq/
			- submission/
		- cloned_repository/
      - environment/
			- graphics/
			- inputs/
			- scripts/
			- setup/
			- MPRAcount.wdl
			- MPRAmatch.wdl
			- README.md
```
## Running the WDL

  * Validate the file
    `womtool validate <pipeline_name>.wdl`

    **NB: use the version number for your version of Womtool downloaded above, but you should not need to validate the pipelines included here**

  * Generate inputs file
    `womtool inputs <pipeline_name>.wdl > <your_projects_name>_inputs.json`

    **NB: see the "Filling in the json" section below for detailed description of input needed, scripts and sample inputs are provided in this repository to make this step easier**

  * Run the pipeline with your inputs
    `cromwell run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json`

**To submit to slurm** Make sure that you give the pipeline enough memory to run, if the pipeline fails the first time you run it, look at the end of the slurm output file to determine whether you need to give it more time or more memory
  * `sbatch -p compute -q batch -t 24:00:00 --mem=45GB -c 8 --wrap "cromwell run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json"` <br>

  * **OR** you can use the runscript and submission template below which utilizes the singularity container: <br>
  Runscript:
  ```
  echo "Running Cromwell"

  cromwell run /path/to/cloned_repository/MPRAcount.wdl --inputs path/to/cloned_repository/inputs/MPRAcount/MPRAcount_<your_project>_inputs.json

  echo "Finished Cromwell"
  ```
  Submission template (for a SLURM based scheduler):
  ```
  #!/bin/bash
  #SBATCH --job-name= <your_projects_name>
  #SBATCH -p compute # partition(this is the standard)
  #SBATCH -q batch
  #SBATCH -N 1 # number of nodes
  #SBATCH -n 45 # number of cores
  #SBATCH --mem 200GB # memory pool for all cores
  #SBATCH -t 3-00:00 # time (D-HH:MM)
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user= <your_email_here>

  echo "Loading Singularity Module"

  module load singularity

  echo "Executing SIF with Code"

  singularity exec /path/to/your/built/container sh /path/to/runscript/MPRAcount_call.sh

  echo "Done"
  ```

## Filling in the json
A generalized filled in example of each .json is below

_MPRAmatch.wdl_
 ```
     {
       "MPRAmatch.read_a": "/full/path/to/read/1.fastq.gz",
       "MPRAmatch.read_b": "/full/path/to/read/2.fastq.gz",
       "MPRAmatch.reference_fasta": "/full/path/to/reference/fasta.fa",
       "MPRAmatch.working_directory": "/full/path/to/MPRA_oligo_barcode_pipeline/scripts",
       "MPRAmatch.out_directory": "/full/path/to/output/directory/"
       "MPRAmatch.id_out": "Your_Project_ID"
     }
 ```

 **NB: If the sequencing for read1 and read2 results in multiple files for both (i.e. multiple lanes or split between plates) when concatenating them into a single read1 file and read2 file make sure paired files are in the same order between both fastqs.**
 
_MPRAcount.wdl_
 ```
     {
       "MPRAcount.parsed": "/full/path/to/MPRAmatch/output.merged.match.enh.mapped.barcode.ct.parsed",
       "MPRAcount.acc_id": "File (tab separated file mapping accession numbers to replicates and associated cell type)",
       "MPRAcount.working_directory": "/full/path/to/MPRA_oligo_barcode_pipeline/scripts",
       "MPRAcount.out_directory": "/full/path/to/output/directory/",
       "MPRAcount.id_out": "Your_Project_ID",
       "MPRAcount.replicate_fastq": ["/full/path/to/celltype1/rep1.fastq.gz", "/full/path/to/celltype1/rep2.fastq.gz", "/full/path/to/celltype1/rep3.fastq.gz", "/full/path/to/celltype2/rep1.fastq.gz",...],
       "MPRAcount.replicate_id": ["Celltype1_r1", "Celltype1_r2", "Celltype1_r3", "Celltype2_r1", ...]

     }
 ```

**NB: For the reference fasta, make sure it contains _only_ the CRE sequences and not the flanking sequence. For the replicate fastq files, a single fastq.gz file is required for each replicate, if there are more than one they should be concatenated together. This can be done quickly using the `merge_*.sh` files in the `setup` folder**

There are several optional inputs designed to make use easier.
 - `MPRAmatch.barcode_orientation` / `MPRAcount.barcode_orientation` : Integer, default to 2. Pulls barcode 5' to 3' end. If set to 1 barcode will be pulled in 3' to 5' orientation.
 - `MPRAmatch.thread` : Integer, default to 30. Number of threads to be passed to FLASH2 and MiniMap2. Each gets this number of threads.
 - `MPRAmatch.mem` : Integer, default to 30. Amount of memory (in GB) to be used when sorting.
 - `MPRAmatch.read_len` : Integer, default to 250. Maximum length (in bp) of reads to be flashed
 - `MPRAmatch.seq_min` : Integer, default to 100. Minimum sequence length to pull for barcode
 - `MPRAmatch.enh_min` : Integer, default to 50. Minimum enhancer length to pull
 - `MPRAmatch.enh_max` : Integer, default to 210. Maximum enhancer length to pull
 - `MPRAcount.bc_len` : Integer, default to 20. Length of barcodes to be pulled from replicate files
 - `MPRAcount.flags` : String, default to "-ECSM -A 0.05" Any combination of these flags or none can be used.
 - `MPRAmatch.barcode_link` : String, default to "TCTAGA". 6 bases at the barcode end of the sequence linking the barcode and oligo
 - `MPRAmatch.oligo_link` : String, default to "AGTG". 4 bases at the oligo end of the sequence linking the barcode and oligo
 - `MPRAmatch.end_oligo_link` : String, default to "CGTC". 4 bases indicating the oligo is no longer being sequenced


It is suggested to use `fill_<pipeline>_json.sh` scripts to help fill in the jsons above instead of trying to fill them out manually.

## Outputs Needed at later steps
It is suggested that you note the job id generated within cromwell for assistance finding these files at a later date.

The output file from `MPRAmatch` needed as input for the `MPRAcount` pipeline can be found at:
  * Parsed File      : `/specified/output/directory/<id_out>.merged.match.enh.mapped.barcode.ct.parsed`

The following files can then be input into the R pipeline [here](https://github.com/tewhey-lab/MPRA_tag_analysis) for analysis:
  * Count File       : `/specified/output/directory/<id_out>.count`
  * Attributes Table : The output of `make_project_list.pl` and `make_attributes_oligo.pl`. If you want to bypass the use of `make_project_list.pl`, you can pass `make_attributes_oligo.pl` a tab delimited file with two columns, the first column should be the oligo names and the second column should be the project(s) that the associated oligo belongs to. If an oligo belongs to multiple projects they should be separated by commas.
  * Condition file   : `/specified/output/directory/<id_out>_condition.txt`

## How the Pipelines work

_MPRAmatch_

![Graphical Pipeline](graphics/MPRAmatch_pipeline.svg?raw=true "MPRAmatch Graphical Pipeline")

The above image is a graphical representation of the MPRAmatch pipeline. Green objects represent files and information provided to the pipeline which are passed directly to a script or program, blue objects are the calls to the modules called for above, yellow objects refer to scripts written for the pipeline, and the barcode-oligo dictionary is in red.

The two fastq files from the initial barcode-oligo sequencing are fed into FLASH2 in order to merge them into a single fastq. The merged fastq is then passed to a script which pulls the barcode and oligo sequences for each record in the fastq based on the linker sequences between the barcode and oligo, and at the end of the oligo. The barcode/oligo pair information is rearranged into a FASTA format and passed to MiniMap2 along with the reference fasta. The resulting SAM file is parsed for the Oligo name, barcode sequence, CIGAR, and error information for each mapped record. The number of times each barcode appears for each oligo is then counted; the output is passed to preseq to determine sequencing depth, and parsed to resolve barcodes which map to multiple oligos.

_MPRAcount_

![Graphical Pipeline](graphics/MPRAcount_pipeline.svg?raw=true "MPRAcount Graphical Pipeline")

The above image is a graphical representation of the MPRAcount pipeline. The tiled green object represents arrays of files and information passed to the pipeline, the turquoise object represents the output of the MPRAmatch pipeline, yellow objects refer to scripts written for the pipeline, and the final count table is in red.

The array of fastq files are processed in parallel and the barcodes, along with the barcode-oligo dictionary, are passed to a script which associates the barcodes with the matching oligo. These files are then passed to a script which organises the matched barcodes and oligos into a count table.
