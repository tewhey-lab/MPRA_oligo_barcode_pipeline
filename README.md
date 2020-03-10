# Tag Analysis Pipeline (two parts)
### Written in the Written Description Language (WDL) version 1.0 more info [here](https://github.com/openwdl/wdl)

## Before running the pipeline
* Have a working version of Java in your workspace
  * Download the installer [here](https://www.java.com/en/download/manual.jsp) and copy to the machine

* Have the latest version of Cromwell and Womtool in your workspace
  * Download [here](https://github.com/broadinstitute/cromwell/releases/tag/48)
  
* Have modules for FLASH2, minimap2, pandas, and B available
  * `conda install -c bioconda flash2 `
  * `conda install -c bioconda minimap2`
  * `conda install -c anaconda pandas`
  * `conda install -c anaconda biopython`

* Make sure all the available scripts (except for the WDL itself) are in a known directory (you will need to provide the path to this directory) in the working directory where the Cromwell and Womtool `.jar`s are present

## Running the WDL
* Validate the file
  `java -jar womtool-48.jar validate <pipeline_name>.wdl`
  
  **NB: use the version number for your version of Womtool downloaded above**

* Generate inputs file
  `java -jar womtool-48.jar inputs <pipeline_name>.wdl > <your_projects_name>_inputs.json`
  
  **NB: see the "Filling in the json" section below for detailed description of input needed**
 
* Run the pipeline with your inputs
  `java -jar cromwell-48.jar run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json`
  
## Filling in the json
A generalized filled in example of each .json is below

_MPRAMatch.wdl_
 ```
     {
       "MPRAMatch.read_a": "full/path/to/read/1.fastq.gz",
       "MPRAMatch.read_b": "full/path/to/read/2.fastq.gz",
       "MPRAMatch.reference_fasta": "/full/path/to/reference/fasta.fa",
       "MPRAMatch.read_b_number": "2",
       "MPRAMatch.seq_min": "100",
       "MPRAMatch.working_directory": "full/path/to/script/location",
       "MPRAMatch.id_out": "Your_Project_ID",
       "MPRAMatch.barcode_link": "6 bases at the barcode end of the sequence linking the barcode and oligo",
       "MPRAMatch.oligo_link": "4 bases at the oligo end of the sequence linking the barcode and oligo",
       "MPRAMatch.end_oligo_link": "4 bases indicating the oligo is no longer being sequenced"

     }
 ```
_ReplicateCount.wd_
 ```
     {
       "ReplicateCount.parsed": "full/path/to/MPRAMatch/output.merged.match.enh.mapped.barcode.ct.parsed",
       "ReplicateCount.read_b_number": "2 (same as used for MPRAMatch)",
       "ReplicateCount.working_directory": "full/path/to/script/location",
       "ReplicateCount.id_out": "Your_Project_ID",
       "ReplicateCount.flags": "-ECSM -A 0.05 (suggested)",
       "ReplicateCount.replicate_fastq": ["full/path/to/celltype1/rep1.fastq.gz", "full/path/to/celltype1/rep2.fastq.gz", "full/path/to/celltype1/rep3.fastq.gz", "full/path/to/celltype2/rep1.fastq.gz",...],
       "ReplicateCount.replicate_id": ["Celltype1_r1", "Celltype1_r2", "Celltype1_r3", "Celltype2_r1", ...]

     }
 ```

**NB: For the replicate fastq files, a single fastq.gz file is required for each replicate, if there are more than one they should be concatenated together**

## Outputs Needed at later steps
It is suggested that you note the job id generated within cromwell for assistance finding these files at a later date.

The output file from MPRAMatch needed as input for the ReplicateCount pipeline can be found at:
  * Parsed File    : `cromwell-exectutions/MPRAMatch/<job_id>/call-Parse/execution/<id_out>.merged.match.enh.mapped.barcode.ct.parsed`
The following file can then be input into the R pipeline [here](https://github.com/tewhey-lab/MPRA_tag_analysis) for analysis:
  * Count File     : `cromwell-exectutions/ReplicateCount/<job_id>/call-make_count_table/execution/<id_out>.count`
