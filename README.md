# Tag Analysis Pipeline
### Written in the Written Description Language (WDL) version 1.0 more info [here](https://github.com/openwdl/wdl)

## Before running the pipeline
* Have a working version of Java in your workspace
  * Download the installer [here](https://www.java.com/en/download/manual.jsp) and copy to the machine

* Have the latest version of Cromwell and Womtool in your workspace
  * Download [here](https://github.com/broadinstitute/cromwell/releases/tag/48)
  
* Have modules for FLASH2 and minimap2 available
  * `conda install -c bioconda flash2 `
  * `conda install -c bioconda minimap2`

* Make sure all the available scripts (except for the WDL itself) are in a folder called scripts in the working directory where the Cromwell and Womtool `.jar`s are present

## Running the WDL
* Validate the file
  `java -jar womtool-48.jar validate count_pipeline_2.0.`
  **NB: use the version number for your version of Womtool downloaded above**

* Generate inputs file
  `java -jar womtool-48.jar inputs count_pipeiline_2.0.wdl > <your_projects_name>_inputs.json`
  **NB: see the "Filling in the json" section below for detailed description of input needed**
 
* Run the pipeline with your inputs
  `java -jar cromwell-48.jar run count_pipeline_2.0.wdl --inputs <your_projects_name>_inputs.json`
  
## Filling in the json
The inputs file will look like this:
  ```
     {
       "MPRACount.seq_min": "Int",
       "MPRACount.id_out": "String",
       "MPRACount.flags": "String",
       "MPRACount.end_oligo_link": "String",
       "MPRACount.oligo_type": "File",
       "MPRACount.oligo_link": "String",
       "MPRACount.replicate_fastq": "Array[File]",
       "MPRACount.associate_tags": "File",
       "MPRACount.read_number": "Int",
       "MPRACount.replicate_id": "Array[String]",
       "MPRACount.read_b": "File",
       "MPRACount.read_a": "File",
       "MPRACount.barcode_link": "String",
       "MPRACount.reference_fasta": "File",
     }
 ```

Below is a description of what kind of file should be input in the .json file above. Script names (italicized in the description table) are the same as those provided in the repository.


**File** | **Description of File**
-------- | -----------------------
MPRACount.seq_min         | Minimum acceptable sequence length for barcode-oligo sequence
MPRACount.id_out          | Project Name to be used for all files
MPRACount.flags           | Flags to be used with the _compile_bc.pl_ script, detailed within the script (Suggested Use: "-ECSM -A 0.05")
MPRACount.end_oligo_link  | 4 base sequence to indicate the end of the oligo
MPRACount.oligo_type      | _make_attributes_oligo_v*.pl_ Here use the script that is associated with the type of oligo you have, you only need to give the name of the script within the `/scripts/` folder
MPRACount.oligo_link      | 4 base sequence on the oligo end of the linker sequence between oligo and barcode
MPRACount.replicate_fastq | Array of the fastq files for each replicate. Each replicate should have 1 composite file.
MPRACount.read_number     | 1 if the reverse complement needs to be taken, otherwise 2
MPRACount.replicate_id    | Array of the replicate ids, to be used as column headers in count table
MPRACount.read_b          | fastq file to be flashed
MPRACount.read_a          | fastq file to be flashed
MPRACount.barcode_link    | 6 base sequence on the barcode end of the linker sequence between oligo and barcode
MPRACount.reference_fasta | Reference fasta of all oligos in experiment


## Outputs
The following files can then be input into the R pipeline [here](https://github.com/tewhey-lab/MPRA_tag_analysis) for analysis:
  * Attributes File: `cromwell-exectutions/MPRACount/<job_id>/call-make_attribute_file/<id_out>.attributes`
  * Count File     : `cromwell-exectutions/MPRACount/<job_id>/call-make_count_table/<id_out>.count`
