# Tag Analysis Pipeline
### Written in the Written Description Language (WDL) version 1.0 more info [here](https://github.com/openwdl/wdl)

## Before running the pipeline
* Have a working version of Java in your workspace
  * `conda install -c cyclus java-jdk`
  * Or download the installer [here](https://www.java.com/en/download/manual.jsp) and copy to the machine

* Have the latest version of Cromwell and Womtool in your workspace
  * Download [here](https://github.com/broadinstitute/cromwell/releases/tag/48)
  
* Have modules for FLASH2 and minimap2 available
  * `conda install -c bioconda flash2 `
  * `conda install -c bioconda minimap2`

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
  ```{
       "MPRACount.count": "File",
       "MPRACount.parse": "File",
       "MPRACount.seq_min": "Int",
       "MPRACount.id_out": "String",
       "MPRACount.sam_convert": "File",
       "MPRACount.end_oligo_link": "String",
       "MPRACount.oligo_link": "String",
       "MPRACount.read_number": "Int",
       "MPRACount.read_b": "File",
       "MPRACount.read_a": "File",
       "MPRACount.barcode_link": "String",
       "MPRACount.reference_fasta": "File",
       "MPRACount.pull": "File"
     }```
     
  * Description of files needed (.pl files are available in this repository):
    * MPRACount.count           : Ct_seq.pl
    * MPRACount.parse           : parse_map.pl
    * MPRACount.seq_min         : Minimum acceptable sequence length for barcode-oligo sequence
    * MPRACount.id_out          : Project Name to be used for all files
    * MPRACount.sam_convert     : SAM2MPRA.pl
    * MPRACount.end_oligo_link  : 4 base sequence to indicate the end of the oligo
    * MPRACount.oligo_link      : 4 base sequence on the oligo end of the linker sequence between oligo and barcode
    * MPRACount.read_number     : 1 if the reverse complement needs to be taken, otherwise 2
    * MPRACount.read_b          : fastq file to be flashed
    * MPRACount.read_a          : fastq file to be flashed
    * MPRACount.barcode_link    : 6 base sequence on the barcode end of the linker sequence between oligo and barcode
    * MPRACount.reference_fasta : Reference fasta of all oligos in experiment
    * MPRACount.pull            : pull_barcodes.pl
