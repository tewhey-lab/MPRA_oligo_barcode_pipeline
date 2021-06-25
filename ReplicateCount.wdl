# Pipeline for counting MPRA Replicates
# Requires the same read_b_number as MPRAMatch
# Requires the parsed file output from MPRAMatch

workflow ReplicateCount {
  Array[File] replicate_fastq #Array of replicate fastq files. Each replicate should have one file.
  Array[String] replicate_id #Identifier for each replicate listed in the replicate fastq list - should be in the same order as replicate fastqs
  Array[Pair[File,String]] fastq_id = zip(replicate_fastq, replicate_id) #Pairing the fastq file to the id
  File parsed #Output of MPRAMatch pipeline
  Int read_b_number #2 if MPRAMatch was run with read_a as R1 and read_b as R2, otherwise it should be 1
  String flags #-ECSM -A 0.05 (Error, Cigar, Start/Stop, MD, Error Cutoff)
  String id_out #Overall project id for the final count table
  String working_directory #directory relative to the WDL where the scripts live
  String out_directory #directory relative to the WDL where relevant files will be moved

  scatter (replicate in fastq_id) {
    call prep_counts { input:
                          working_directory=working_directory,
                          sample_fastq=replicate.left,
                          read_b_number=read_b_number,
                          parsed=parsed,
                          sample_id=replicate.right
                        }
    call associate { input:
                        working_directory=working_directory,
                        matched=prep_counts.out,
                        parsed=parsed,
                        read_b_number=read_b_number,
                        sample_id=replicate.right
                      }
                    }
  call make_infile { input:
                        working_directory=working_directory,
                        tag_files=associate.outF,
                        tag_ids=associate.outS,
                        id_out=id_out
                      }
  call make_count_table { input:
                            working_directory=working_directory,
                            list_inFile=make_infile.out,
                            flags=flags,
                            id_out=id_out
                          }
  call relocate { input:
                    matched = prep_counts.out,
                    tag_files = associate.outF,
                    count_out = make_count_table.out,
                    out_directory = out_directory
                  }
}

task prep_counts {
  # Grab the barcodes and check that they exist in the dictionary - if they don't exist write them to a seqparate fastq
  File sample_fastq
  File parsed
  Int read_b_number
  String working_directory
  String sample_id
  command {
    python ${working_directory}/make_counts.py ${sample_fastq} ${parsed} ${sample_id} ${read_b_number}
    }
  output {
    File out="${sample_id}.match"
    }
  }
task associate {
  # Associate the matched barcodes with the associated oligos
  File matched
  File parsed
  Int read_b_number
  String working_directory
  String sample_id
  command {
    perl ${working_directory}/associate_tags.pl ${matched} ${parsed} ${sample_id}.tag ${read_b_number}
    }
  output {
    File outF="${sample_id}.tag"
    String outS="${sample_id}"
    }
  }
task make_infile {
  # make a list of the association output and the tag ids to pass to the barcode compilation function
  Array[File] tag_files
  Array[String] tag_ids
  String working_directory
  String id_out
  command <<<
    python ${working_directory}/make_infile.py ${sep=',' tag_ids} ${sep=',' tag_files} ${id_out}
  >>>
  output {
    File out="${id_out}_samples.txt"
    }
  }
task make_count_table {
  # Compile barcodes into a count table - columns from left to right: barcode oligo (Error CIGAR MD Aln_Start:Stop) [replicate names]
  File list_inFile
  String working_directory
  String? flags = ""
  String id_out
  command {
    perl ${working_directory}/compile_bc.pl ${flags} ${list_inFile} ${id_out}.count
    }
  output {
    File out="${id_out}.count"
    }
  }
task relocate {
  Array[File] matched
  Array[File] tag_files
  File count_out
  String out_directory
  command {
    mv ${matched} ${tag_files} ${count_out} ${out_directory}
    }
  }
