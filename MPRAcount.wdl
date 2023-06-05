# Pipeline for counting MPRA Replicates
# Requires the same barcode_orientation as MPRAmatch
# Requires the parsed file output from MPRAmatch

workflow MPRAcount {
  Array[File] replicate_fastq #Array of replicate fastq files. Each replicate should have one file.
  Array[String] replicate_id #Identifier for each replicate listed in the replicate fastq list - should be in the same order as replicate fastqs
  Array[Pair[File,String]] fastq_id = zip(replicate_fastq, replicate_id) #Pairing the fastq file to the id
  File parsed #Output of MPRAMatch pipeline
  File acc_id
  Int? barcode_orientation = 2 #2 if MPRAMatch was run with read_a as R1 and read_b as R2, otherwise it should be 1
  Int? bc_len = 20 #If this is not changed a barcode of length 20 will be defaulted.
  String? flags = "-ECSM -A 0.05" #(Error, Cigar, Start/Stop, MD, Error Cutoff)
  String id_out #Overall project id for the final count table
  String working_directory #directory relative to the WDL where the scripts live
  String out_directory #directory relative to the WDL where relevant files will be moved

  scatter (replicate in fastq_id) {
    call prep_counts { input:
                          working_directory=working_directory,
                          sample_fastq=replicate.left,
                          barcode_orientation=barcode_orientation,
                          bc_len=bc_len,
                          sample_id=replicate.right,
                        }
    call associate { input:
                        working_directory=working_directory,
                        matched=prep_counts.out,
                        parsed=parsed,
                        barcode_orientation=barcode_orientation,
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
                            out_directory=out_directory,
                            list_inFile=make_infile.out,
                            flags=flags,
                            id_out=id_out,
                            acc_id=acc_id
                          }
  call count_QC { input:
                    count_out = make_count_table.count,
                    out_directory = out_directory,
                    working_directory = working_directory,
                    id_out = id_out,
                    acc_id = acc_id
                }
  call countRaw { input:
                    count_out = make_count_table.count,
                    cond_out = count_QC.out,
                    id_out = id_out,
                    out_directory = out_directory,
                    working_directory = working_directory
                }
  call relocate { input:
                    matched = prep_counts.out,
                    tag_files = associate.outF,
                    count_out = make_count_table.count,
                    count_log = make_count_table.log,
                    count_stats = make_count_table.stats,
                    cond_out = count_QC.out,
                    out_directory = out_directory
                  }
}

task prep_counts {
  # Grab the barcodes and check that they exist in the dictionary - if they don't exist write them to a seqparate fastq
  File sample_fastq
  Int barcode_orientation
  Int bc_len
  String working_directory
  String sample_id

  command {
    python ${working_directory}/make_counts.py ${sample_fastq} ${sample_id} ${barcode_orientation} ${bc_len}
    }
  output {
    File out="${sample_id}.match"
    }
  }
task associate {
  # Associate the matched barcodes with the associated oligos
  File matched
  File parsed
  Int barcode_orientation
  String working_directory
  String sample_id
  command {
    perl ${working_directory}/associate_tags.pl ${matched} ${parsed} ${sample_id}.tag ${barcode_orientation}
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
  File acc_id
  String working_directory
  String out_directory
  String? flags = ""
  String id_out
  command <<<
    perl ${working_directory}/compile_bc_cs.pl ${flags} ${list_inFile} ${id_out}.count > ${id_out}.log
    awk '{if(NR%7==1){sum=0;good=0;bc=0;over=0;}
      if(NR%7==1){printf "%s\t",$3; printf "%s\t", ${id_out};}
      if(NR%7==3){sum+=3;bc+=$2;over+=$3;}
      if(NR%7==4){printf "%0.f\t", $2; printf "%0.f\t", $3;good+=$3;sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==5){sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==6){sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==0){printf "%0.f\t", sum; printf "%.2f\t", good/(sum)*100; printf "%0.f\t", bc; printf "%0.f\t", over; printf "%.2f\n", good/(over)*100}
      }' ${id_out}.log > ${id_out}.stats
    Rscript ${working_directory}/read_stats.R ${id_out}.stats ${acc_id} ${id_out} ${out_directory}
    >>>
  output {
    File count="${id_out}.count"
    File log="${id_out}.log"
    File stats="${id_out}.stats"
    }
  }
task count_QC {
  File acc_id
  File count_out
  String id_out
  String out_directory
  String working_directory
  command {
    Rscript ${working_directory}/count_QC.R ${acc_id} ${count_out} ${id_out} ${out_directory}
    }
  output {
    File out="${id_out}_condition.txt"
    }
  }
task countRaw {
  File count_out
  File cond_out
  String id_out
  String out_directory
  String working_directory
  command {
    Rscript ${working_directory}/bc_raw.R ${cond_out} ${count_out} ${id_out} ${out_directory}
    }
  }
task relocate {
  Array[File] matched
  Array[File] tag_files
  File count_out
  File count_log
  File count_stats
  File cond_out
  String out_directory
  command <<<
    mv ${sep=' ' matched} ${sep=' ' tag_files} ${count_out} ${count_log} ${count_stats} ${cond_out} ${out_directory}
    >>>
  }
