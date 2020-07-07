# Pipeline for matching barcodes and oligos for MPRA data
# output of the Parse task should be used as the input for the ReplicateCount pipeline

workflow MPRAMatch {
  File read_a #R1 fastq
  File read_b #R2 fastq
  File reference_fasta #Oligo sequences with names (can be the oligo order sheet)
  Int read_b_number #2 if you followed the method above
  Int seq_min #Minimum acceptable sequence length when separating the barcodes and oligos
  String working_directory #String of the directory relative to the WDL where the other required scripts live
  String id_out #Project identifier - all files will have this as the prefix for their name
  String barcode_link #6 base sequence on the barcode end of the link between the barcode and oligo - orientation barcode to oligo
  String oligo_link #4 base sequence on the oligo end of the link between the barcode and oligo - orientation barcode to oligo
  String end_oligo_link #4 base sequence at the very end of the oligo

  call Flash { input:
                  read_a=read_a,
                  read_b=read_b,
                  id_out=id_out
                }
  call Pull_Barcodes { input:
                          #pull = pull,
                          working_directory=working_directory,
                          merged_fastq=Flash.out,
                          read_number=read_b_number,
                          id_out=id_out,
                          barcode_link=barcode_link,
                          seq_min=seq_min,
                          oligo_link=oligo_link,
                          end_oligo_link=end_oligo_link
                        }
  call Rearrange { input:
                      matched_barcodes=Pull_Barcodes.out1,
                      id_out=id_out
                    }
  call MiniMap { input:
                    reference_fasta=reference_fasta,
                    organized_fasta=Rearrange.out,
                    id_out=id_out
                  }
  call SAM2MPRA { input:
                      #sam=sam,
                      working_directory=working_directory,
                      sam_file=MiniMap.out1,
                      id_out=id_out
                    }
  call Sort { input:
                  MPRA_out=SAM2MPRA.out,
                  id_out=id_out
                }
  call Ct_Seq { input:
                    #count=count,
                    working_directory=working_directory,
                    sorted=Sort.out,
                    id_out=id_out
                  }
  call Parse { input:
                #  parse=parse,
                  working_directory=working_directory,
                  counted=Ct_Seq.out,
                  id_out=id_out
                }
  call preseq { input:
                 counted=Ct_Seq.out,
                 id_out=id_out
              }

  output {
    File out=Parse.out
     }
  }

task Flash {
  # Flashing raw fastq files together
  File read_a
  File read_b
  String id_out
  command {
    flash2 -r 150 -f 274 -s 20 -o ${id_out}.merged -t 10 ${read_a} ${read_b}
    }
  output {
    File out="${id_out}.merged.extendedFrags.fastq"
    }
  }
task Pull_Barcodes {
  # Pull barcodes from the barcode oligo sequences
  File merged_fastq
  Int read_number
  String working_directory
  String id_out
  String barcode_link
  String oligo_link
  String end_oligo_link
  Int seq_min
  command {
    perl ${working_directory}/pull_barcodes.pl ${merged_fastq} ${read_number} ${id_out}.merged ${barcode_link} ${oligo_link} ${end_oligo_link} ${seq_min}
    }
  output {
    File out1="${id_out}.merged.match"
    File out2="${id_out}.merged.reject"
    }
  }
task Rearrange {
  # Rearrange the ouptut of the pull_barcodes task to be in a fasta format
  File matched_barcodes
  String id_out
  command <<<
    awk '{print ">"$1"#"$3"\n"$4}' ${matched_barcodes} > ${id_out}.merged.match.enh.fa
    >>>
  output {
    File out="${id_out}.merged.match.enh.fa"
    }
  }
task MiniMap {
  # Map the oligos to the reference to get the oligo names
  File reference_fasta
  File organized_fasta
  String id_out
  command {
    minimap2 --for-only -Y --secondary=no -m 10 -n 1 --end-bonus 12 -O 5 -E 1 -k 10 -2K50m --MD --eqx --cs=long -c -a ${reference_fasta} ${organized_fasta} > ${id_out}.merged.match.enh.sam 2> ${id_out}.merged.match.enh.log
    }
  output {
    File out1="${id_out}.merged.match.enh.sam"
    File out2="${id_out}.merged.match.enh.log"
    }
  }
task SAM2MPRA {
  # Convert the output of the SAM file into a format closer to the dictionary. Pulls the CIGAR and error information
  File sam_file
  String working_directory
  String id_out
  command {
    perl ${working_directory}/SAM2MPRA.pl -C ${sam_file} ${id_out}.merged.match.enh.mapped
    }
  output {
    File out="${id_out}.merged.match.enh.mapped"
    }
  }
task Sort {
  File MPRA_out
  String id_out
  command {
    sort -S30G -k2 ${MPRA_out} > ${id_out}.merged.match.enh.mapped.barcode.sort
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.sort"
    }
  }
task Ct_Seq {
  # Counts the number of times a barcode-oligo pair occurs
  File sorted
  String working_directory
  String id_out
  command {
    perl ${working_directory}/Ct_seq.pl ${sorted} 2 4 > ${id_out}.merged.match.enh.mapped.barcode.ct
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct"
    }
  }
task Parse {
  # Parses the barcode oligo pairs to resolve multimapping barcodes
  File counted
  String working_directory
  String id_out
  command {
    perl ${working_directory}/parse_map.pl ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    }
  }

task preseq {
  # Determine sequencing depth
  File counted
  String id_out
  command <<<
    awk '{ct[$4]++}END{for (i in ct)print i "\t" ct[i]}' ${counted} | sort -k1n > ${id_out}.merged.match.enh.mapped.barcode.ct.hist
    preseq lc_extrap -H ${id_out}.merged.match.enh.mapped.barcode.ct.hist -o ${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq -s 25000000 -n 1000 -e 1000000000
    >>>
  output {
    File hist="${id_out}.merged.match.enh.mapped.barcode.ct.hist"
    File res="${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq"
   }
 }
