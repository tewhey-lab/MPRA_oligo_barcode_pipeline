# Pipeline for matching barcodes and oligos for MPRA data
# output of the Parse task should be used as the input for the ReplicateCount pipeline

workflow MPRAmatch {
  File read_a #R1 fastq
  File read_b #R2 fastq
  File reference_fasta #Oligo sequences with names (can be the oligo order sheet)
  File? attributes #Optional attributes file used for saturation mutagenesis libraries
  Int? barcode_orientation = 2 #2 if you followed the method above, otherwise 1. Default to 2
  Int? thread = 30 #Number of threads to be passed to FLASH2 and MiniMap2. Default to 30
  Int? mem = 30 #Memory to be passed to the sort function. Default to 30G
  Int? read_len = 250 #Length of reads that are being flashed. If mixed lengths use max.
  Int? seq_min = 100 #Minimum acceptable sequence length when separating the barcodes and oligos
  Int? enh_min = 50 #Minimum acceptable length for an oligo
  Int? enh_max = 210 #Maximum acceptable length for an oligo
  Int? bc_len = 20 #Length of barcodes used for project
  String working_directory #String of the directory relative to the WDL where the other required scripts live
  String out_directory #String of the directory that all files will be copied to
  String id_out #Project identifier - all files will have this as the prefix for their name
  String? barcode_link = "TCTAGA" #6 base sequence on the barcode end of the link between the barcode and oligo - orientation barcode to oligo
  String? oligo_link = "AGTG" #4 base sequence on the oligo end of the link between the barcode and oligo - orientation barcode to oligo
  String? end_oligo_link = "CGTC" #4 base sequence at the very end of the oligo

  call Flash { input:
                  read_a=read_a,
                  read_b=read_b,
                  flash_thread=thread,
                  read_len=read_len,
                  id_out=id_out
                }
  call Pull_Barcodes { input:
                          #pull = pull,
                          working_directory=working_directory,
                          merged_fastq=Flash.out,
                          read_number=barcode_orientation,
                          id_out=id_out,
                          barcode_link=barcode_link,
                          seq_min=seq_min,
                          enh_min=enh_min,
                          enh_max=enh_max,
                          bc_len=bc_len,
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
                    map_thread=thread,
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
                  sort_mem=mem,
                  id_out=id_out
                }
  call Ct_Seq { input:
                    #count=count,
                    working_directory=working_directory,
                    sorted=Sort.out,
                    id_out=id_out
                  }
  if (defined(attributes)) {
    call Parse_sat_mut { input:
                  #  parse=parse,
                    attributes=attributes,
                    working_directory=working_directory,
                    counted=Ct_Seq.out,
                    id_out=id_out
                  }
  }
  if (!defined(attributes)) {
    call Parse { input:
                  #  parse=parse,
                    working_directory=working_directory,
                    counted=Ct_Seq.out,
                    id_out=id_out
                  }
  }

  call preseq { input:
                 counted=Ct_Seq.out,
                 id_out=id_out
              }
  call qc_plot_t { input:
                    parsed=Parse.out_parsed,
                    hist=Parse.out_hist,
                    preseq_out=preseq.res,
                    preseq_in=preseq.hist,
                    reference_fasta=reference_fasta,
                    working_directory=working_directory,
                    id_out=id_out
              }
  call relocate { input:
                    flashed=Flash.out,
                    matched=Pull_Barcodes.out1,
                    rejected=Pull_Barcodes.out2,
                    organized_fasta=Rearrange.out,
                    sam_file=MiniMap.out3,
                    map_log=MiniMap.out2,
                    MPRA_out=SAM2MPRA.out,
                    sorted=Sort.out,
                    counted=Ct_Seq.out,
                    parsed=Parse.out_parsed,
                    plothist=Parse.out_hist,
                    preseq_hist=preseq.hist,
                    preseq_res=preseq.res,
                    qc_plot=qc_plot_t.plots,
                    out_directory=out_directory
                  }
  }

task Flash {
  # Flashing raw fastq files together
  File read_a
  File read_b
  Int read_len
  Int flash_thread
  String id_out
  command {
    flash2 -r ${read_len} -f 274 -s 25 -o ${id_out}.merged -t ${flash_thread} ${read_a} ${read_b}
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
  Int enh_min
  Int enh_max
  Int bc_len
  command {
    perl ${working_directory}/pull_barcodes.pl ${merged_fastq} ${read_number} ${id_out}.merged ${barcode_link} ${oligo_link} ${end_oligo_link} ${seq_min} ${enh_min} ${enh_max} ${bc_len}
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
    gzip ${id_out}.merged.match.enh.fa
    >>>
  output {
    File out="${id_out}.merged.match.enh.fa.gz"
    }
  }
task MiniMap {
  # Map the oligos to the reference to get the oligo names
  File reference_fasta
  File organized_fasta
  Int map_thread
  String id_out
  command {
    minimap2 --for-only -Y --secondary=no -m 10 -n 1 -t ${map_thread} --end-bonus 12 -O 5 -E 1 -k 10 -2K50m --eqx --cs=short -c -a ${reference_fasta} ${organized_fasta} > ${id_out}.merged.match.enh.sam 2> ${id_out}.merged.match.enh.log
    samtools view -S -b ${id_out}.merged.match.enh.sam > ${id_out}.merged.match.enh.bam
    }
  output {
    File out1="${id_out}.merged.match.enh.sam"
    File out2="${id_out}.merged.match.enh.log"
    File out3="${id_out}.merged.match.enh.bam"
    }
  }
task SAM2MPRA {
  # Convert the output of the SAM file into a format closer to the dictionary. Pulls the CIGAR and error information
  File sam_file
  String working_directory
  String id_out
  command {
    perl ${working_directory}/SAM2MPRA_cs.pl -C ${sam_file} ${id_out}.merged.match.enh.mapped
    }
  output {
    File out="${id_out}.merged.match.enh.mapped"
    }
  }
task Sort {
  File MPRA_out
  Int sort_mem
  String id_out
  command {
    sort -S${sort_mem}G -k2 ${MPRA_out} > ${id_out}.merged.match.enh.mapped.barcode.sort
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
  command <<<
    perl ${working_directory}/parse_map.pl ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed

    awk '($5 == 0)' ${counted} | awk '{ct[$2]++;cov[$2]+=$4;} END {for(i in ct) print i"\t"ct[i]"\t"cov[i]}' > ${id_out}.merged.match.enh.mapped.barcode.ct.plothist
    >>>
  output {
    File out_parsed="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    File out_hist="${id_out}.merged.match.enh.mapped.barcode.ct.plothist"
    }
  }

task Parse_sat_mut {
  # Parses the barcode oligo pairs to resolve multimapping barcodes
  File? attributes
  File counted
  String working_directory
  String id_out
  Boolean sat_mut=defined(attributes)
  command <<<
    perl ${working_directory}/parse_map.pl -S -A ${attributes} ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed
    awk '($5 == 0)' ${counted} | awk '{ct[$2]++;cov[$2]+=$4;} END {for(i in ct) print i"\t"ct[i]"\t"cov[i]}' > ${id_out}.merged.match.enh.mapped.barcode.ct.plothist
    >>>
  output {
    File out_parsed="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    File out_hist="${id_out}.merged.match.enh.mapped.barcode.ct.plothist"
    }
  }
task qc_plot_t {
  File parsed
  File hist
  File preseq_out
  File preseq_in
  File reference_fasta
  String working_directory
  String id_out
  command {
    Rscript ${working_directory}/mapping_QC_plots.R ${parsed} ${hist} ${preseq_out} ${preseq_in} ${reference_fasta} ${id_out}
    }
  output {
    File plots="${id_out}_barcode_qc.pdf"
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
task relocate{
 File flashed
 File matched
 File rejected
 File organized_fasta
 File sam_file
 File map_log
 File MPRA_out
 File sorted
 File counted
 File parsed
 File plothist
 File preseq_hist
 File preseq_res
 File qc_plot
 String out_directory
 command {
     mv ${flashed} ${matched} ${rejected} ${organized_fasta} ${sam_file} ${map_log} ${MPRA_out} ${sorted} ${counted} ${parsed} ${plothist} ${preseq_hist} ${preseq_res} ${qc_plot} ${out_directory}
   }
 }
