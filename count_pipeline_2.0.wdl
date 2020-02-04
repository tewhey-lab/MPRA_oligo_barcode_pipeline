###Wdl version of pipeline for organizing MPRA count data
workflow MPRACount {
  File read_a
  File read_b
  File pull
  File reference_fasta
  File sam_convert
  File count
  File parse
  File proj_list
  Int read_number
  Int seq_min
  String id_out
  String barcode_link
  String oligo_link
  String end_oligo_link

  call Flash { input:
                  read_a=read_a,
                  read_b=read_b,
                  id_out=id_out
                }
  call Pull_Barcodes { input:
                          pull=pull,
                          merged_fastq=Flash.out,
                          read_number=read_number,
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
                      sam_convert=sam_convert,
                      sam_file=MiniMap.out1,
                      id_out=id_out
                    }
  call Sort { input:
                  MPRA_out=SAM2MPRA.out,
                  id_out=id_out
                }
  call Ct_Seq { input:
                    count=count,
                    sorted=Sort.out,
                    id_out=id_out
                  }
  call Parse { input:
                  counted=Ct_Seq.out,
                  parse=parse,
                  id_out=id_out
                }
  call make_attr_list { input:
                          proj_list=proj_list,
                          reference_fasta=reference_fasta,
                          id_out=id_out
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
  File pull
  File merged_fastq
  Int read_number
  String id_out
  String barcode_link
  String oligo_link
  String end_oligo_link
  Int seq_min
  command {
    perl ${pull} ${merged_fastq} ${read_number} ${id_out}.merged ${barcode_link} ${oligo_link} ${end_oligo_link} ${seq_min}
  }
  output {
    File out1="${id_out}.merged.match"
    File out2="${id_out}.merged.reject"
  }
}

task Rearrange {
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
  File sam_convert
  File sam_file
  String id_out
  command {
    perl ${sam_convert} -C ${sam_file} ${id_out}.merged.match.enh.mapped
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
  File count
  File sorted
  String id_out
  command {
    perl ${count} ${sorted} 2 4 > ${id_out}.merged.match.enh.mapped.barcode.ct
  }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct"
  }
}

task Parse {
  File counted
  File parse
  String id_out
  command {
    perl ${parse} ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed
  }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
  }
}

task make_attr_list {
  File proj_list
  File reference_fasta
  String id_out
  command {
    perl ${proj_list} ${reference_fasta} ${id_out}
  }
  output {
    File out="${id_out}.proj_list"
  }
}
