###Wdl version of pipeline for organizing MPRA count data
workflow MPRACount {
  Array[String] replicate_fastq
  Array[String] replicate_id
  Array[Pair[String,String]] fastq_id = zip(replicate_fastq, replicate_id)
  File read_a
  File read_b
  File pull
  File reference_fasta
  File sam_convert
  File count
  File parse
  File list_proj
  File oligo_type
  File make_counts
  File associate_tags
  File compile
  Int read_number
  Int seq_min
  String id_out
  String barcode_link
  String oligo_link
  String end_oligo_link
  String flags

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
                          list_proj=list_proj,
                          reference_fasta=reference_fasta,
                          id_out=id_out
                        }
  call make_attr_file { input:
                          proj_list=make_attr_list.out,
                          oligo_type=oligo_type,
                          id_out=id_out
                        }
  scatter (replicate in fastq_id) {
    call prep_counts { input:
                          make_counts=make_counts,
                          sample_fastq=replicate.left,
                          parsed=Parse.out,
                          sample_id=replicate.right
                        }
    call associate { input:
                        associate_tags=associate_tags,
                        matched=prep_counts.out,
                        parsed=Parse.out,
                        sample_id=replicate.right
                      }
                    }
  call make_infile { input:
                        tag_files=associate.outF,
                        tag_ids=associate.outS,
                        id_out=id_out
                      }
  call make_count_table { input:
                            list_inFile=make_infile.out,
                            compile=compile,
                            flags=flags,
                            id_out=id_out,
                            replicate_id=replicate_id
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
    awk '{print ">"$1"#"$3"\n"$4}' ~{matched_barcodes} > ~{id_out}.merged.match.enh.fa
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
  File list_proj
  File reference_fasta
  String id_out
  command {
    perl ${list_proj} ${reference_fasta} ${id_out}
  }
  output {
    File out="${id_out}.proj_list"
  }
  }
task make_attr_file {
  File proj_list
  File oligo_type
  String id_out
  command {
    perl ${oligo_type} ${proj_list} ${id_out}
  }
  output {
    File out="${id_out}.attributes"
  }
  }
task prep_counts {
  File make_counts
  File sample_fastq
  File parsed
  String sample_id
  command {
    python ${make_counts} ${sample_fastq} ${parsed} ${sample_id}
  }
  output {
    File out="${sample_id}.match"
  }
  }
task associate {
  File associate_tags
  File matched
  File parsed
  String sample_id
  command {
    perl ${associate_tags} ${matched} ${parsed} ${sample_id}.tag
  }
  output {
    File outF="${sample_id}.tag"
    String outS="${sample_id}"
  }
  }
task make_infile {
  Array[File] tag_files
  Array[String] tag_ids
  String id_out
  Array[Pair[String, File]] tags=zip(tag_ids,tag_files)
  command <<<
    for pair in ~{sep=' ' tags}; do
      printf "%s\t%s\n" "~{pair.left}" "~{pair.right}" >> ~{id_out}_samples.txt
    done
  >>>
  output {
    File out="${id_out}_samples.txt"
  }
  }
task make_count_table {
  File list_inFile
  File compile
  Array[String] replicate_id
  String? flags = ""
  String id_out
  command {
    perl ${compile} ${flags} ${list_inFile} ${id_out}.count ${replicate_id}
  }
  output {
    File out="$id_out.count"
  }
  }
