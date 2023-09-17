version 1.0
## Copyright CMG@KIGM, Ales Maver & Peter Juvan

##############################
task ConvertIntervalListToBed {
  input {
    File interval_list
  }

  # Command section where the conversion is performed using Picard
  command <<<
    # Convert an interval list to BED 
    java  -Xmx14g -jar /usr/picard/picard.jar IntervalListToBed \
      I=~{interval_list} \
      O=interval.bed

    # Format the scatter regions
    # awk '{print $1":"$2"-"$3}' interval.bed |awk 'NR % 50 == 0' > regions.txt # FOR TESTING - This will subset every 50th row in the regions
    awk '{print $1":"$2"-"$3}' interval.bed > regions.txt

    #output_regions=$(cat output_regions.txt)
  >>>

  # Specify the runtime parameters for the task
  runtime {
    docker: "broadinstitute/picard:2.26.0"  # Use the appropriate Picard Docker image
    #cpu: 1
    memory: "8G"
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_bed = "interval.bed"
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
  }
}

##############################
task SplitRegions {
  input {
    File input_bed
    Int? thinning_parameter
    Int scatter_region_size
  }

  # Command section where the conversion is performed using Picard
  command <<<
    # Convert an interval list to BED 
    window=~{scatter_region_size}
    step=$(($window + 1))
    #step=$window

    # FOR TESTING - This will subset every n-th row in the regions
    bedtools makewindows -b ~{input_bed} -w $window -s $step |awk '{print $1":"$2"-"$3}' |awk 'NR % ~{default="1" thinning_parameter} == 0' > regions.txt 
    # bedtools makewindows -b ~{input_bed} -w 3000000 |awk '{print $1":"$2"-"$3}' > regions.txt
  >>>

  # Specify the runtime parameters for the task
  runtime {
    docker: "pegi3s/bedtools"  # Use the appropriate Picard Docker image
    #cpu: 1
    memory: "8G"
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
  }
}

##############################
## bcftools -r -t -S | norm -m-any -f ~{referenceFasta}
task VCFsplitter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File samplesFile
    String chromosome
    File referenceFasta
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    bcftools view -r ~{chromosome} -t ~{chromosome} -S ~{samplesFile} ~{input_vcf} | bcftools norm -m-any -f ~{referenceFasta} --threads ~{threads} -Oz -o ~{chromosome_filename}.~{vcf_basename}.vcf.gz
    bcftools index -t ~{chromosome_filename}.~{vcf_basename}.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{chromosome_filename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.~{vcf_basename}.vcf.gz.tbi"
  }
}

##############################
## bcftools -r -t -S
task VCFsplit {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File samplesFile
    String chromosome
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    bcftools view -r ~{chromosome} -t ~{chromosome} -S ~{samplesFile} ~{input_vcf} --threads ~{threads} -Oz -o ~{chromosome_filename}.~{vcf_basename}.vcf.gz
    bcftools index -t ~{chromosome_filename}.~{vcf_basename}.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{chromosome_filename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.~{vcf_basename}.vcf.gz.tbi"
  }
}

##############################
## bcftools +setGT -- -t q -n . -i 'FORMAT/GQ<20' | annotate -x FORMAT/PGT,FORMAT/PID | --types snps,indels | -i 'F_MISSING<1' |
##          +fill-tags | filter -e 'INFO/AC=0' | -i "QUAL>100" 
task VCFfilter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  command {
    set -e
    #zcat ~{input_vcf} | bcftools view -Oz -o input.vcf.gz
    #bcftools index input.vcf.gz
    bcftools view ~{input_vcf} | bcftools +setGT -- -t q -n . -i 'FORMAT/GQ<20' | bcftools annotate -x FORMAT/PGT,FORMAT/PID | bcftools view --types snps,indels | bcftools +fill-tags | bcftools view -i 'F_MISSING<1' | bcftools filter -e 'INFO/AC=0' | bcftools filter --threads ~{threads} -i "QUAL>100" -Oz -o ~{vcf_basename}_flt.vcf.gz
    bcftools index -t ~{vcf_basename}_flt.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{vcf_basename}_flt.vcf.gz"
    File output_vcf_index = "~{vcf_basename}_flt.vcf.gz.tbi"
  }
}

##############################
## bcftools norm -m-any -f ~{referenceFasta}
task VCFnorm {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File referenceFasta
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  command {
    set -e
    bcftools view ~{input_vcf} | bcftools norm -m-any -f ~{referenceFasta} --threads ~{threads} -Oz -o ~{vcf_basename}_norm.vcf.gz
    bcftools index -t ~{vcf_basename}_norm.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{vcf_basename}_norm.vcf.gz"
    File output_vcf_index = "~{vcf_basename}_norm.vcf.gz.tbi"
  }
}


##############################
task RemoveReportedVariants {
  input {
    # Command parameters
    File input_vcf
    File reported_variants
  }

  String output_vcf_filename = "cleaned.vcf.gz"

  command {
    set -e
    wget https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/region_scatter/removeReportedVariants.py
    python removeReportedVariants.py -i ~{input_vcf} -o ~{output_vcf_filename} -v ~{reported_variants}
  }
  runtime {
    docker: "amancevice/pandas"
    requested_memory_mb_per_core: 2000
    cpu: 8
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
  }
}

##############################
## Called after RemoveReportedVariants
## Applies bcftools +fill-tags to fix AN and AC after removal of reported variants
task VCFindex {
  input {
    File input_vcf
    String chromosome = "chromosome"
    Int threads
  }

  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    zcat ~{input_vcf} | bcftools +fill-tags | bcftools view --threads ~{threads} -Oz -o ~{chromosome_filename}.indexed.vcf.gz
    bcftools index  -t ~{chromosome_filename}.indexed.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{chromosome_filename}.indexed.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.indexed.vcf.gz.tbi"
  }
}


##############################
task concatVcf {
    input {
      Array[File] input_vcfs
      Array[File] input_vcfs_indices
      String output_name
      Int threads
    }
  
  command <<<
  set -e
    mkdir $PWD/sort_tmp
    bcftools concat --threads ~{threads} -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}.vcf.gz
    bcftools index -t ~{output_name}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 90
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
task concatCrams {
    input {
      Array[File] input_crams
      Array[File] input_cram_indices

      String chromosome
      File referenceFasta
    }
  
  command <<<
    # Ensure output files are in the executions dir to allow continuation in case of an empty cram
    touch ~{chromosome}.cram
    touch ~{chromosome}.cram.crai
    touch chromosome.cram.list
    touch chromosome.cram.non_empty.list

    cat ~{write_lines(input_crams)} | grep "~{chromosome}__" > chromosome.cram.list
    cat chromosome.cram.list | xargs -I {} find {} -type f -empty -prune -o -print > chromosome.cram.non_empty.list

    if [ -s "chromosome.cram.non_empty.list" ]; then
        echo "At least one input CRAM file found, therefore merging!"
        samtools merge -b chromosome.cram.list -O CRAM ~{chromosome}.cram --reference ~{referenceFasta} -f
        samtools index ~{chromosome}.cram
    else
        echo "No input CRAM files, therefore leaving empty final crams!"
    fi
  >>>

  runtime {
    docker: "alesmaver/bravo-pipeline-sgp:latest"
    #requested_memory_mb_per_core: 5000
    #cpu: 1
    #runtime_minutes: 90
  }
  output {
    File output_cram = "~{chromosome}.cram"
    File output_cram_index = "~{chromosome}.cram.crai"
  }
}

##############################
## Misc
##############################
task GenerateTable {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String chromosome
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    bcftools view -r ~{chromosome} -t ~{chromosome} ~{input_vcf} | \
    bcftools +fill-tags | \
    bcftools +split-vep \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%IMPACT\t%MANE_SELECT\t%CANONICAL\t%EXON\t%HGVSc\t%clinvar_clnsig\t%clinvar_review\t%INFO/AF\t%AC\t%AC_Hom\t%LoF\t%LoF_filter\t%LoF_flags\t%LoF_info\t%gnomADe_NFE_AF\t[%SAMPLE,]\n' \
      -s primary \
      -i'(clinvar_clnsig ~ "pathogenic/i") && (clinvar_clnsig !~ "conflicting/i") && GT="alt"' \
      -d \
    > ~{chromosome_filename}.~{vcf_basename}.tab
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_tab = "~{chromosome_filename}.~{vcf_basename}.tab"
  }
}

##############################
task ConcatenateTabFiles {
  input {
    Array[File] input_files
  }
  command {
    cat $(cat ~{write_lines(input_files)}) > MergedVariantTable.tab
  }
    runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_tab = "MergedVariantTable.tab"
  }
}

##############################
task PrintStringToStdout {
  # Define the input string
  input {
    String input_string
  }

  # Command section where the input string is printed to stdout
  command {
    echo "${input_string}"
  }

  # Specify the output declaration (optional in this case)
  output {
    String output_string = read_string(stdout())
  }
}
