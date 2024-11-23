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
    cpu: 4
    memory: "8G"
    runtime_minutes: 10
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
    cpu: 4
    memory: "8G"
    runtime_minutes: 10
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
  }
}

##############################
## bcftools view -r -t -S | norm -m-any -f ~{referenceFasta}
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
    bcftools index -t ~{chromosome_filename}.~{vcf_basename}.vcf.gz --threads ~{threads}
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{chromosome_filename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.~{vcf_basename}.vcf.gz.tbi"
  }
}

##############################
task VCFindex {
  input {
    # Command parameters
    File input_vcf
    Int threads
  }

  command {
    bcftools index -t ~{input_vcf} --threads ~{threads}
  }

  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    runtime_minutes: 360
  }

  output {
    File output_vcf = input_vcf
    File output_vcf_index = input_vcf + ".tbi"
  }
}

##############################
## bcftools view -r -t -S --force-samples
task VCFsplitSubset {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File? samplesFile
    String region
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String region_filename = sub(sub(region, "-", "_"), ":", "__")

  command {
    set -e
    #bcftools index -t ~{input_vcf} --threads ~{threads}
    bcftools view -r ~{region} -t ~{region} ~{"--force-samples -S " + samplesFile} ~{input_vcf} --threads ~{threads} -Oz -o ~{region_filename}.~{vcf_basename}.vcf.gz
    bcftools index -t ~{region_filename}.~{vcf_basename}.vcf.gz --threads ~{threads}
  }

  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    runtime_minutes: 360
  }

  output {
    File output_vcf = "~{region_filename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{region_filename}.~{vcf_basename}.vcf.gz.tbi"
  }
}

##############################
## bcftools 
##  +setGT -- -t q -n . -i 'FORMAT/GQ<20'     # phred-scaled probability that the call is incorrect
##  annotate -x FORMAT/PGT,FORMAT/PID         # physical phasing haplotype information + physical phasing ID information
##  view --types snps,indels
##  +fill-tags
##  view -i 'F_MISSING<1'                     # Fraction of missing genotypes: include sites with with at least one genotypeany, i.e. not all missing
##  filter -e 'INFO/AC=0'                     # allele count in genotypes, for each ALT (alternative) allele, in the same order as listed: exclude sites with no alelles
##  view -i "QUAL>100"                        # phred-scaled probability that the site has no variant
task VCFfilter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    Int threads
    Float F_MISSING_upper_bounds = 1
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  command {
    set -e
    bcftools view ~{input_vcf} | \
      bcftools +setGT -- -t q -n . -i 'FORMAT/GQ<20' | \
      bcftools annotate -x FORMAT/PGT,FORMAT/PID | \
      bcftools view --types snps,indels | \
      bcftools +fill-tags | \
      bcftools view -i 'F_MISSING<~{F_MISSING_upper_bounds}' | \
      bcftools filter -e 'INFO/AC=0' | \
      bcftools filter --threads ~{threads} -i "QUAL>100" -Oz -o ~{vcf_basename}_flt~{F_MISSING_upper_bounds}.vcf.gz
    bcftools index -t ~{vcf_basename}_flt~{F_MISSING_upper_bounds}.vcf.gz --threads ~{threads}
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 120
  }
  output {
    File output_vcf = "~{vcf_basename}_flt~{F_MISSING_upper_bounds}.vcf.gz"
    File output_vcf_index = "~{vcf_basename}_flt~{F_MISSING_upper_bounds}.vcf.gz.tbi"
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
    bcftools index -t ~{vcf_basename}_norm.vcf.gz --threads ~{threads}
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    runtime_minutes: 10
  }
  output {
    File output_vcf = "~{vcf_basename}_norm.vcf.gz"
    File output_vcf_index = "~{vcf_basename}_norm.vcf.gz.tbi"
  }
}

##############################
## bcftools merge -Oz vcf1 vcf2 ... > vcf_merged
## --force-samples: if the merged files contain duplicate samples names, duplicate sample names will be resolved by prepending the index of the file as it appeared on the command line to the conflicting sample name.
task VCFmerge {
  input {
    # Command parameters
    Array [File] input_vcfs
    Array [File] input_vcfs_indices
    String output_name
    Int threads
  }

  command <<<
    set -e
    bcftools merge --threads ~{threads} --force-samples -Oz -l ~{write_lines(input_vcfs)} > ~{output_name}.vcf.gz
    bcftools index -t ~{output_name}.vcf.gz --threads ~{threads}
  >>>

  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    runtime_minutes: 10
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
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
    wget https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/kigm-dev/removeReportedVariants.py
    python removeReportedVariants.py -i ~{input_vcf} -o ~{output_vcf_filename} -v ~{reported_variants}
  }
  runtime {
    docker: "amancevice/pandas"
    requested_memory_mb_per_core: 2000
    cpu: 8
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
  }
}

##############################
## Called after RemoveReportedVariants
## Applies bcftools +fill-tags to fix AN and AC after removal of reported variants
task VCFfillTags {
  input {
    File input_vcf
    String chromosome = "chromosome"
    Int threads
  }

  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    zcat ~{input_vcf} | bcftools +fill-tags | bcftools view --threads ~{threads} -Oz -o ~{chromosome_filename}.indexed.vcf.gz
    bcftools index  -t ~{chromosome_filename}.indexed.vcf.gz --threads ~{threads}
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 60
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
    bcftools concat --threads ~{threads} -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}_unsorted.vcf.gz
    bcftools index -t ~{output_name}_unsorted.vcf.gz --threads ~{threads}
  >>>

  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: >11h
  }
  output {
    File output_vcf = "~{output_name}_concat.vcf.gz"
    File output_vcf_index = "~{output_name}_concat.vcf.gz.tbi"
  }
}

##############################
## bcftools norm can affect the order of variants in a VCF file; thus we need to sort
## mem is scaled for largemem partition @ Vega (8G per core)
task sortVcf {
    input {
      File input_vcf
      File input_vcf_index
      String output_name
      Int threads
      Float max_mem_scale_factor = 7.5 
    }
  
  command <<<
    set -e
    mkdir $PWD/sort_tmp
    bcftools sort ~{input_vcf} -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m "~{max_mem_scale_factor * threads}G"
    #bcftools sort ~{input_vcf} -Oz -o ~{output_name}.vcf.gz -m "~{max_mem_scale_factor * threads}G"
    bcftools index -t ~{output_name}.vcf.gz --threads ~{threads}
  >>>

  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 8000
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
