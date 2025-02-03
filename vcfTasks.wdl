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
    requested_memory_mb_per_core: 2000
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
    requested_memory_mb_per_core: 2000
    runtime_minutes: 10
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
  }
}

##############################
## DEPRECATED for VCFsplitSubset + VCFnorm
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
    bcftools view --threads ~{threads} -r ~{chromosome} -t ~{chromosome} -S ~{samplesFile} ~{input_vcf} | \
      bcftools norm -m-any -f ~{referenceFasta} --threads ~{threads} -Oz -o ~{chromosome_filename}.~{vcf_basename}.vcf.gz --write-index=tbi
  }
  runtime {
    docker: "peterjuv/bcftools"
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
## Alternative to VCFindex
## For cases like the following: 
##  Contig 'chr3' is not defined in the header. (Quick workaround: index the file with tabix.)
##  Undefined tags in the header, cannot proceed in the sample subset mode.
task VCFtabix {
  input {
    # Command parameters
    File input_vcf
    Int threads # not used, we keep it here to be consistent with inputs of task VCFindex
  }

  String vcf_basename = basename(input_vcf)

  command {
    tabix --force --preset vcf ~{input_vcf}
    ln ~{input_vcf}.tbi ~{vcf_basename}.tbi
  }

  runtime {
    docker: "ensemblorg/ensembl-vep:latest"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 360
  }

  output {
    File output_vcf = input_vcf
    File output_vcf_index = "~{vcf_basename}.tbi"
  }
}

##############################
task VCFindex {
  input {
    # Command parameters
    File input_vcf
    Int threads
  }

  String vcf_basename = basename(input_vcf)

  command {
    bcftools index -t ~{input_vcf} --threads ~{threads} -o $PWD/~{vcf_basename}.tbi
  }

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 360
  }

  output {
    File output_vcf = input_vcf
    File output_vcf_index = "~{vcf_basename}.tbi"
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
    bcftools view -r ~{region} -t ~{region} ~{"--force-samples -S " + samplesFile} ~{input_vcf} --threads ~{threads} -Oz -o ~{region_filename}.~{vcf_basename}.vcf.gz --write-index=tbi
  }

  runtime {
    docker: "peterjuv/bcftools"
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
    bcftools view --threads ~{threads} ~{input_vcf} | \
      bcftools +setGT --threads ~{threads} -- -t q -n . -i 'FORMAT/GQ<20' | \
      bcftools annotate --threads ~{threads} -x FORMAT/PGT,FORMAT/PID | \
      bcftools view --threads ~{threads} --types snps,indels | \
      bcftools +fill-tags --threads ~{threads} | \
      bcftools view --threads ~{threads} -i 'F_MISSING<~{F_MISSING_upper_bounds}' | \
      bcftools filter --threads ~{threads} -e 'INFO/AC=0' | \
      bcftools filter --threads ~{threads} -i "QUAL>100" -Oz -o ~{vcf_basename}_flt~{F_MISSING_upper_bounds}.vcf.gz --write-index=tbi
  }
  runtime {
    docker: "peterjuv/bcftools"
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
    bcftools view ~{input_vcf} | bcftools norm -m-any -f ~{referenceFasta} --threads ~{threads} -Oz -o ~{vcf_basename}_norm.vcf.gz --write-index=tbi
  }
  runtime {
    docker: "peterjuv/bcftools"
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
    bcftools merge --threads ~{threads} --force-samples --force-single -l ~{write_lines(input_vcfs)} -Oz -o ~{output_name}.vcf.gz
    bcftools index -t ~{output_name}.vcf.gz --threads ~{threads}
  >>>

  runtime {
    # bcftools v1.21 causes segmentation fault, thus we use v.1.20 here
    docker: "peterjuv/bcftools:v.1.20"
    requested_memory_mb_per_core: 2000
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
    Int threads = 8
  }

  String output_vcf_filename = basename(input_vcf, ".vcf.gz") + "_remReported.vcf.gz"

  command {
    set -e
    wget https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/kigm-dev/removeReportedVariants.py
    python3 removeReportedVariants.py -i ~{input_vcf} -o ~{output_vcf_filename} -v ~{reported_variants}
    bcftools index -t ~{output_vcf_filename} --threads ~{threads}
  }
  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 120
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

##############################
## Called after RemoveReportedVariants
## Applies bcftools +fill-tags to fix AN and AC after removal of reported variants
task VCFfillTags {
  input {
    File input_vcf
    File input_vcf_index
    Int threads
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  command {
    set -e
    zcat ~{input_vcf} | bcftools +fill-tags --threads ~{threads} | bcftools view --threads ~{threads} -Oz -o ~{vcf_basename}_fillTags.vcf.gz --write-index=tbi
  }

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{vcf_basename}_fillTags.vcf.gz"
    File output_vcf_index = "~{vcf_basename}_fillTags.vcf.gz.tbi"
  }
}

##############################
## bcftools -G to remove individual genotype information
task VCFdropGeno {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String output_name
    Int threads
  }
  
  command <<<
    set -e
    bcftools view ~{input_vcf} -G --threads ~{threads} -Oz -o ~{output_name}.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: ?
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
## Set all IDs to the following format
## bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT'
task VCFsetIDs {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String output_name
    Int threads
  }
  
  command <<<
    set -e
    bcftools annotate ~{input_vcf} --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' --threads ~{threads} -Oz -o ~{output_name}.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: ?
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
## Set all IDs to the following format
## bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT'
task VCFsetMissingIDs {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String output_name
    Int threads
  }
  
  command <<<
    set -e
    bcftools annotate ~{input_vcf} --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' --threads ~{threads} -Oz -o ~{output_name}.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: ?
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
## DEPRECATED for concatIdxVcf
## Note that we do not index here due to time limits of individual tasks
task concatVcf {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indices
    String output_name
    Int threads
  }
  
  command <<<
    set -e
    bcftools concat --threads ~{threads} -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}.vcf.gz
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: >11h
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    #File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
task concatIdxVcf {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indices
    String output_name
    Int threads
  }
  
  command <<<
    set -e
    bcftools concat --threads ~{threads} -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: >11h
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
## bcftools concat --allow-overlaps --remove-duplicates
## Note that records duplicate within one file are not removed 
## Every line that finds a matching record in another file will be printed only once
## Alias --rm-dups exact
task concatOverlapsIdxVcf {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indices
    String output_name
    Int threads
    Int memory_mb_per_core = 2000
  }
  
  command <<<
    set -e
    bcftools concat --allow-overlaps --remove-duplicates --threads ~{threads} -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools:v.1.20"
    requested_memory_mb_per_core: memory_mb_per_core
    cpu: threads
    #runtime_minutes: >11h
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}


##############################
## bcftools norm can affect the order of variants in a VCF file; thus we need to sort
## mem can be scaled for largemem partition @ Vega by setting memory_mb_per_core = 8000 
## we do not use --temp-dir because we want to use /scratch/slurm/$SLURM_JOB_ID @ Vega
## Note that we do not index here due to time limits of individual tasks
## CONSIDER using sortIdxVcf instead
task sortVcf {
  input {
    File input_vcf
    File input_vcf_index
    String output_name
    Int threads
    Int memory_mb_per_core = 2000
  }
  
  command <<<
    set -e
    mkdir $PWD/sort_tmp
    bcftools sort ~{input_vcf} -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m "~{memory_mb_per_core/1000*threads-1}G"
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: memory_mb_per_core
    cpu: threads
    #runtime_minutes: 2880
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
  }
}

##############################
## bcftools norm can affect the order of variants in a VCF file; thus we need to sort
## mem can be scaled for largemem partition @ Vega by setting memory_mb_per_core = 8000 
## we do not use --temp-dir because we want to use /scratch/slurm/$SLURM_JOB_ID @ Vega
task sortIdxVcf {
  input {
    File input_vcf
    File input_vcf_index
    String output_name
    Int threads
    Int memory_mb_per_core = 2000
  }
  
  command <<<
    set -e
    mkdir $PWD/sort_tmp
    bcftools sort ~{input_vcf} -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m "~{memory_mb_per_core/1000*threads-1}G" --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: memory_mb_per_core
    cpu: threads
    #runtime_minutes: 2880
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
## DEPRECATED for concatIdxVcf + sortIdxVcf
## FOR TESTING PURPOSES
## Note that sort is much slower than concat, thus makes sense to run these in separate tasks
task concatSortIdxVcf {
    input {
      Array[File] input_vcfs
      Array[File] input_vcfs_indices
      String output_name
      Int threads
      Int memory_mb_per_core = 2000
    }
  
  command <<<
    set -e
    mkdir $PWD/sort_tmp
    bcftools concat --threads ~{threads} -f ~{write_lines(input_vcfs)} | \
      bcftools sort -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m "~{memory_mb_per_core/1000*threads-1}G" --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: threads
    #runtime_minutes: >11h
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}


##############################
## Get samples names from a VCF, optionally append prefix (path) and suffix (extension), and write to a file
task VCFquerySamples {
  input {
    File input_vcf
    File input_vcf_index
    String? prefix
    String? suffix
    String output_name
  }
  
  command <<<
    set -e
    bcftools query -l ~{input_vcf} | awk '{print prefix$1suffix}' prefix=~{prefix + "/"} suffix=~{suffix} > ~{output_name}.tab
  >>>

  runtime {
    docker: "peterjuv/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 1
    runtime_minutes: 10
  }
  output {
    File out = "~{output_name}.tab"
  }
}

##############################
## Read a non-random proportion of lines from input file
task read_lines_proportion {
  input {
    File in
    Float proportion
  }
  
  command <<<
    set -e
    cat ~{in} | awk 'rand()<proportion' proportion=~{proportion}
  >>>

  runtime {
    docker: "bashell/alpine-bash:latest"
    requested_memory_mb_per_core: 2000
    cpu: 1
    runtime_minutes: 10
  }
  output {
    Array[File] out = read_lines(stdout())
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
    mkdir -p cram
    # Ensure output files are in cram dir to allow continuation in case of an empty cram
    touch cram/~{chromosome}.cram
    touch cram/~{chromosome}.cram.crai
    touch chromosome.cram.list
    touch chromosome.cram.non_empty.list

    cat ~{write_lines(input_crams)} | grep "~{chromosome}__" > chromosome.cram.list
    cat chromosome.cram.list | xargs -I {} find {} -type f -empty -prune -o -print > chromosome.cram.non_empty.list

    if [ -s "chromosome.cram.non_empty.list" ]; then
        echo "At least one input CRAM file found, therefore merging!"
        samtools merge -b chromosome.cram.non_empty.list -O CRAM cram/~{chromosome}.cram --reference ~{referenceFasta} -f
        samtools index cram/~{chromosome}.cram
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
    File output_cram = "cram/~{chromosome}.cram"
    File output_cram_index = "cram/~{chromosome}.cram.crai"
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
