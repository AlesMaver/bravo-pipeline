version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow PopulationStatistics {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 3000000

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  
    # Reference FASTA file - hg37/38
    File referenceFasta
  }

  call ConvertIntervalListToBed {
    input:
      interval_list = interval_list
  }

  call SplitRegions {
    input:
      input_bed = ConvertIntervalListToBed.converted_bed,
      thinning_parameter = thinning_parameter,
      scatter_region_size = scatter_region_size
  }

  scatter (chromosome in SplitRegions.scatter_regions ) {
  	call GenerateTable {
  		input:
  			input_vcf = input_vcf,
  			input_vcf_index = input_vcf_index,
  			chromosome = chromosome,
        referenceFasta = referenceFasta
  	    }
    }

  call ConcatenateTabFiles {
    input:
      input_files = GenerateTable.output_tab
  }

  output {
    File output_tab = ConcatenateTabFiles.output_tab
  }

} # Close workflow



# Tasks
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
    cpu: 1
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

    bedtools makewindows -b ~{input_bed} -w $window -s $step |awk '{print $1":"$2"-"$3}' |awk 'NR % ~{default="1" thinning_parameter} == 0' > regions.txt # FOR TESTING - This will subset every n-th row in the regions
    # bedtools makewindows -b ~{input_bed} -w 3000000 |awk '{print $1":"$2"-"$3}' > regions.txt
  >>>

  # Specify the runtime parameters for the task
  runtime {
    docker: "pegi3s/bedtools"  # Use the appropriate Picard Docker image
    cpu: 1
    memory: "8G"
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
  }
}

##############################
task GenerateTable {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index

    String chromosome
    File referenceFasta
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