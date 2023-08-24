version 1.0
## Copyright CMG@KIGM, Ales Maver

# Subworkflows
import "./vcfPercentilesPreparation.wdl" as vcfPercentilesPreparation
import "./cramPreparation.wdl" as cramPreparation

# WORKFLOW DEFINITION 
workflow BravoDataPreparation {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 1000000

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    # File for samples
    File samplesFile
    # List of sampls and cram file locations for CRAM generation step
    String sampleLocationPath
    #File sampleLocationFile
    
    # Generate CRAMs optionally (if only update of the frequencies is needed)
    Boolean generate_crams = true

    # Reference FASTA file - hg37/38
    File referenceFasta

    # CAD score files and associated index files
    File cadScores
    File cadScoresIndex

    ### Prepare percentiles ###
    Array[String] infoFields
    Int threads = 40
    Int numberPercentiles = 10
    String description = "Description"

    # Remove reported variants
    File reported_variants

    String vcf_basename = basename(input_vcf, ".vcf.gz")
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
  	call VCFsplitter {
  		input:
  			input_vcf = input_vcf,
  			input_vcf_index = input_vcf_index,
        samplesFile = samplesFile,
  			chromosome = chromosome,
        referenceFasta = referenceFasta
  	}

    call RemoveReportedVariants {
      input:
        input_vcf = VCFsplitter.output_vcf,
        reported_variants = reported_variants
    }

    call VCFindex {
      input:
        input_vcf = RemoveReportedVariants.output_vcf,
        chromosome = chromosome,
        threads = threads
    }

    call VCFfilter {
  		input:
  			input_vcf = VCFindex.output_vcf,
        input_vcf_index = VCFindex.output_vcf_index,
        referenceFasta = referenceFasta,
        threads = threads
  	}

  	call vcfPercentilesPreparation.prepareVCFPercentiles as prepareVCFs {
  		input:
  			input_vcf = VCFfilter.output_vcf,
  			input_vcf_index = VCFfilter.output_vcf_index,
  			samplesFile = samplesFile,
  			#referenceDir = referenceDir, 
  			#lofteeDir = lofteeDir,
  			referenceFasta = referenceFasta,
  			cadScores = cadScores,
  			cadScoresIndex = cadScoresIndex,
  			infoFields = infoFields
		}

    if (generate_crams) {
      call cramPreparation.prepareCram as prepareCRAMs {
        input:
          chromosome = chromosome,
          chromosomeVCF = VCFsplitter.output_vcf,
          chromosomeVCFIndex = VCFsplitter.output_vcf_index,
          samplesFile = samplesFile,
          referenceFasta = referenceFasta,
          sampleLocationPath = sampleLocationPath
          #sampleLocationFile = sampleLocationFile
      }

    }
  } # Close per chromosome scatter

  # Concatenate VCFs from prepare percentiles task
  call concatVcf as concatVcf_RemoveReportedVariants{
    input:
      input_vcfs = VCFindex.output_vcf,
      input_vcfs_indices = VCFindex.output_vcf_index,
      output_name = "output_RemoveReportedVariants",
      threads = threads
  }

  call concatVcf as concatVcf_RemoveReportedVariants_filtered{
    input:
      input_vcfs = VCFfilter.output_vcf,
      input_vcfs_indices = VCFfilter.output_vcf_index,
      output_name = "output_RemoveReportedVariants_filtered",
      threads = threads
  }

  # Concatenate VCFs with removed reported variants
  call concatVcf {
    input:
      input_vcfs = prepareVCFs.output_annotated_vcf,
      input_vcfs_indices = prepareVCFs.output_annotated_vcf_index,
      output_name = "output",
      threads = threads
  }

  if (generate_crams) {
    scatter (chromosome in chromosomes) {
      call concatCrams {
        input:
          input_crams = select_all(prepareCRAMs.combined_cram_result),
          input_cram_indices = select_all(prepareCRAMs.combined_cram_result_index),

          chromosome = chromosome,
          referenceFasta = referenceFasta
      }
    }
  }

  scatter (field in infoFields) {
      call vcfPercentilesPreparation.computePercentiles as computePercentiles {
          input: chromosomeVCF = concatVcf.output_vcf,
              infoField = field,
              threads = threads,
              numberPercentiles = numberPercentiles,
              description = description
      }
  }

  call vcfPercentilesPreparation.addPercentiles as addPercentiles {
    input: 
      chromosomeVCF = concatVcf.output_vcf,
      chromosomeVCFIndex = concatVcf.output_vcf_index,
      variantPercentiles = computePercentiles.outVariantPercentile,
      variantPercentilesIndex = computePercentiles.outVariantPercentileIndex,
      metricJSONs = computePercentiles.outAllPercentiles,
      vcf_basename = "all"
  }

  output {
    File output_vcf = addPercentiles.out
    File output_vcfs_indices = addPercentiles.out_index
    File output_metrics_json = addPercentiles.metrics_json
    #Array[Array[File]] out_metrics_files = prepareVCFs.out_metrics
    Array[File] out_metrics_file = computePercentiles.outAllPercentiles
    File RemoveReportedVariants_output_vcf = concatVcf_RemoveReportedVariants.output_vcf
    File RemoveReportedVariants_output_vcf_index = concatVcf_RemoveReportedVariants.output_vcf_index
    File RemoveReportedVariants_filtered_output_vcf = concatVcf_RemoveReportedVariants_filtered.output_vcf
    File RemoveReportedVariants_filtered_output_vcf_index = concatVcf_RemoveReportedVariants_filtered.output_vcf_index
    Array[File]? out_crams = concatCrams.output_cram
    Array[File]? out_crais = concatCrams.output_cram_index
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
    #cpu: 1
    #memory: "8G"
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
    #cpu: 1
    #memory: "8G"
  }

  # Specify the output declaration to capture the output BED file
  output {
    File converted_regions = "regions.txt"
    Array[String] scatter_regions = read_lines("regions.txt")
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

##############################
task VCFsplitter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File samplesFile

    String chromosome
    File referenceFasta
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    bcftools view -r ~{chromosome} -t ~{chromosome} -S ~{samplesFile} ~{input_vcf} | bcftools norm -m-any -f ~{referenceFasta} -Oz -o ~{chromosome_filename}.~{vcf_basename}.vcf.gz
    bcftools index -t ~{chromosome_filename}.~{vcf_basename}.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    #requested_memory_mb_per_core: 2000
    #cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{chromosome_filename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.~{vcf_basename}.vcf.gz.tbi"
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
    #requested_memory_mb_per_core: 2000
    #cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
  }
}

##############################
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
task VCFclean {
  input {
    File input_vcf
    File input_vcf_index
    String chromosome = "chromosome"
  }

  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command {
    set -e
    bcftools view ~{input_vcf} | bcftools filter -i "QUAL>100" | bcftools +setGT -- -t q -n . -i 'FORMAT/GQ<20' | bcftools view -i 'F_MISSING<1' | bcftools +fill-tags | bcftools filter -e 'INFO/AC=0' --threads 10 -Oz -o ~{chromosome_filename}.clean.vcf.gz
    bcftools index -t ~{chromosome_filename}.clean.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: 10
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "~{chromosome_filename}.clean.vcf.gz"
    File output_vcf_index = "~{chromosome_filename}.clean.vcf.gz.tbi"
  }
}

##############################
task VCFfilter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    File referenceFasta
    Int threads
  }

  command {
    set -e
    #zcat ~{input_vcf} | bcftools view -Oz -o input.vcf.gz
    #bcftools index input.vcf.gz
    bcftools view ~{input_vcf} | bcftools +setGT -- -t q -n . -i 'FORMAT/GQ<20' | bcftools annotate -x FORMAT/PGT,FORMAT/PID | bcftools view --types snps,indels | bcftools view -i 'F_MISSING<1' | bcftools +fill-tags | bcftools filter -e 'INFO/AC=0' | bcftools filter --threads ~{threads} -i "QUAL>100" -Oz -o filtered.vcf.gz
    bcftools index -t filtered.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    requested_memory_mb_per_core: 1000
    cpu: threads
    #runtime_minutes: 180
  }
  output {
    File output_vcf = "filtered.vcf.gz"
    File output_vcf_index = "filtered.vcf.gz.tbi"
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
    # bcftools concat -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}_unsorted.vcf.gz
    # bcftools sort ~{output_name}_unsorted.vcf.gz -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m 9000000000
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