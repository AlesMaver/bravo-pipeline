version 1.0
## Copyright CMG@KIGM, Peter Juvan

# Subworkflows
import "./vcfPercentilesPreparation.wdl" as vcfPercentilesPreparation

# WORKFLOW DEFINITION 
workflow vcfFilterNorm {
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

    # Reference FASTA file - hg37/38
    File referenceFasta
    Int threads = 40
    String vcf_basename = basename(input_vcf, ".vcf.gz")
  }

  call vcfTasks.ConvertIntervalListToBed {
    input:
      interval_list = interval_list
  }

  call vcfTasks.SplitRegions {
    input:
      input_bed = ConvertIntervalListToBed.converted_bed,
      thinning_parameter = thinning_parameter,
      scatter_region_size = scatter_region_size
  }

  scatter (chromosome in SplitRegions.scatter_regions ) {
#    call vcfTasks.VCFsplitter {
#      input:
#        input_vcf = input_vcf,
#        input_vcf_index = input_vcf_index,
#        samplesFile = samplesFile,
#        chromosome = chromosome,
#        referenceFasta = referenceFasta,
#        threads = threads
#    }

    call vcfTasks.VCFsplit {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        samplesFile = samplesFile,
        chromosome = chromosome,
        threads = threads
    }

    call vcfTasks.VCFfilter {
      input:
        input_vcf = VCFsplit.output_vcf,
        input_vcf_index = VCFsplit.output_vcf_index,
        threads = threads
    }

    call vcfTasks.VCFnorm {
      input:
        input_vcf = VCFfilter.output_vcf,
        input_vcf_index = VCFfilter.output_vcf_index,
        referenceFasta = referenceFasta,
        threads = threads
    }


  } # Close per chromosome scatter

  call vcfTasks.concatVcf {
    input:
      input_vcfs = VCFnorm.output_vcf,
      input_vcfs_indices = VCFnorm.output_vcf_index,
      output_name = vcf_basename + "_fltNorm"
      threads = threads
  }

  output {
    File output_vcf = concatVcf.output_vcf
    File output_vcf_index = concatVcf.output_vcf_index
  }
} # Close workflow

