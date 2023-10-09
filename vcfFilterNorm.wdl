version 1.0
## Copyright CMG@KIGM, Peter Juvan

# Subworkflows
import "./vcfTasks.wdl" as vcfTasks

workflow vcfFilterNorm {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 1000000

    # File for samples
    File samplesFile

    # Reference FASTA file - hg37/38
    File referenceFasta
    Int threads = 5
    String output_vcf_basename = basename(input_vcf, ".vcf.gz") + "_nrmFlt"
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

    call vcfTasks.VCFsplit {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        samplesFile = samplesFile,
        chromosome = chromosome,
        threads = threads
    }

    call vcfTasks.VCFnorm {
      input:
        input_vcf = VCFsplit.output_vcf,
        input_vcf_index = VCFsplit.output_vcf_index,
        referenceFasta = referenceFasta,
        threads = threads
    }

    call vcfTasks.VCFfilter {
      input:
        input_vcf = VCFnorm.output_vcf,
        input_vcf_index = VCFnorm.output_vcf_index,
        threads = threads
    }

  } # Close per chromosome scatter

  call vcfTasks.concatVcf {
    input:
      input_vcfs = VCFfilter.output_vcf,
      input_vcfs_indices = VCFfilter.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

  output {
    File output_vcf = concatVcf.output_vcf
    File output_vcf_index = concatVcf.output_vcf_index
  }

} # Close workflow

