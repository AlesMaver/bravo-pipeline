version 1.0
## Copyright CMG@KIGM, Peter Juvan

# Subworkflows
import "./vcfTasks.wdl" as vcfTasks

workflow vcfMerge {
  input {
    Array [File] input_vcfs

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 1000000

    Int threads = 5
    String output_vcf_basename = basename(select_first(input_vcfs), ".vcf.gz") + "_merged"
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

    scatter (input_vcf in input_vcfs) {
      
      call vcfTasks.VCFsplit {
        input:
          input_vcf = input_vcf,
          input_vcf_index = basename(input_vcf, ".vcf.gz") + ".vcf.gz.tbi",
          chromosome = chromosome,
          threads = threads
      }

    } # close vcf scatter

    call vcfTasks.VCFmerge {
      input:
        input_vcfs = VCFsplit.output_vcf,
        output_name = output_vcf_basename
        threads = threads
    }

  } # Close per chromosome scatter

  call vcfTasks.concatVcf {
    input:
      input_vcfs = VCFmerge.output_vcf,
      input_vcfs_indices = VCFmerge.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

  output {
    File output_vcf = concatVcf.output_vcf
    File output_vcf_index = concatVcf.output_vcf_index
  }

} # Close workflow

