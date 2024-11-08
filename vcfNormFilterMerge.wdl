## Copyright CMG@KIGM, Peter Juvan
##
## Uses bcftools with multiple VCF files to:
## - subset samples (skipping non-existing), 
## - normalize VCF (left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into biallelic -m-any),
## - filter/annotate (+setGT ./. GQ<20, annotate PGT & PID, --types snps,indels, +fill-tags, include F_MISSING<..., exclude AC=0, include QUAL>100)
## - merge resulting VCFs

version 1.0

# Subworkflows
import "./vcfTasks.wdl" as vcfTasks
import "./VEP.wdl" as VEP

#struct VcfAndIndex {
#  File vcf
#  File vcf_index
#}

workflow vcfNormFilterMerge {
  input {
 #   Array [VcfAndIndex] input_vcfAndInds
    Array [File] input_vcfs

    File interval_list
    Int? thinning_parameter
    Int scatter_region_size = 1000000

    # File for samples
    File samplesFile

    # Reference FASTA file - hg37/38
    File referenceFasta
    Int threads = 5

    # Filter
    Float F_MISSING_upper_bounds = 1

    # Output
    String output_vcf_basename
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

  scatter (region in SplitRegions.scatter_regions) {

    scatter (input_vcf in input_vcfs) {

      call vcfTasks.VCFsplitSubset {
        input:
          input_vcf = input_vcf,
          samplesFile = samplesFile,
          region = region,
          threads = threads
      }

      call vcfTasks.VCFnorm {
        input:
          input_vcf = VCFsplitSubset.output_vcf,
          input_vcf_index = VCFsplitSubset.output_vcf_index,
          referenceFasta = referenceFasta,
          threads = threads
      }

      call vcfTasks.VCFfilter {
        input:
          input_vcf = VCFnorm.output_vcf,
          input_vcf_index = VCFnorm.output_vcf_index,
          threads = threads,
          F_MISSING_upper_bounds = F_MISSING_upper_bounds
      }

    } # Close per input vcf scatter

    call vcfTasks.VCFmerge {
      input:
        input_vcfs = VCFfilter.output_vcf,
        input_vcfs_indices = VCFfilter.output_vcf_index,
        output_name = sub(sub(region, "-", "_"), ":", "__") + "." + output_vcf_basename,
        threads = threads
    }

  } # Close per region scatter

  call vcfTasks.concatVcf {
    input:
      input_vcfs = VCFmerge.output_vcf,
      input_vcfs_indices = VCFmerge.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

  call vcfTasks.sortVcf {
    input:
      input_vcf = concatVcf.output_vcf,
      input_vcf_index = concatVcf.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

  output {
    File output_vcf = sortVcf.output_vcf
    File output_vcf_index = sortVcf.output_vcf_index
  }

} # Close workflow

