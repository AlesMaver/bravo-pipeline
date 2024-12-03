version 1.0
## Copyright CMG@KIGM, Ales Maver

# Subworkflows
import "./vcfTasks.wdl" as vcfTasks
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

    # File for sample names
    File? samplesFile
    # Cram/crai file locations for CRAM generation step
    String sampleLocationPath
    
    # Generate CRAMs optionally (if only update of the frequencies is needed)
    Boolean generate_crams = true

    # Reference FASTA file - hg37/38
    File referenceFasta

    # CAD score files and associated index files
    File cadScores
    File cadScoresIndex

    ### Prepare percentiles ###
    Array[String] infoFields
    Int threads = 4
    Int numberPercentiles = 10
    String description = "Description"

    # Remove reported variants
    File reported_variants

    # VEP
    VEPReferences vep_ref
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz") 

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

  if (!defined(samplesFile)) {
    call vcfTasks.VCFquerySamples {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        output_name = vcf_basename + "_smplNames"
    }
  }

  scatter (region in SplitRegions.scatter_regions ) {

    call vcfTasks.VCFsplitSubset {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        samplesFile = samplesFile,
        region = region,
        threads = threads
    }

    call vcfTasks.RemoveReportedVariants {
      input:
        input_vcf = VCFsplitSubset.output_vcf,
        reported_variants = reported_variants,
        threads = threads
    }

    call vcfTasks.VCFfillTags {
      input:
        input_vcf = RemoveReportedVariants.output_vcf,
        input_vcf_index = RemoveReportedVariants.output_vcf_index,
        threads = threads
    }

  	call vcfPercentilesPreparation.prepareVCFPercentiles as prepareVCFs {
  		input:
        input_vcf = VCFfillTags.output_vcf,
        input_vcf_index = VCFfillTags.output_vcf_index,
        samplesFile = select_first([samplesFile, VCFquerySamples.out]),
        cadScores = cadScores,
        cadScoresIndex = cadScoresIndex,
        threads = threads,
        numberPercentiles = numberPercentiles,
        description = description,
        vep_ref = vep_ref
    }

    if (generate_crams) {
      call cramPreparation.prepareCram as prepareCRAMs {
        input:
          chromosome = region,
          chromosomeVCF = VCFsplitSubset.output_vcf,
          chromosomeVCFIndex = VCFsplitSubset.output_vcf_index,
          samplesFile = select_first([samplesFile, VCFquerySamples.out]),
          referenceFasta = referenceFasta,
          sampleLocationPath = sampleLocationPath,
          threads = threads
      }

    }
  } # Close per region scatter

  # Concatenate VCFs with removed reported variants
  call vcfTasks.concatIdxVcf as concatIdxVcf_RemoveReportedVariants {
    input:
      input_vcfs = VCFfillTags.output_vcf,
      input_vcfs_indices = VCFfillTags.output_vcf_index,
      output_name = vcf_basename + "_repVarRem",
      threads = threads
  }

  # Concatenate VCFs (annotated with VEP) from prepare percentiles task
  call vcfTasks.concatIdxVcf {
    input:
      input_vcfs = prepareVCFs.output_annotated_vcf,
      input_vcfs_indices = prepareVCFs.output_annotated_vcf_index,
      output_name = vcf_basename + "_repVarRem_annotated",
      threads = threads
  }

  if (generate_crams) {
    scatter (chromosome in chromosomes) {
      call vcfTasks.concatCrams {
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
      input: 
        chromosomeVCF = concatIdxVcf.output_vcf,
        infoField = field,
        threads = threads,
        numberPercentiles = numberPercentiles,
        description = description
    }
  }

  call vcfPercentilesPreparation.addPercentiles as addPercentiles {
    input: 
      chromosomeVCF = concatIdxVcf.output_vcf,
      chromosomeVCFIndex = concatIdxVcf.output_vcf_index,
      variantPercentiles = computePercentiles.outVariantPercentile,
      variantPercentilesIndex = computePercentiles.outVariantPercentileIndex,
      metricJSONs = computePercentiles.outAllPercentiles,
      vcf_basename = "all"
  }

  output {
    File output_vcf = addPercentiles.out
    File output_vcfs_indices = addPercentiles.out_index
    File output_metrics_json = addPercentiles.metrics_json
    Array[File] out_metrics_file = computePercentiles.outAllPercentiles
    File RemoveReportedVariants_output_vcf = concatIdxVcf_RemoveReportedVariants.output_vcf
    File RemoveReportedVariants_output_vcf_index = concatIdxVcf_RemoveReportedVariants.output_vcf_index
    Array[File]? out_crams = concatCrams.output_cram
    Array[File]? out_crais = concatCrams.output_cram_index
  }
} # Close workflow
