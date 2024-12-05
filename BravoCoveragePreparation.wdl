version 1.0
## Copyright CMG@KIGM, Peter Juvan

# Subworkflows
import "./vcfTasks.wdl" as vcfTasks
import "./coveragePreparation.wdl" as coveragePreparation

workflow BravoCoveragePreparation {
  input {
    # Multi-sample VCF that we use to get sample names in order to generate cram names
    File input_vcf
    File input_vcf_index
    # Cram/crai file locations for CRAM generation step
    String sampleLocationPath
    # Proportion of cram/crai files for coverage preparation
    Float proportion

    # Alternative inputs: if defined, we use these for coverage preparation
    Array[File] inputCramFiles
    Array[File] inputCraiFiles

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    # Reference FASTA file - hg37/hg38
    File referenceFasta
    # Get reference fasta cache using: wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz
    File referenceFastaCache
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  if (length(select_all(inputCramFiles)) == 0) {
    call vcfTasks.VCFquerySamples as VCFqueryCrams {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        prefix = sampleLocationPath,
        suffix = ".cram",
        output_name = vcf_basename + "_smplNames"
    }

    call vcfTasks.read_lines_proportion as subsetCramFiles {
      input:
        in = VCFqueryCrams.out,
        proportion = proportion
    }

    call vcfTasks.VCFquerySamples as VCFqueryCrais {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        prefix = sampleLocationPath,
        suffix = ".cram.crai",
        output_name = vcf_basename + "_smplNames"
    }

    call vcfTasks.read_lines_proportion as subsetCraiFiles {
      input:
        in = VCFqueryCrais.out,
        proportion = proportion
    }

  }

  scatter (chromosome in chromosomes) {
    call coveragePreparation.prepareCoverage {
      input: 
        inputCramFiles = select_first([inputCramFiles,subsetCramFiles.out]),
        inputCraiFiles = select_first([inputCraiFiles,subsetCraiFiles.out]),
        chromosome = chromosome, 
        referenceFasta = referenceFasta,
        referenceFastaCache = referenceFastaCache
    }
  }
    
  output {
    Array[Array[File]] extractDepthFiles = prepareCoverage.extractDepthFiles
    Array[Array[File]] extractDepthIndices = prepareCoverage.extractDepthIndices
    Array[File] aggrBasePair_output = prepareCoverage.aggrBasePair_output
    Array[File] aggrBasePair_outPruneCov0_25 = prepareCoverage.aggrBasePair_outPruneCov0_25
    Array[File] aggrBasePair_outPruneCov0_50 = prepareCoverage.aggrBasePair_outPruneCov0_50
    Array[File] aggrBasePair_outPruneCov0_75 = prepareCoverage.aggrBasePair_outPruneCov0_75
    Array[File] aggrBasePair_outPruneCov1_00 = prepareCoverage.aggrBasePair_outPruneCov1_00
    Array[File] aggrBasePair_output_index = prepareCoverage.aggrBasePair_output_index
    Array[File] aggrBasePair_outPruneCov0_25_index = prepareCoverage.aggrBasePair_outPruneCov0_25_index
    Array[File] aggrBasePair_outPruneCov0_50_index = prepareCoverage.aggrBasePair_outPruneCov0_50_index
    Array[File] aggrBasePair_outPruneCov0_75_index = prepareCoverage.aggrBasePair_outPruneCov0_75_index
    Array[File] aggrBasePair_outPruneCov1_00_index = prepareCoverage.aggrBasePair_outPruneCov1_00_index        
  }
}

