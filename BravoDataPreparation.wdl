version 1.0
## Copyright CMG@KIGM, Ales Maver

# Subworkflows
import "https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/vcfPercentilesPreparation.wdl" as vcfPercentilesPreparation
import "https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/cramPreparation.wdl" as cramPreparation

# WORKFLOW DEFINITION 
workflow BravoDataPreparation {
  input {
    File input_vcf
    File input_vcf_index

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

    # File for samples
    File samplesFile
    # List of sampls and cram file locations for CRAM generation step
    File sampleLocationFile

    # Directory for reference data
    File referenceDir

    # Directory for loftee
    File lofteeDir

    # Reference FASTA file - hg37/38
    File referenceFasta

    # CAD score files and associated index files
    File cadScores
    File cadScoresIndex

    ### Prepare percentiles ###
    Array[String] infoFields

    String vcf_basename = basename(input_vcf, ".vcf.gz")
  }

scatter (chromosome in chromosomes ) {
	call VCFsplitter {
		input:
			input_vcf = input_vcf,
			input_vcf_index = input_vcf_index,
			chromosome = chromosome
		}

	call vcfPercentilesPreparation.prepareVCFPercentiles as prepareVCFs {
		input:
			input_vcf = VCFsplitter.output_vcf,
			input_vcf_index = VCFsplitter.output_vcf_index,
			samplesFile = samplesFile,
			referenceDir = referenceDir, 
			lofteeDir = lofteeDir,
			referenceFasta = referenceFasta,
			cadScores = cadScores,
			cadScoresIndex = cadScoresIndex,
			infoFields = infoFields
		}

    call cramPreparation.prepareCram as prepareCRAMs {
      input:
        chromosome = chromosome,
        chromosomeVCF = prepareVCFs.output_vcf,
        samplesFile = samplesFile,
        referenceFasta = referenceFasta,
        sampleLocationFile = sampleLocationFile
    }
  # Close per chromosome scatter
  }

  output {
    Array[File] output_vcfs = prepareVCFs.output_vcf
    Array[File] output_vcfs_indices = prepareVCFs.output_vcf_index
    Array[File] out_metrics_files = prepareVCFs.out_metrics
    Array[File] out_crams = prepareCRAMs.combined_cram_result
    Array[File] out_crais = prepareCRAMs.combined_cram_result_index
  }
# Close workflow
}  


# Tasks
task VCFsplitter {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index

    String chromosome
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")

  command {
    set -e
    bcftools view -r ~{chromosome}:10000000-11000000 ~{input_vcf} -Oz -o ~{chromosome}.~{vcf_basename}.vcf.gz
    bcftools index -t ~{chromosome}.~{vcf_basename}.vcf.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    maxRetries: 3
    requested_memory_mb_per_core: 2000
    cpu: 3
    runtime_minutes: 180
  }
  output {
    File output_vcf = "~{vcf_basename}.~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{vcf_basename}.~{vcf_basename}.vcf.gz.tbi"
  }
}