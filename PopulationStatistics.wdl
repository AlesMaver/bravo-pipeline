version 1.0
## Copyright CMG@KIGM, Ales Maver

import "./vcfTasks.wdl" as vcfTasks

# WORKFLOW DEFINITION 
workflow PopulationStatistics {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 1000000

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  
    # Reference FASTA file - hg37/38
    File referenceFasta
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
  	call vcfTasks.GenerateTable {
  		input:
  			input_vcf = input_vcf,
  			input_vcf_index = input_vcf_index,
  			chromosome = chromosome
  	    }
    }

  call vcfTasks.ConcatenateTabFiles {
    input:
      input_files = GenerateTable.output_tab
  }

  output {
    File output_tab = ConcatenateTabFiles.output_tab
  }

} # Close workflow
