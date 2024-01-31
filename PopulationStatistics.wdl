version 1.0
## Copyright CMG@KIGM, Ales Maver

#import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/common/annotation/CreateGenesBed.wdl" as CreateGenesBed
import "http://wdl_server/CMGpipeline/common/annotation/CreateGenesBed.wdl" as CreateGenesBed


# WORKFLOW DEFINITION 
workflow PopulationStatistics {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 3000000

    Boolean annotate_with_clinvar = true

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
  
    Array[String]? panel_genes
    String panel_name = "PANEL"

    # Reference FASTA file - hg37/38
    File gencode_gff3
    File referenceFasta
  }

  String output_file_name = "MergedVariantTable_~{panel_name}.tab"

  call GetClinVarVCF

  if ( defined(panel_genes) ) {
    call CreateGenesBed.DownloadAndPrepareBed as CreateGenesBed {
      input:
        gencode_gff3 = gencode_gff3
    }

    call CreateGenesBed.FilterGenesBED as FilterGenesBED {
      input:
        genes_bed = CreateGenesBed.genes_bed,
        filter_genes = select_first([panel_genes])
      }
  }

  call ConvertIntervalListToBed {
    input:
      interval_list = interval_list
  }

  call SplitRegions {
    input:
      input_bed = select_first([FilterGenesBED.filtered_genes_bed, ConvertIntervalListToBed.converted_bed]),
      thinning_parameter = thinning_parameter,
      scatter_region_size = scatter_region_size
  }

  scatter (chromosome in SplitRegions.scatter_regions) {
    if ( annotate_with_clinvar ) {
      call AnnotateWithVCF {
        input:
          input_vcf = input_vcf,
          input_vcf_index = input_vcf_index,
          annotation_vcf = GetClinVarVCF.output_vcf,
          annotation_vcf_index = GetClinVarVCF.output_vcf_index,
          chromosome = chromosome,
          annotation_fields ="CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS"
      }
    }

  	call GenerateTable {
  		input:
  			input_vcf = select_first([AnnotateWithVCF.output_vcf, input_vcf]),
  			input_vcf_index = select_first([AnnotateWithVCF.output_vcf_index, input_vcf_index]),
  			chromosome = chromosome,
        referenceFasta = referenceFasta
  	    }

    call GenerateGeneFilteredTable {
        input:
          input_vcf = select_first([AnnotateWithVCF.output_vcf, input_vcf]),
          input_vcf_index = select_first([AnnotateWithVCF.output_vcf_index, input_vcf_index]),
          filter_genes = select_first([panel_genes, ["ACTB"]]),
          chromosome = chromosome,
          referenceFasta = referenceFasta
        }
    }


  call ConcatenateTabFiles {
    input:
      input_files = GenerateTable.output_tab,
      output_name = output_file_name
  }

  if ( defined(panel_genes) ) {
    call ConcatenateTabFiles as ConcatenateTabFiles2 {
      input:
        input_files = select_first([GenerateGeneFilteredTable.output_tab, GenerateTable.output_tab]),
        output_name = output_file_name
    }
  }

  if ( false ) {
    # Alternative gene filtration
    call SplitRegions as SplitRegions2 {
      input:
        input_bed = ConvertIntervalListToBed.converted_bed,
        thinning_parameter = thinning_parameter,
        scatter_region_size = scatter_region_size
    }

    scatter (chromosome in SplitRegions2.scatter_regions) {
      if ( annotate_with_clinvar ) {
        call AnnotateWithVCF as AnnotateWithVCF2  {
          input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            annotation_vcf = GetClinVarVCF.output_vcf,
            annotation_vcf_index = GetClinVarVCF.output_vcf_index,
            chromosome = chromosome,
            annotation_fields ="CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS"
        }
      }

      call GenerateGeneFilteredTable as GenerateGeneFilteredTable2 {
        input:
          input_vcf = select_first([AnnotateWithVCF2.output_vcf, input_vcf]),
          input_vcf_index = select_first([AnnotateWithVCF2.output_vcf_index, input_vcf_index]),
          filter_genes = select_first([panel_genes, ["ACTB"]]),
          chromosome = chromosome,
          referenceFasta = referenceFasta
          }
      }
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
task GetClinVarVCF {
  command <<<
    set -e
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/c691a9607e33337084afcc11372ba69ed5870178/references/rename_chrs
    bcftools annotate --rename-chrs rename_chrs clinvar.vcf.gz --write-index -Oz -o clinvar_fixed.vcf.gz
  >>>

  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }

  output {
    File output_vcf = "clinvar_fixed.vcf.gz"
    File output_vcf_index = "clinvar_fixed.vcf.gz.csi"
  }
} 

##############################
task AnnotateWithVCF {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index

    File annotation_vcf
    File annotation_vcf_index

    String chromosome

    String annotation_fields ="CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS"
  }

  command <<<
    set -e
    bcftools annotate -r ~{chromosome} -a ~{annotation_vcf} -c ~{annotation_fields} ~{input_vcf} -Oz -o vcf_annotated.vcf.gz --write-index
  >>>

  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }

  output {
    File output_vcf = "vcf_annotated.vcf.gz"
    File output_vcf_index = "vcf_annotated.vcf.gz.csi"
  }
}

##############################
task GenerateTable {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index

    Boolean Pathogenic = false
    Boolean ModerateOrHigh = true

    String chromosome
    File referenceFasta
  }

  String vcf_basename = basename(input_vcf, ".vcf.gz")
  String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

  command <<<
    set -e
    bcftools +split-vep -l ~{input_vcf} > available_info_fields.txt

    # bcftools view -r ~{chromosome} -t ~{chromosome} ~{input_vcf} | \
    # bcftools norm -m -any | \
    # bcftools +fill-tags | \
    # bcftools +split-vep -c gnomADg_AF:Float -s primary | \
    # bcftools +split-vep \
    #   -f '[%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%IMPACT\t%MANE_SELECT\t%CANONICAL\t%EXON\t%HGVSc\t%HGVSp\t%clinvar_clnsig\t%clinvar_review\t%CLNSIG\t%CLNDN\t%CLNSIGCONF\t%CLNREVSTAT\t%MutationTaster_pred\t%Polyphen2_HDIV_pred\t%SIFT_pred\t%PROVEAN_pred\t%LRT_pred\t%CADD_phred\t%REVEL_rankscore\t%REVEL_score\t%MetaSVM_pred\t%MetaRNN_pred\t%MutationAssessor_pred\t%GERP___RS\t%Interpro_domain\t%LoF\t%LoF_filter\t%LoF_flags\t%LoF_info\t%INFO/AF\t%INFO/AC\t%AC_Hom\t%AN\t%INFO/gnomADg_AF\t%SAMPLE\t%GT\n]' \
    #   -s primary \
    #   -i'( ( CLIN_SIG ~ "pathogenic/i" || ((clinvar_clnsig ~ "pathogenic/i") && (clinvar_clnsig !~ "conflicting/i"))) || ((IMPACT ~ "MODERATE/i") || (IMPACT ~ "HIGH/i")) ) && (gnomADg_AF<=0.01 || gnomADg_AF=".") && (AF<=0.01 || AF=".") && GT="alt"' \
    #   -d -H \
    # > ~{chromosome_filename}.~{vcf_basename}.tab

    bcftools view -r ~{chromosome} -t ~{chromosome} ~{input_vcf} | \
    bcftools norm -m -any | \
    bcftools +fill-tags | \
    bcftools +split-vep -c SYMBOL,HGVSc,HGVSp -s primary -p "primary_gene_vep_" | \
    bcftools +split-vep -c SYMBOL,HGVSc,HGVSp -s worst -p "worst_vep_" |
    bcftools +split-vep -c gnomADg_AF:Float -s primary | \
    bcftools +split-vep \
      -f '[%CHROM\t%POS\t%REF\t%ALT\t%primary_gene_vep_SYMBOL\t%primary_gene_vep_HGVSc\t%SYMBOL\t%IMPACT\t%MANE_SELECT\t%CANONICAL\t%EXON\t%HGVSc\t%HGVSp\t%worst_vep_SYMBOL\t%clinvar_clnsig\t%clinvar_review\t%CLNSIG\t%CLNDN\t%CLNSIGCONF\t%CLNREVSTAT\t%MutationTaster_pred\t%Polyphen2_HDIV_pred\t%SIFT_pred\t%PROVEAN_pred\t%LRT_pred\t%CADD_phred\t%REVEL_rankscore\t%REVEL_score\t%MetaSVM_pred\t%MetaRNN_pred\t%MutationAssessor_pred\t%GERP___RS\t%Interpro_domain\t%LoF\t%LoF_filter\t%LoF_flags\t%LoF_info\t%INFO/AF\t%INFO/AC\t%AC_Hom\t%AN\t%INFO/gnomADg_AF\t%SAMPLE\t%GT\t%AD\n]' \
      -s primary \
      -i'CLNSIG ~ "pathogenic/i" && GT="alt"' \
      -d -H \
    | uniq >> ~{chromosome_filename}.~{vcf_basename}.tab

  >>>
  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_tab = "~{chromosome_filename}.~{vcf_basename}.tab"
  }
}

##############################
task GenerateGeneFilteredTable {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index

    Array[String] filter_genes

    String chromosome
    File referenceFasta
  }

  command <<<
    set -e
    # Split gene list into smaller chunks
    split ~{write_lines(filter_genes)} genes_ -l 500 --numeric-suffixes -a 5

    # Process the chunks
    for file in genes_* 
      do
        genes=$(tr '\n' ',' < $file)
        echo Working on the gene list: $genes
        genefilter='primary_gene_vep_SYMBOL=="'$genes'"'

        bcftools view -r ~{chromosome} -t ~{chromosome} ~{input_vcf} | \
        bcftools norm -m -any | \
        bcftools +fill-tags | \
        bcftools +split-vep -i ${genefilter[@]} -c SYMBOL,HGVSc,HGVSp -s primary -p "primary_gene_vep_" | \
        bcftools +split-vep -c SYMBOL,HGVSc,HGVSp -s worst -p "worst_vep_" |
        bcftools +split-vep -c gnomADg_AF:Float -s primary | \
        bcftools +split-vep \
          -f '[%CHROM\t%POS\t%REF\t%ALT\t%primary_gene_vep_SYMBOL\t%primary_gene_vep_HGVSc\t%SYMBOL\t%IMPACT\t%MANE_SELECT\t%CANONICAL\t%EXON\t%HGVSc\t%HGVSp\t%worst_vep_SYMBOL\t%clinvar_clnsig\t%clinvar_review\t%CLNSIG\t%CLNDN\t%CLNSIGCONF\t%CLNREVSTAT\t%MutationTaster_pred\t%Polyphen2_HDIV_pred\t%SIFT_pred\t%PROVEAN_pred\t%LRT_pred\t%CADD_phred\t%REVEL_rankscore\t%REVEL_score\t%MetaSVM_pred\t%MetaRNN_pred\t%MutationAssessor_pred\t%GERP___RS\t%Interpro_domain\t%LoF\t%LoF_filter\t%LoF_flags\t%LoF_info\t%INFO/AF\t%INFO/AC\t%AC_Hom\t%AN\t%INFO/gnomADg_AF\t%SAMPLE\t%GT\n]' \
          -s primary \
          -i'( ( CLIN_SIG ~ "pathogenic/i" || ((clinvar_clnsig ~ "pathogenic/i") && (clinvar_clnsig !~ "conflicting/i"))) || ((IMPACT ~ "MODERATE/i") || (IMPACT ~ "HIGH/i")) ) && (gnomADg_AF<=0.2 || gnomADg_AF=".") && (AF<=0.01 || AF=".") && GT="alt"' \
          -d -H \
        | uniq >> output_rare_pathogenic.tab

        bcftools view -r ~{chromosome} -t ~{chromosome} ~{input_vcf} | \
        bcftools norm -m -any | \
        bcftools +fill-tags | \
        bcftools +split-vep -i ${genefilter[@]} -c SYMBOL,HGVSc,HGVSp -s primary -p "primary_gene_vep_" | \
        bcftools +split-vep -c SYMBOL,HGVSc,HGVSp -s worst -p "worst_vep_" |
        bcftools +split-vep -c gnomADg_AF:Float -s primary | \
        bcftools +split-vep \
          -f '[%CHROM\t%POS\t%REF\t%ALT\t%primary_gene_vep_SYMBOL\t%primary_gene_vep_HGVSc\t%SYMBOL\t%IMPACT\t%MANE_SELECT\t%CANONICAL\t%EXON\t%HGVSc\t%HGVSp\t%worst_vep_SYMBOL\t%clinvar_clnsig\t%clinvar_review\t%CLNSIG\t%CLNDN\t%CLNSIGCONF\t%CLNREVSTAT\t%MutationTaster_pred\t%Polyphen2_HDIV_pred\t%SIFT_pred\t%PROVEAN_pred\t%LRT_pred\t%CADD_phred\t%REVEL_rankscore\t%REVEL_score\t%MetaSVM_pred\t%MetaRNN_pred\t%MutationAssessor_pred\t%GERP___RS\t%Interpro_domain\t%LoF\t%LoF_filter\t%LoF_flags\t%LoF_info\t%INFO/AF\t%INFO/AC\t%AC_Hom\t%AN\t%INFO/gnomADg_AF\t%SAMPLE\t%GT\n]' \
          -s primary \
          -i'CLNSIG ~ "pathogenic/i" && GT="alt"' \
          -d -H \
        | uniq >> output.tab
    done
    
  >>>
  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_tab = "output.tab"
  }
}

##############################
task ConcatenateTabFiles {
  input {
    Array[File] input_files
    String output_name = "MergedVariantTable"
  }
  command <<<
    cat $(cat ~{write_lines(input_files)}) | awk '!seen[$0]++ || !/^#/' | uniq > ~{output_name}
  >>>

  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    #runtime_minutes: 180
  }
  output {
    File output_tab = "~{output_name}"
  }
}