## Copyright CMG@KIGM, Peter Juvan & Ales Maver

## CONSIDER:
## zcat vcfFilterNormMerge/e4eea4d9-2ea9-4187-b789-c91355bc3877/call-RunVEP/shard-10/execution/chr1__6746294_7746294.SGP9427_nrmFlt0.1Mrgd_ClinVar_vep.vcf.gz | grep dbNSFP
## cat vcfFilterNormMerge/e5698468-c4c8-44b2-a89f-f4174d2d1078/call-RunVEP/shard-145/execution/chr1__155184599_156184599.CMG14137_nrmFlt0.1Mrgd_ClinVar_vep.vcf.gz_warnings.txt
##  WARNING: Transcript-assembly mismatch in rs914616
##  WARNING: Transcript-assembly mismatch in chr1_155324912_T/G
##  WARNING: Transcript-assembly mismatch in chr1_155324912_T/G
##  WARNING: Transcript-assembly mismatch in rs145411349
##  WARNING: Transcript-assembly mismatch in rs145411349
## FIX:
##  bcftools +fixref $PVCF -- -m flip -f wgs_reference/Homo_sapiens_assembly38.fasta -i DPSNP/All_20180418_chr.vcf.gz

version 1.0

import "./vcfTasks.wdl" as vcfTasks

struct VEPReferences {
  String cache_dir
  String plugins_dir
  File dbNSFP_vcf
  File dbNSFP_vcf_index
  File dbNSFP_vcf_readme
  String loftee_data_dir
  String AlphaMissense_data_dir
}

struct AnnotationFields {
  String ClinVar
  String dbNSFP
}

# WORKFLOW DEFINITION 
workflow vcfAnnotate {
  input {
    File input_vcf
    File input_vcf_index

    File interval_list
    Int? thinning_parameter
    Int  scatter_region_size = 300000

    Int threads = 4   # use even numbers because slurm floors cpu to even numbers, but not total memory

    Boolean annotate_with_clinvar = true

    #### TODO IMPLEMENT
    Boolean annotate_with_dbnsfp = true
    Boolean annotate_with_alphamissense = true
    Boolean annotate_with_loftee = true

    # VEP references
    VEPReferences vep_ref

    ## VEP annotations
    # String DBNSFP_ANNFIELDS_DEFAULT="1000Gp3_AC,1000Gp3_EUR_AC,CADD_phred,ESP6500_AA_AC,ESP6500_EA_AC,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate"
    # String DBNSFP_ANNFIELDS_PRED=   "MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence"
    # String DBNSFP_ANNFIELDS_GNOMAD= "gnomAD_exomes_AC,gnomAD_exomes_nhomalt,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_nhomalt,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_NFE_nhomalt"
    # #String DBNSFP_ANNFIELDS_CLINVAR="clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id"
    AnnotationFields annotation_fields = {
      "ClinVar": "CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS",
      "dbNSFP": "1000Gp3_AC,1000Gp3_EUR_AC,CADD_phred,ESP6500_AA_AC,ESP6500_EA_AC,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,gnomAD_exomes_AC,gnomAD_exomes_nhomalt,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_nhomalt,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_NFE_nhomalt"
    }

    # Output
    String output_vcf_basename = basename(input_vcf, ".vcf.gz")
  }

  call GetClinVarVCF

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

    if ( annotate_with_clinvar ) {
      call AnnotateWithClinVarVCF {
        input:
          input_vcf = input_vcf,
          input_vcf_index = input_vcf_index,
          threads = threads,
          annotation_vcf = GetClinVarVCF.output_vcf,
          annotation_vcf_index = GetClinVarVCF.output_vcf_index,
          region = region,
          annotation_fields = annotation_fields.ClinVar,
          output_basename = sub(sub(region, "-", "_"), ":", "__") + ".anClinVar." + output_vcf_basename
      }
    }

    call VEP {
      input:
        input_vcf = select_first([AnnotateWithClinVarVCF.output_vcf, input_vcf]),
        input_vcf_index = select_first([AnnotateWithClinVarVCF.output_vcf_index, input_vcf_index]),
        cpus = if threads < 6 then 6 else threads,
        vep_ref = vep_ref,
        #vep_ref = vep_ref_split,
        annotation_fields = annotation_fields.dbNSFP,
        output_basename = sub(sub(region, "-", "_"), ":", "__") + ".anVEP." + output_vcf_basename
      }

  } # Close scatter region

  call vcfTasks.concatVcf {
    input:
      input_vcfs = VEP.output_vcf,
      input_vcfs_indices = VEP.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

  output {
    File output_vcf = concatVcf.output_vcf
    File output_vcfs_indices = concatVcf.output_vcf_index
  }

} # Close workflow


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
    cpu: 4
    runtime_minutes: 20
  }

  output {
    File output_vcf = "clinvar_fixed.vcf.gz"
    File output_vcf_index = "clinvar_fixed.vcf.gz.csi"
  }
} 

##############################
task AnnotateWithClinVarVCF {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    Int threads = 4
    File annotation_vcf
    File annotation_vcf_index

    String region
    String annotation_fields
    String output_basename = basename(input_vcf, ".vcf.gz")
  }

  command <<<
    set -e
    bcftools annotate -r ~{region} -a ~{annotation_vcf} -c ~{annotation_fields} ~{input_vcf} -Oz -o ~{output_basename}_ClinVar.vcf.gz
    bcftools index -t ~{output_basename}_ClinVar.vcf.gz --threads ~{threads}

  >>>

  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 2
    runtime_minutes: 10
  }

  output {
    File output_vcf = "~{output_basename}_ClinVar.vcf.gz"
    File output_vcf_index = "~{output_basename}_ClinVar.vcf.gz.tbi"
  }
}

##############################
## Plugin LoF requires input_vcf_index to be in .tbi format
## --fork chould not be used @ Vega
## recommended runtime cpu: 6 (estimated by mem usage for SGP VCF consisting 9425 samples and scatter_region_size = 300000)
## To consider: --buffer_size 50 (default 5000) will use less memory
task VEP {
  input {
    File input_vcf
    File input_vcf_index
    Int cpus = 6
    VEPReferences vep_ref
    String annotation_fields
    String output_basename = basename(input_vcf, ".vcf.gz")
  }

  command <<<
    vep -i ~{input_vcf} \
      -o ~{output_basename}_VEP.vcf.gz \
      --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
      --cache --merged --dir_cache ~{vep_ref.cache_dir} \
      --assembly GRCh38 \
      --everything \
      --flag_pick \
      --allele_number \
      --nearest symbol \
      --no_stats \
      --dir_plugins ~{vep_ref.plugins_dir} \
      --plugin dbNSFP,~{vep_ref.dbNSFP_vcf},~{annotation_fields} \
      --plugin LoF,loftee_path:~{vep_ref.plugins_dir},human_ancestor_fa:~{vep_ref.loftee_data_dir}/human_ancestor.fa.gz,conservation_file:~{vep_ref.loftee_data_dir}/loftee.sql,gerp_bigwig:~{vep_ref.loftee_data_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
      --plugin AlphaMissense,file=~{vep_ref.AlphaMissense_data_dir}/AlphaMissense_hg38.tsv.gz
      #--use_given_ref \

    tabix --force --preset vcf ~{output_basename}_VEP.vcf.gz
  >>>

  runtime {
      docker: "ensemblorg/ensembl-vep:latest"
      requested_memory_mb_per_core: 2000
      cpu: cpus
      runtime_minutes: 60
  }

  output {
      File output_vcf = "~{output_basename}_VEP.vcf.gz"
      File output_vcf_index = "~{output_basename}_VEP.vcf.gz.tbi"
  }
}
