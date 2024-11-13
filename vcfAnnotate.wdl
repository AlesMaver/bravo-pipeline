version 1.0
## Copyright CMG@KIGM, Peter Juvan & Ales Maver

import "./vcfTasks.wdl" as vcfTasks

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
    Boolean annotate_with_dbnsfp = true
    Boolean annotate_with_alphamissense = true
    Boolean annotate_with_loftee = true

    String ClinVar_annotation_fields = "CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS"

    # VEP config
#    String VEP_plugins_dir = "/plugins"
#    File dbNSFP_vcf
#    File dbNSFP_vcf_index
#    String dbNSFP_annotation_fields
#
#      DBNSFP_ANNFIELDS_DEFAULT_VEP="1000Gp3_AC,1000Gp3_EUR_AC,CADD_phred,ESP6500_AA_AC,ESP6500_EA_AC,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate"
#      DBNSFP_ANNFIELDS_PRED_VEP="MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence"
#      DBNSFP_ANNFIELDS_GNOMAD_VEP="gnomAD_exomes_AC,gnomAD_exomes_nhomalt,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_nhomalt,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_NFE_nhomalt"
#      #DBNSFP_ANNFIELDS_CLINVAR_VEP="clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id"
#      #DBNSFP_ANNFIELDS_VEP="$DBNSFP_ANNFIELDS_DEFAULT_VEP,$DBNSFP_ANNFIELDS_PRED_VEP,$DBNSFP_ANNFIELDS_GNOMAD_VEP,$DBNSFP_ANNFIELDS_CLINVAR_VEP"
#      DBNSFP_ANNFIELDS_VEP="$DBNSFP_ANNFIELDS_DEFAULT_VEP,$DBNSFP_ANNFIELDS_PRED_VEP,$DBNSFP_ANNFIELDS_GNOMAD_VEP"
#
#    File AlphaMissense_tsv_gz = "/plugins_data/AlphaMissense_hg38.tsv.gz"
#    File AlphaMissense_tsv_gz_tbi = "/plugins_data/AlphaMissense_hg38.tsv.gz.tbi"
#
#    File loftee_gerp_conservation_scores_bw = "/plugins_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
#    File loftee_human_ancestor_fa_gz = "/plugins_data/human_ancestor.fa.gz"
#    File loftee_human_ancestor_fa_gz_fai = "/plugins_data/human_ancestor.fa.gz.fai"
#    File loftee_human_ancestor_fa_gz_gzi = "/plugins_data/human_ancestor.fa.gz.gzi"
#    File loftee_sql = "/plugins_data/loftee.sql"

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
      call AnnotateWithVCF {
        input:
          input_vcf = input_vcf,
          input_vcf_index = input_vcf_index,
          annotation_vcf = GetClinVarVCF.output_vcf,
          annotation_vcf_index = GetClinVarVCF.output_vcf_index,
          region = region,
          annotation_fields = ClinVar_annotation_fields,
          output_basename = output_vcf_basename
      }
    }

    call RunVEP {
      input:
        input_vcf = select_first([AnnotateWithVCF.output_vcf, input_vcf]),
        input_vcf_index = select_first([AnnotateWithVCF.output_vcf_index, input_vcf_index]),
        cpus = if threads < 24 then 24 else threads,
        output_basename = output_vcf_basename
      }

  } # Close scatter region

  output {
    Array[File] output_vcf = RunVEP.output_vcf
    Array[File] output_vcfs_indices = RunVEP.output_vcf_index
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
    cpu: 3
    runtime_minutes: 20
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

    String region

    String annotation_fields
  
    String output_basename = basename(input_vcf, ".vcf.gz")
  }


  command <<<
    set -e
    bcftools annotate -r ~{region} -a ~{annotation_vcf} -c ~{annotation_fields} ~{input_vcf} -Oz -o ~{output_basename}_ClinVar.vcf.gz --write-index
  >>>

  runtime {
    docker: "alesmaver/bcftools"
    requested_memory_mb_per_core: 2000
    cpu: 3
    runtime_minutes: 59
  }

  output {
    File output_vcf = "~{output_basename}_ClinVar.vcf.gz"
    File output_vcf_index = "~{output_basename}_ClinVar.vcf.gz.csi"
  }
}


##############################
task RunVEP {
    input {
      File input_vcf
      File input_vcf_index
      Int cpus = 12
      String output_basename = basename(input_vcf, ".vcf.gz")
    }


    command <<<
      PERL5LIB=:\$PERL5LIB:/opt/vep/.vep/Plugins/loftee

      DBNSFP_ANNFIELDS_DEFAULT_VEP="1000Gp3_AC,1000Gp3_EUR_AC,CADD_phred,ESP6500_AA_AC,ESP6500_EA_AC,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate"
      DBNSFP_ANNFIELDS_PRED_VEP="MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence"
      DBNSFP_ANNFIELDS_GNOMAD_VEP="gnomAD_exomes_AC,gnomAD_exomes_nhomalt,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_nhomalt,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_NFE_nhomalt"
      #DBNSFP_ANNFIELDS_CLINVAR_VEP="clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id"
      #DBNSFP_ANNFIELDS_VEP="$DBNSFP_ANNFIELDS_DEFAULT_VEP,$DBNSFP_ANNFIELDS_PRED_VEP,$DBNSFP_ANNFIELDS_GNOMAD_VEP,$DBNSFP_ANNFIELDS_CLINVAR_VEP"
      DBNSFP_ANNFIELDS_VEP="$DBNSFP_ANNFIELDS_DEFAULT_VEP,$DBNSFP_ANNFIELDS_PRED_VEP,$DBNSFP_ANNFIELDS_GNOMAD_VEP"

      vep -i ~{input_vcf} \
        -o ~{output_basename}_vep.vcf.gz \
        --fork "~{cpus}" --cache --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
        --assembly GRCh38 \
        --everything \
        --flag_pick \
        --allele_number \
        --dir_cache /opt/vep/.vep \
        --merged \
        --nearest symbol \
        --no_stats \
        --plugin dbNSFP,/opt/vep/.vep/dbNSFP/dbNSFPv4.9a_custombuild.gz,$DBNSFP_ANNFIELDS_VEP \
        --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg38.tsv.gz \
        --plugin LoF,loftee_path:/opt/vep/.vep/Plugins/loftee/,human_ancestor_fa:/opt/vep/.vep/Plugins/loftee/data/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/Plugins/loftee/data/loftee.sql,gerp_bigwig:/opt/vep/.vep/Plugins/loftee/data/gerp_conservation_scores.homo_sapiens.GRCh38.bw \

      tabix -p vcf ~{output_basename}_vep.vcf.gz
    >>>

    runtime {
        docker: "peterjuv/vep_docker:latest"
        requested_memory_mb_per_core: 2000
        cpu: cpus
        runtime_minutes: 360
    }

    output {
        File output_vcf = "~{output_basename}_vep.vcf.gz"
        File output_vcf_index = "~{output_basename}_vep.vcf.gz.tbi"
    }
}
