## Copyright CMG@KIGM, Peter Juvan & Ales Maver

## CONSIDER:
##  bcftools +fixref $PVCF -- -m flip -f wgs_reference/Homo_sapiens_assembly38.fasta -i DPSNP/All_20180418_chr.vcf.gz
##  update AnnotateWithClinVarVCF: split to VCFsplitSubset + AnnotateWithClinVarVCF
##  VCFindex after VEP (don't use tabix within VEP)

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
    Int  scatter_region_size = 1000000 # smaller regions result in empty VCFs that will fail at concatenate

    Int threads = 4   # use even numbers because slurm floors cpu to even numbers, but not total memory

    Boolean annotate_with_clinvar = true
    #### TODO IMPLEMENT
    Boolean annotate_with_dbnsfp = true
    Boolean annotate_with_alphamissense = true
    Boolean annotate_with_loftee = true

    # VEP
    VEPReferences vep_ref
    String assembly = "GRCh38"
    Int buffer_size = 5000

    ## Annotations
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
        assembly = assembly,
        buffer_size = buffer_size,
        annotate_with_dbnsfp = annotate_with_dbnsfp,
        annotate_with_alphamissense = annotate_with_alphamissense,
        annotate_with_loftee = annotate_with_loftee,
        annotation_fields = annotation_fields.dbNSFP,
        output_basename = sub(sub(region, "-", "_"), ":", "__") + ".anVEP." + output_vcf_basename
      }

    # sort after VEP to avoid indexing error after concat, e.g.:
    #  [E::hts_idx_push] Unsorted positions on sequence #9: 133220600 followed by 133220598
    # TODO: REMOVE, NOT NEEDED, WE STILL GET THE SAME ERROR
    call vcfTasks.sortIdxVcf {
      input:
        input_vcf = VEP.output_vcf,
        input_vcf_index = VEP.output_vcf_index,
        output_name = sub(sub(region, "-", "_"), ":", "__") + ".sorted." + output_vcf_basename,
        threads = threads
    }

#    call vcfTasks.VCFindex as sortVcf_index {
#      input:
#        input_vcf = sortVcf.output_vcf,
#        threads = threads
#    }

  } # Close scatter region

  call vcfTasks.concatIdxVcf {
    input:
      input_vcfs = sortIdxVcf.output_vcf,
      input_vcfs_indices = sortIdxVcf.output_vcf_index,
      output_name = output_vcf_basename,
      threads = threads
  }

#  call vcfTasks.VCFindex as concatVcf_index{
#    input:
#      input_vcf = concatVcf.output_vcf,
#      threads = threads
#  }

  output {
    File output_vcf = concatIdxVcf.output_vcf
    File output_vcfs_indices = concatIdxVcf.output_vcf_index
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
    docker: "peterjuv/bcftools"
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
    bcftools annotate --threads ~{threads} -r ~{region} -a ~{annotation_vcf} -c ~{annotation_fields} ~{input_vcf} -Oz -o ~{output_basename}_ClinVar.vcf.gz --write-index=tbi
  >>>

  runtime {
    docker: "peterjuv/bcftools"
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
## recommended runtime cpu: 6 (estimated by mem usage for SGP VCF consisting 9425 samples and scatter_region_size = 300000)
## 
## Options disabled:
## --fork: should not be used @ Vega, makes it crash with ERROR: Forked process(es) died: read-through of cross-process communication detected
## Options used:
##  --flag_pick: Instead of choosing one block and removing the others, this option adds a flag "PICK=1" to picked annotation block, allowing you to easily filter on this
## Options to consider:
##  --use_given_ref: Using --bam or a BAM-edited RefSeq cache by default enables --use_transcript_ref; add this flag to override this behaviour and use the provided reference allele from the input. 
##  --shift_hgvs 0: Was used in Bravo vcfPercentilesPreparation.wdl 
task VEP {
  input {
    File input_vcf
    File input_vcf_index
    Int cpus = 6
    VEPReferences vep_ref
    String assembly = "GRCh38"
    Int buffer_size = 5000
    Boolean annotate_with_dbnsfp = true
    Boolean annotate_with_loftee = true
    Boolean annotate_with_alphamissense = true
    String? annotation_fields
    String output_basename = basename(input_vcf, ".vcf.gz")
  }

  String arg_dbnsfp = if annotate_with_dbnsfp && defined(annotation_fields) then '--plugin dbNSFP,' + vep_ref.dbNSFP_vcf + ',' + annotation_fields else ''
  String arg_loftee = if annotate_with_loftee then '--plugin LoF,loftee_path:' + vep_ref.plugins_dir + ',human_ancestor_fa:' + vep_ref.loftee_data_dir + '/human_ancestor.fa.gz,conservation_file:' + vep_ref.loftee_data_dir + '/loftee.sql,gerp_bigwig:' + vep_ref.loftee_data_dir + '/gerp_conservation_scores.homo_sapiens.GRCh38.bw' else ''
  String arg_am = if annotate_with_alphamissense then '--plugin AlphaMissense,file=' + vep_ref.AlphaMissense_data_dir + '/AlphaMissense_hg38.tsv.gz' else ''

  command <<<
    vep -i ~{input_vcf} \
      -o ~{output_basename}_VEP.vcf.gz \
      --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
      --cache --merged --dir_cache ~{vep_ref.cache_dir} \
      --assembly ~{assembly} \
      --everything \
      --flag_pick \
      --allele_number \
      --nearest symbol \
      --no_stats \
      --dir_plugins ~{vep_ref.plugins_dir} \
      ~{arg_dbnsfp} \
      ~{arg_loftee} \
      ~{arg_am} \
      --buffer_size ~{buffer_size}

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

###################################
## VEP for Bravo
task VEP_Bravo {
  input {
    File input_vcf
    File input_vcf_index
    Int cpus = 6
    VEPReferences vep_ref
    String assembly = "GRCh38"
    Int buffer_size = 5000
    Boolean annotate_with_dbnsfp = true
    Boolean annotate_with_loftee = true
    Boolean annotate_with_alphamissense = true
    String? annotation_fields
    String output_basename = basename(input_vcf, ".vcf.gz")
  }

  String arg_dbnsfp = if annotate_with_dbnsfp && defined(annotation_fields) then '--plugin dbNSFP,' + vep_ref.dbNSFP_vcf + ',' + annotation_fields else ''
  String arg_loftee = if annotate_with_loftee then '--plugin LoF,loftee_path:' + vep_ref.plugins_dir + ',human_ancestor_fa:' + vep_ref.loftee_data_dir + '/human_ancestor.fa.gz,conservation_file:' + vep_ref.loftee_data_dir + '/loftee.sql,gerp_bigwig:' + vep_ref.loftee_data_dir + '/gerp_conservation_scores.homo_sapiens.GRCh38.bw' else ''
  String arg_am = if annotate_with_alphamissense then '--plugin AlphaMissense,file=' + vep_ref.AlphaMissense_data_dir + '/AlphaMissense_hg38.tsv.gz' else ''

  command <<<
    vep -i ~{input_vcf} \
      -o ~{output_basename}_VEP.vcf.gz \
      --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
      --cache --dir_cache ~{vep_ref.cache_dir} \
      --assembly ~{assembly} \
      --allele_number \
      --no_stats \
      --dir_plugins ~{vep_ref.plugins_dir} \
      ~{arg_dbnsfp} \
      ~{arg_loftee} \
      ~{arg_am} \
      --buffer_size ~{buffer_size} \
      --sift b \
      --polyphen b \
      --ccds \
      --uniprot \
      --hgvs \
      --symbol \
      --numbers \
      --domains \
      --regulatory \
      --canonical \
      --protein \
      --biotype \
      --af \
      --af_1kg \
      --pubmed \
      --shift_hgvs 0

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
