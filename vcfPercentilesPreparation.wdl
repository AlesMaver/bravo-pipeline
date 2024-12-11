version 1.0

import "./vcfAnnotate.wdl" as vcfAnnotate
import "./vcfTasks.wdl" as vcfTasks


workflow prepareVCFPercentiles {
    ### Prepare VCF Inputs ###
    input {
        # Chromosome VCF file
        File input_vcf
        File input_vcf_index

        # File for samples
        File samplesFile

        # Size of memory buffer to use
        Int buffer_size  = 100000

        # VEP
        VEPReferences vep_ref
        String assembly = "GRCh38"

        # CAD score files and associated index files
        File cadScores
        File cadScoresIndex

        Int threads = 4
        Int numberPercentiles = 10
        String description = "Description"
    }

    String vcf_basename = basename(input_vcf, ".vcf.gz")

    ###############
    # Prepare VCF #
    ###############

    call computeAlleleCountsAndHistograms {
        input: 
            chromosomeVCF = input_vcf,
            samplesFile = samplesFile
    }

#    call AddOriginalVCFAnnotations {
#        input:
#            original_input_vcf = input_vcf,
#            original_input_vcf_index = input_vcf_index,
#            analysed_vcf = computeAlleleCountsAndHistograms.out,
#            threads = threads
#    }

#    call variantEffectPredictor {
#        input: 
#            chromosomeVCF = AddOriginalVCFAnnotations.output_vcf,
#            assembly = assembly,
#            buffer_size = buffer_size,
#            #referenceDir = referenceDir,
#            referenceFasta = referenceFasta,
#            #lofteeDir = lofteeDir
#            forks = threads
#    }

    call vcfTasks.VCFindex {
        input:
            input_vcf = computeAlleleCountsAndHistograms.out,
            threads = threads
    }

    call vcfTasks.VCFdropGeno {
        input: 
            input_vcf = VCFindex.output_vcf,
            input_vcf_index = VCFindex.output_vcf_index,
            output_name = vcf_basename + "_droppedGeno",
            threads = threads
    }

    call vcfAnnotate.VEP {
        input: 
            input_vcf = VCFdropGeno.output_vcf,
            input_vcf_index = VCFdropGeno.output_vcf_index,
            cpus = if threads < 6 then 6 else threads,
            vep_ref = vep_ref,
            assembly = assembly,
            buffer_size = buffer_size,
            annotate_with_dbnsfp = false,
            annotate_with_alphamissense = false,
            annotate_with_loftee = true
    }

    call addCaddScores {
        input: 
            chromosomeVCF = VEP.output_vcf,
            cadScores = cadScores,
            cadScoresIndex = cadScoresIndex
    }

    #######################
    # Prepare percentiles #
    #######################

    #scatter (field in infoFields) {
    #    call computePercentiles {
    #        input: chromosomeVCF = addCaddScores.out,
    #            infoField = field,
    #            threads = threads,
    #            numberPercentiles = numberPercentiles,
    #            description = description
    #    }
    #}
    #call addPercentiles {
    #    input: 
    #        chromosomeVCF = addCaddScores.out,
    #        chromosomeVCFIndex = addCaddScores.out_index,
    #        variantPercentiles = computePercentiles.outVariantPercentile,
    #        vcf_basename = vcf_basename
    #
    #}

  # Outputs that will be retained when execution is complete
  output {
    File output_annotated_vcf = addCaddScores.out
    File output_annotated_vcf_index = addCaddScores.out_index
    #File output_vcf = addPercentiles.out
    #File output_vcf_index = addPercentiles.out_index
    #Array[File] out_metrics = computePercentiles.outAllPercentiles
  }

}

task computeAlleleCountsAndHistograms {
    input {
        File chromosomeVCF
        File samplesFile
    }

    command {
        /srv/data/bravo_data_prep/data_prep/cpp_tools/cget/bin/ComputeAlleleCountsAndHistograms -i ~{chromosomeVCF} -s ~{samplesFile} -o computeAlleleCtHst.vcf.gz  #--fields AS_InbreedingCoeff AS_QD ExcessHet FS InbreedingCoeff MLEAF MQ QD RAW_MQ SOR VQSLOD
    }
    output {
        File out = "computeAlleleCtHst.vcf.gz"
    }
    runtime {
        # Docker "alesmaver/bravo-pipeline-sgp:latest" prepared using the following steps
        # This docker was necessary for the carry-over of field annotations to the resulting VCF files
        # docker run -it statgen/bravo-pipeline
        # git clone https://github.com/statgen/bravo_data_prep.git
        # cd bravo_data_prep/data_prep/cpp_tools/
        # cget install .
        #
        # Also installing python3
        # apt-get install python3-pip
        # pip3 install pysam

        docker: "alesmaver/bravo-pipeline-sgp:latest"
        runtime_minutes: 60
    }
}

## Generate table of variants for interpretation
#task AddOriginalVCFAnnotations {
#    input {
#      File original_input_vcf
#      File original_input_vcf_index
#      File analysed_vcf
#      Int threads
#    }
#  
#  command <<<
#  set -e
#    bcftools index -t ~{analysed_vcf}
#    bcftools annotate  --threads ~{threads} -a ~{original_input_vcf} -c +INFO ~{analysed_vcf} -Oz -o output.vcf.gz
#    bcftools index -t output.vcf.gz
#  >>>
#
#  runtime {
#    docker: "dceoy/bcftools"
#    requested_memory_mb_per_core: 2000
#    cpu: threads
#    runtime_minutes: 60
#  }
#  output {
#    File output_vcf = "output.vcf.gz"
#    File output_vcf_index = "output.vcf.gz.tbi"
#  }
#}

#task variantEffectPredictor {
#    input {
#        File chromosomeVCF
#        String assembly
#        Int buffer_size
#        #File referenceDir
#        File referenceFasta
#        #File lofteeDir
#        Int forks
#
#    }
#
#    command <<<
#        PERL5LIB=$PERL5LIB:/opt/vep/plugins/loftee/
#
#        vep -i ~{chromosomeVCF} \
#        --plugin LoF,loftee_path:/opt/vep/plugins/loftee/,human_ancestor_fa:/opt/vep/plugins/loftee/data/human_ancestor.fa.gz,conservation_file:/opt/vep/plugins/loftee/data/phylocsf_gerp.sql  \
#        --dir_cache /opt/vep/.vep/ \
#        --fasta ~{referenceFasta} \
#        --assembly ~{assembly} \
#        --cache \
#        --offline \
#        --vcf \
#        --sift b \
#        --polyphen b \
#        --ccds \
#        --uniprot \
#        --hgvs \
#        --symbol \
#        --numbers \
#        --domains \
#        --regulatory \
#        --canonical \
#        --protein \
#        --biotype \
#        --af \
#        --af_1kg \
#        --pubmed \
#        --shift_hgvs 0 \
#        --allele_number \
#        --format vcf \
#        --force \
#        --buffer_size ~{buffer_size} \
#        --compress_output bgzip \
#        --no_stats \
#        --fork ~{forks} \
#        --dir_plugins /opt/vep/plugins/loftee/ \
#        -o variantEP.vcf.gz
#    >>>
#    output {
#        File out = "variantEP.vcf.gz"
#    }
#    runtime {
#        #docker: "ensemblorg/ensembl-vep:release_95.1"
#        #docker: "ensemblorg/ensembl-vep:release_106.1"
#        docker: "alesmaver/vep:testing"
#        cpu: forks # "1" # changed in order to increase memory, see https://github.com/Ensembl/ensembl-vep/issues/150
#        runtime_minutes: 60
#    }
#
#}

task addCaddScores {
    input {
        File chromosomeVCF
        File cadScores
        File cadScoresIndex
    }

    command {
        add_cadd_scores.py -i ~{chromosomeVCF} -c ~{cadScores} -o annotated.vcf.gz
        tabix annotated.vcf.gz
    }
    output {
        File out = "annotated.vcf.gz"
        File out_index = "annotated.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "2"
        runtime_minutes: 10
    }
}

task computePercentiles {
    input {
        File chromosomeVCF
        String infoField
        Int threads
        Int numberPercentiles
        String description
    }

    command <<<
        mkdir -p metrics
        ComputePercentiles -i ~{chromosomeVCF} \
        -m ~{infoField} \
        -t ~{threads} \
        -p ~{numberPercentiles} \
        -d ~{description}_~{infoField} \
        -o metrics/~{infoField}

        tabix metrics/~{infoField}.variant_percentile.vcf.gz
    >>>
    output {
        File outAllPercentiles = "metrics/~{infoField}.all_percentiles.json.gz"
        File outVariantPercentile = "metrics/~{infoField}.variant_percentile.vcf.gz"
        File outVariantPercentileIndex = "metrics/~{infoField}.variant_percentile.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: threads
        runtime_minutes: 120
    }
}

task addPercentiles {
    input {
        File chromosomeVCF
        File chromosomeVCFIndex
        Array[File] variantPercentiles
        Array[File] variantPercentilesIndex
        Array[File] metricJSONs
        String vcf_basename
    }

    command <<<
        mkdir -p vcf
        mkdir -p metrics
        add_percentiles.py -i ~{chromosomeVCF} -p ~{sep=' ' variantPercentiles} -o vcf/~{vcf_basename}.percentiles.vcf.gz
        tabix vcf/~{vcf_basename}.percentiles.vcf.gz
        ## create mertics.json
        echo -n "[" > metrics/metrics.json
        zcat ~{sep=" " metricJSONs} | tr "\n" ","  >> metrics/metrics.json
        sed '$ s/.$//' -i metrics/metrics.json
        echo -n "]" >> metrics/metrics.json
        ## fix metrics.json: inf -> "flaot('inf')"
        sed 's/inf/"float(\'inf\')"/g' metrics/metrics.json | \
        sed 's/nan/"float(\'nan\')"/g' | \
        sed 's/"Description_AVGDP"/"AVGDP: Average depth per sample"/g' | \
        sed 's/"Description_VQSLOD"/"VQSLOD: Log odds of being a true variant versus being false under the trained gaussian mixture model"/g' | \
        sed 's/"Description_SOR"/"SOR: Allele specific strand Odds Ratio of 2x|Alts| contingency table to detect allele specific strand bias"/g' | \
        sed 's/"Description_AS_InbreedingCoeff"/"AS_InbreedingCoeff: Allele-specific inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"/g' | \
        sed 's/"Description_AS_QD"/"AS_QD: Allele-specific variant confidence (Quality by Depth)"/g' | \
        sed 's/"Description_ExcessHet"/"ExcessHet: Phred-scaled p-value for exact test of excess heterozygosity"/g' | \
        sed 's/"Description_FS"/"FS: Phred-scaled p-value using Fisher exact test to detect strand bias"/g' | \
        sed 's/"Description_InbreedingCoeff"/"InbreedingCoeff: Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"/g' | \
        sed 's/"Description_MLEAF"/"MLEAF: Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"/g' | \
        sed 's/"Description_MQ"/"MQ: RMS Mapping Quality"/g' | \
        sed 's/"Description_QD"/"QD: Variant confidence (Quality by Depth)"/g' | \
        sed 's/"Description_RAW_MQ"/"RAW_MQ: Raw data for RMS Mapping Quality (deprecated, use RAW_MQandDP instead)"/g' | \
        sed 's/"Description_BaseQRankSum"/"BaseQRankSum: Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"/g' | \
        sed 's/"Description_MQRankSum"/"MQRankSum: Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"/g' > metrics/metrics-fix.json
    >>>

    output {
        File out = "vcf/~{vcf_basename}.percentiles.vcf.gz"
        File out_index = "vcf/~{vcf_basename}.percentiles.vcf.gz.tbi"
        File metrics_json = "metrics/metrics.json"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "2"
        runtime_minutes: 360
    }
}
