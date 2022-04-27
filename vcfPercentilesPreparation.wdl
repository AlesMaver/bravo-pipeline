version 1.0

workflow prepareVCFPercentiles {
    ### Prepare VCF Inputs ###
    input {
        # Chromosome VCF file
        File input_vcf
        File input_vcf_index

        # File for samples
        File samplesFile

        # Size of memory buffer to use
        Int bufferSize  = 100000

        # HG37/H38
        String assembly = "GRCh38"

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
        Int threads = 10
        Int numberPercentiles = 10
        String description = "Description"
    }

    String vcf_basename = basename(input_vcf, ".vcf.gz")

    ###############
    # Prepare VCF #
    ###############

    call FilterVCF {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            referenceFasta = referenceFasta
    }

    call computeAlleleCountsAndHistograms {
        input: chromosomeVCF = FilterVCF.output_vcf,
            samplesFile = samplesFile,
    }
    call variantEffectPredictor {
        input: chromosomeVCF = computeAlleleCountsAndHistograms.out,
            assembly = assembly,
            bufferSize = bufferSize,
            referenceDir = referenceDir,
            referenceFasta = referenceFasta,
            lofteeDir = lofteeDir
    }
    call addCaddScores {
        input: chromosomeVCF = variantEffectPredictor.out,
            cadScores = cadScores,
            cadScoresIndex = cadScoresIndex
    }

    #######################
    # Prepare percentiles #
    #######################

    scatter (field in infoFields) {
        call computePercentiles {
            input: chromosomeVCF = addCaddScores.out,
                infoField = field,
                threads = threads,
                numberPercentiles = numberPercentiles,
                description = description
        }
    }
    call addPercentiles {
        input: 
            chromosomeVCF = addCaddScores.out,
            chromosomeVCFIndex = computePercentiles.outVariantPercentileIndex,
            variantPercentiles = computePercentiles.outVariantPercentile,
            vcf_basename = vcf_basename

    }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = addPercentiles.out
    File output_vcf_index = addPercentiles.out_index
  }

}

# Generate table of variants for interpretation
task FilterVCF {
    input {
      File input_vcf
      File input_vcf_index
      File referenceFasta
    }
  
  command <<<
  set -e
    bcftools +setGT ~{input_vcf} -- -t q -n . -i' FORMAT/GQ<=90' | bcftools norm -m-any -f ~{referenceFasta} | bcftools view --types snps,indels | bcftools +fill-tags -Oz -o output.vcf.gz
    bcftools index -t output.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    maxRetries: 3
    requested_memory_mb_per_core: 5000
    cpu: 1
    runtime_minutes: 90
  }
  output {
    File output_vcf = "output.vcf.gz"
    File output_vcf_index = "output.vcf.gz.tbi"
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
        # docker run -it statgen/bravo-pipeline
        # git clone https://github.com/statgen/bravo_data_prep.git
        # cd bravo_data_prep/data_prep/cpp_tools/
        # cget install .
        docker: "alesmaver/bravo-pipeline-sgp:latest"
    }
}

task  variantEffectPredictor {
    input {
        File chromosomeVCF
        String assembly
        Int bufferSize
        File referenceDir
        File referenceFasta
        File lofteeDir
    }

    command {
        vep -i ~{chromosomeVCF} \
        --plugin LoF,loftee_path:~{lofteeDir} \
        --dir_cache ~{referenceDir} \
        --fasta ~{referenceFasta} \
        --assembly ~{assembly} \
        --cache \
        --offline \
        --vcf \
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
        --shift_hgvs 0 \
        --allele_number \
        --format vcf \
        --force \
        --buffer_size ~{bufferSize} \
        --compress_output bgzip \
        --no_stats \
        --dir_plugins ~{lofteeDir} \
        -o variantEP.vcf.gz
    }
    output {
        File out = "variantEP.vcf.gz"
    }
    runtime {
        #docker: "ensemblorg/ensembl-vep:release_95.1"
        docker: "ensemblorg/ensembl-vep:release_106.1"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}

task addCaddScores {
    input {
        File chromosomeVCF
        File cadScores
        File cadScoresIndex
    }

    command {
        add_cadd_scores.py -i ~{chromosomeVCF} -c ~{cadScores} -o annotated.vcf.gz
    }
    output {
        File out = "annotated.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
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
        ComputePercentiles -i ~{chromosomeVCF} \
        -m ~{infoField} \
        -t ~{threads} \
        -p ~{numberPercentiles} \
        -d ~{description}_~{infoField} \
        -o ~{infoField}

        tabix ~{infoField}.variant_percentile.vcf.gz
    >>>
    output {
        File outAllPercentiles = "~{infoField}.all_percentiles.json.gz"
        File outVariantPercentile = "~{infoField}.variant_percentile.vcf.gz"
        File outVariantPercentileIndex = "~{infoField}.variant_percentile.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: threads
        bootDiskSizeGb: "150"
    }
}

task addPercentiles {
    input {
        File chromosomeVCF
        Array[File] chromosomeVCFIndex
        Array[File] variantPercentiles
        String vcf_basename
    }

    command {
        add_percentiles.py -i ~{chromosomeVCF} -p ~{sep=' ' variantPercentiles} -o ~{vcf_basename}.percentiles.vcf.gz
        tabix ~{vcf_basename}.percentiles.vcf.gz
    }
    output {
        File out = "~{vcf_basename}.percentiles.vcf.gz"
        File out_index = "~{vcf_basename}.percentiles.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}
