workflow prepareVCFPercentiles {
    ### Prepare VCF Inputs ###
    # Chromosome VCF file
    File input_vcf
    File input_vcf_index

    # File for samples
    File samplesFile

    # Size of memory buffer to use
    Int bufferSize

    # HG37/H38
    String assembly

    # Optional input for VEP lof plugin
    String lofOptions

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
    Int threads
    Int numberPercentiles
    String description

    ###############
    # Prepare VCF #
    ###############

    call FilterVCF {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index
    }

    call computeAlleleCountsAndHistograms {
        input: chromosomeVCF = FilterVCF.output_vcf,
            samplesFile = samplesFile,
    }
    call variantEffectPredictor {
        input: chromosomeVCF = computeAlleleCountsAndHistograms.out,
            assembly = assembly,
            lofOptions = lofOptions,
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
        input: chromosomeVCF = addCaddScores.out,
            chromosomeVCFIndex = computePercentiles.outVariantPercentileIndex,
            variantPercentiles = computePercentiles.outVariantPercentile

    }
}

# Generate table of variants for interpretation
task FilterVCF {
  File input_vcf
  File input_vcf_index
  
  command <<<
  set -e
    bcftools +setGT ~{input_vcfgz} -- -t q -n . -i' FORMAT/GQ<=90' | bcftools +fill-tags -Oz -o output.vcf.gz
    bcftools index -t /home/ales/vcf/output.vcf.gz
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
    File chromosomeVCF
    File samplesFile

    command {
        ComputeAlleleCountsAndHistograms -i ${chromosomeVCF} -s ${samplesFile} -o computeAlleleCtHst.vcf.gz
    }
    output {
        File out = "computeAlleleCtHst.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
    }

}

task  variantEffectPredictor {
    File chromosomeVCF
    String assembly
    String lofOptions
    Int bufferSize
    File referenceDir
    File referenceFasta
    File lofteeDir

    command {
        vep -i ${chromosomeVCF} \
        --plugin LoF,loftee_path:${lofteeDir} \
        --dir_cache ${referenceDir} \
        --fasta ${referenceFasta} \
        --assembly ${assembly} \
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
        --buffer_size ${bufferSize} \
        --compress_output bgzip \
        --no_stats \
        --dir_plugins ${lofteeDir} \
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
    File chromosomeVCF
    File cadScores
    File cadScoresIndex

    command {
        add_cadd_scores.py -i ${chromosomeVCF} -c ${cadScores} -o annotated.vcf.gz
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
    File chromosomeVCF
    String infoField
    Int threads
    Int numberPercentiles
    String description

    command <<<
        ComputePercentiles -i ${chromosomeVCF} \
        -m ${infoField} \
        -t ${threads} \
        -p ${numberPercentiles} \
        -d ${description} \
        -o ${infoField}

        tabix ${infoField}.variant_percentile.vcf.gz
    >>>
    output {
        File outAllPercentiles = "${infoField}.all_percentiles.json.gz"
        File outVariantPercentile = "${infoField}.variant_percentile.vcf.gz"
        File outVariantPercentileIndex = "${infoField}.variant_percentile.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: threads
        bootDiskSizeGb: "150"
    }
}

task addPercentiles {
    File chromosomeVCF
    Array[File] chromosomeVCFIndex
    Array[File] variantPercentiles

    command {
        add_percentiles.py -i ${chromosomeVCF} -p ${sep=' ' variantPercentiles} -o percentiles.vcf.gz
    }
    output {
        File out = "percentiles.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}
