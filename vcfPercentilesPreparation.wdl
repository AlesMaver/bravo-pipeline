workflow prepareVCFPercentiles {
    ### Prepare VCF Inputs ###
    # Chromosome VCF file
    File chromosomeVCF

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

    call computeAlleleCountsAndHistograms {
        input: chromosomeVCF = chromosomeVCF,
            samplesFile = samplesFile,
    }
    call variantEffectPredictor {
        input: chromosomeVCF = computeAlleleCountsAndHistograms.out,
            assembly = assembly,
            lofOptions = lofOptions,
            bufferSize = bufferSize,
            referenceDir = referenceDir,
            referenceFasta = referenceFasta
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

    command {
        vep -i ${chromosomeVCF} \
        --plugin LoF,${lofOptions} \
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
        -o variantEP.vcf.gz
    }
    output {
        File out = "variantEP.vcf.gz"
    }
    runtime {
        #docker: "ensemblorg/ensembl-vep:release_95.1 docker pull ensemblorg/ensembl-vep:release_106.1"
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
