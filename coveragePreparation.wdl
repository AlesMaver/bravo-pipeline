workflow prepareCoverage {
    # Might be the same files as samples for Cram Prep step?
    Array[File] inputCramFiles
    Array[File] inputCraiFiles

    # Chromosome name e.g. chr22
    String chromosome

    # Reference FASTA file - hg37/hg38
    File referenceFasta
    # Get reference fasta cache using: wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz
    File referenceFastaCache
    
    scatter (idx in range(length(inputCramFiles))) {
        call extractDepth {
                input: 
                    inputCramFile = inputCramFiles[idx],
                    inputCraiFile = inputCraiFiles[idx],
                    chromosome = chromosome,
                    referenceFasta = referenceFasta,
                    referenceFastaCache = referenceFastaCache
        }
    }
    
    call aggrBasePair {
            input: 
                inputFiles = extractDepth.outDepth,
                inputIndices = extractDepth.outIndex,
                chromosome = chromosome
        }
    
    output {
        File aggrBasePair_output = aggrBasePair.outAggrBasePair
        File aggrBasePair_outPruneCov0_25 = aggrBasePair.outPruneCov0_25
        File aggrBasePair_outPruneCov0_50 = aggrBasePair.outPruneCov0_50
        File aggrBasePair_outPruneCov0_75 = aggrBasePair.outPruneCov0_75
        File aggrBasePair_outPruneCov1_00 = aggrBasePair.outPruneCov1_00
        }
}

task extractDepth {
    File inputCramFile
    File inputCraiFile
    String chromosome
    File referenceFasta
    File referenceFastaCache
    String sample = basename(inputCramFile, ".bam")

    command {
        tar xzf ${referenceFastaCache}
        export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"
    
        samtools view -q 20 -F 0x0704 -uh ${inputCramFile} ${chromosome} | \
        samtools calmd -uAEr - ${referenceFasta} | \
        bam clipOverlap --in -.ubam --out -.ubam | \
        samtools mpileup -f ${referenceFasta} -Q 20 -t DP - | \
        cut -f1-4 | \
        bgzip > ${chromosome}.${sample}.depth.gz \
        && tabix -b 2 ${chromosome}.${sample}.depth.gz
    }
    output {
        File outDepth = "${chromosome}.${sample}.depth.gz"
        File outIndex = "${chromosome}.${sample}.depth.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task aggrBasePair {
    Array[File] inputFiles
    Array[File] inputIndices
    String chromosome
    # Not splitting by BP for now
    Int startBP = 0
    Int endBP = 999999999

    command {
        create_coverage.py -i ${write_lines(inputFiles)} chunk -c ${chromosome} -s 2000000 > commands.list
        bash commands.list
        find $PWD -name "*bgz" > file.list
        merge_coverage.py -i file.list -o ${chromosome}.full.json.gz
        prune_coverage.py -i ${chromosome}.full.json.gz -l 0.25 -o ${chromosome}.bin_0.25.json.gz
        prune_coverage.py -i ${chromosome}.full.json.gz -l 0.50 -o ${chromosome}.bin_0.50.json.gz
        prune_coverage.py -i ${chromosome}.full.json.gz -l 0.75 -o ${chromosome}.bin_0.75.json.gz
        prune_coverage.py -i ${chromosome}.full.json.gz -l 1.00 -o ${chromosome}.bin_1.00.json.gz
    }
    output {
        File outAggrBasePair = "${chromosome}.full.json.gz"
        File outPruneCov0_25 = "${chromosome}.bin_0.25.json.gz"
        File outPruneCov0_50 = "${chromosome}.bin_0.50.json.gz"
        File outPruneCov0_75 = "${chromosome}.bin_0.75.json.gz"
        File outPruneCov1_00 = "${chromosome}.bin_1.00.json.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
