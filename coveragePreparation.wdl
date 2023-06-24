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
        Array[File] extractDepthFiles = extractDepth.outDepth
        Array[File] extractDepthIndices = extractDepth.outIndex
        File aggrBasePair_output = aggrBasePair.outAggrBasePair
        File aggrBasePair_outPruneCov0_25 = aggrBasePair.outPruneCov0_25
        File aggrBasePair_outPruneCov0_50 = aggrBasePair.outPruneCov0_50
        File aggrBasePair_outPruneCov0_75 = aggrBasePair.outPruneCov0_75
        File aggrBasePair_outPruneCov1_00 = aggrBasePair.outPruneCov1_00
        File aggrBasePair_output_index = aggrBasePair.outAggrBasePair_index
        File aggrBasePair_outPruneCov0_25_index = aggrBasePair.outPruneCov0_25_index
        File aggrBasePair_outPruneCov0_50_index = aggrBasePair.outPruneCov0_50_index
        File aggrBasePair_outPruneCov0_75_index = aggrBasePair.outPruneCov0_75_index
        File aggrBasePair_outPruneCov1_00_index = aggrBasePair.outPruneCov1_00_index        
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
        #tar xzf ${referenceFastaCache}
        #export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        #export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"
    
        samtools view -T ${referenceFasta} -q 20 -F 0x0704 -uh ${inputCramFile} ${chromosome} | \
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
        # Filter out any empty files which may break the processing
        cat ${write_lines(inputFiles)} | xargs du | awk '$1>50' | cut -f 2 > files.txt
        create_coverage.py -i files.txt chunk -c ${chromosome} -s 2000000 > commands.list
        bash commands.list
        find $PWD -name "*bgz" > file.list
        mkdir -p full
        merge_coverage.py -i file.list -o full/${chromosome}.full.json.gz
        mkdir -p bin_25e-2
        prune_coverage.py -i full/${chromosome}.full.json.gz -l 0.25 -o bin_25e-2/${chromosome}.bin_0.25.json.gz
        mkdir -p bin_50e-2
        prune_coverage.py -i full/${chromosome}.full.json.gz -l 0.50 -o bin_50e-2/${chromosome}.bin_0.50.json.gz
        mkdir -p bin_75e-2
        prune_coverage.py -i full/${chromosome}.full.json.gz -l 0.75 -o bin_75e-2/${chromosome}.bin_0.75.json.gz
        mkdir -p bin_10e-1
        prune_coverage.py -i full/${chromosome}.full.json.gz -l 1.00 -o bin_10e-1/${chromosome}.bin_1.00.json.gz
    }
    output {
        File outAggrBasePair = "full/${chromosome}.full.json.gz"
        File outAggrBasePair_index = "full/${chromosome}.full.json.gz.tbi"
        File outPruneCov0_25 = "bin_25e-2/${chromosome}.bin_0.25.json.gz"
        File outPruneCov0_25_index = "bin_25e-2/${chromosome}.bin_0.25.json.gz.tbi"
        File outPruneCov0_50 = "bin_50e-2/${chromosome}.bin_0.50.json.gz"
        File outPruneCov0_50_index = "bin_50e-2/${chromosome}.bin_0.50.json.gz.tbi"
        File outPruneCov0_75 = "bin_75e-2/${chromosome}.bin_0.75.json.gz"
        File outPruneCov0_75_index = "bin_75e-2/${chromosome}.bin_0.75.json.gz.tbi"
        File outPruneCov1_00 = "bin_10e-1/${chromosome}.bin_1.00.json.gz"
        File outPruneCov1_00_index = "bin_10e-1/${chromosome}.bin_1.00.json.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
