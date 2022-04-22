workflow prepareCoverage {
    # Might be the same files as samples for Cram Prep step?
    Array[File] inputCramFiles

    # Chromosome name e.g. chr22
    String chromosome

    # Reference FASTA file - hg37/hg38
    File referenceFasta

    scatter (file in inputCramFiles) {
        call extractDepth {
            input: inputCramFile = file,
                chromosome = chromosome,
                referenceFasta = referenceFasta
        }
    }
    call aggrBasePair {
        input: inputFiles = extractDepth.outDepth,
            chromosome = chromosome
    }
}

task extractDepth {
    File inputCramFile
    Int chromosome
    File referenceFasta
    String sample = basename(inputCramFile, ".bam")

    command {
        samtools view -q 20 -F 0x0704 -uh ${inputCramFile} ${chromosome} | \
        samtools calmd -uAEr - ${referenceFasta} | \
        bam clipOverlap --in -.ubam --out -.ubam | \
        samtools mpileup -f ${referenceFasta} -Q 20 -t DP - | \
        cut -f1-4 | \
        bgzip > ${chromosome}.${sample}.depth.gz \
        && tabix ${chromosome}.${sample}.depth.gz
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
    Int chromosome
    # Not splitting by BP for now
    Int startBP = 0
    Int endBP = 999999999

    command {
        create_coverage.py -i ${write_lines(inputFiles)} aggregate -c ${chromosome} -s ${startBP} -e ${endBP} | \
        bgzip -c > ${chromosome}.${startBP}.${endBP}.json.gz
    }
    output {
        File outAggrBasePair = "${chromosome}.${startBP}.${endBP}.json.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
