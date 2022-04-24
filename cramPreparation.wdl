workflow prepareCram {
    # VCF file for chromosome
    File chromosomeVCF

    # File for samples
    File samplesFile

    # .txt files listing the location of BAM/CRAMs and their index files
    File sampleLocationFile
    File sampleIndexLocationFile

    call extractId {
        input: chromosomeVCF = chromosomeVCF,
            samplesFile = samplesFile
    }
    call prepareSequences {
        input: chromosomeVCF = chromosomeVCF,
            sampleLocationFile = sampleLocationFile,
            sampleIndexLocationFile = sampleIndexLocationFile
    }
}

task extractId {
    File chromosomeVCF
    File samplesFile
    String sample = basename(chromosomeVCF, ".vcf.gz")

    command {
        RandomHetHom -k 5 -e 1985 -i ${chromosomeVCF} -s ${samplesFile} -o ${sample}.vcf.gz
    }
    output {
        File out = "${sample}.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task prepareSequences {
    File chromosomeVCF
    File sampleLocationFile
    File sampleIndexLocationFile

    # Need to read index and bam files into WDL task scope
    Array[File] sampleFiles = read_lines(sampleLocationFile)
    Array[File] sampleFilesIndex = read_lines(sampleIndexLocationFile)

    command {
        prepare_sequences.py cram -i ${chromosomeVCF} -c ${sep=' ' sampleFiles} -w 100 -o combined.cram
    }
    output {
        File out = "combined.cram"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
