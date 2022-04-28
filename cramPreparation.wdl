version 1.0

workflow prepareCram {
    input {

    # VCF file for chromosome
    String chromosome
    File chromosomeVCF

    # File for samples
    File samplesFile

    # Reference file
    File referenceFasta

    # .txt files listing the location of BAM/CRAMs and their index files
    File sampleLocationFile
    #File sampleIndexLocationFile
    }

    call extractId {
        input: chromosomeVCF = chromosomeVCF,
            samplesFile = samplesFile
    }
    call prepareSequences {
        input: 
            chromosome = chromosome,
            chromosomeVCF = chromosomeVCF,
            sampleLocationFile = sampleLocationFile,
            #sampleIndexLocationFile = sampleIndexLocationFile,
            referenceFasta = referenceFasta
    }
    output {
        File combined_cram_result = prepareSequences.combined_cram
        File combined_cram_result_index = prepareSequences.combined_cram_index
    }
}

task extractId {
    input {
        File chromosomeVCF
        File samplesFile
        String sample = basename(chromosomeVCF, ".vcf.gz")
    }
    
    command {
        RandomHetHom -k 5 -e 1985 -i ~{chromosomeVCF} -s ~{samplesFile} -o ~{sample}.vcf.gz
    }
    output {
        File out = "~{sample}.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task prepareSequences {
    input {
        String chromosome
        File chromosomeVCF
        File sampleLocationFile
        #File sampleIndexLocationFile
        File referenceFasta
        }

    # Need to read index and bam files into WDL task scope
    # Array[File] sampleFiles = read_lines(sampleLocationFile)
    # Array[File] sampleFilesIndex = read_lines(sampleIndexLocationFile)

    command {
        python3 /srv/data/bravo_data_prep/data_prep/py_tools/prepare_sequences.py cram -i ~{chromosomeVCF} -c ~{sep=' ' sampleFiles} -w 100 -r ~{referenceFasta} -o ~{chromosome}.cram
        samtools index ~{chromosome}.cram
    }
    output {
        File combined_cram = "combined.cram"
        File combined_cram_index = "combined.cram.crai"
    }
    runtime {
        docker: "alesmaver/bravo-pipeline-sgp:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
