version 1.0

workflow prepareCram {
    input {

    # VCF file for chromosome
    String chromosome
    File chromosomeVCF
    File chromosomeVCFIndex

    # File for samples
    File samplesFile

    # Reference file
    File referenceFasta

    # .txt files listing the location of BAM/CRAMs and their index files
    File sampleLocationFile
    #File sampleIndexLocationFile
    }

    call extractId {
        input: 
            chromosomeVCF = chromosomeVCF,
            chromosomeVCFIndex = chromosomeVCFIndex,
            samplesFile = samplesFile
    }
    call prepareSequences {
        input: 
            chromosome = chromosome,
            input_vcf = chromosomeVCF,
            chromosomeVCF = extractId.out,
            chromosomeVCFIndex = extractId.out_index,
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
        File chromosomeVCFIndex
        File samplesFile
        String sample = basename(chromosomeVCF, ".vcf.gz")
    }
    
    command {
        bcftools query -l ~{chromosomeVCF} > samples.txt
        /srv/data/bravo_data_prep/data_prep/cpp_tools/cget/bin/RandomHetHom -k 5 -e 1234 -i ~{chromosomeVCF} -s samples.txt -o ~{sample}.vcf.gz
        tabix ~{sample}.vcf.gz
    }
    output {
        File out = "~{sample}.vcf.gz"
        File out_index = "~{sample}.vcf.gz.tbi"
    }
    runtime {
        docker: "alesmaver/bravo-pipeline-sgp"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task prepareSequences {
    input {
        String chromosome
        File input_vcf
        File chromosomeVCF
        File chromosomeVCFIndex
        File sampleLocationFile
        #File sampleIndexLocationFile
        File referenceFasta
        }

    # Need to read index and bam files into WDL task scope
    # Array[File] sampleFiles = read_lines(sampleLocationFile)
    # Array[File] sampleFilesIndex = read_lines(sampleIndexLocationFile)

    command <<<
        bcftools query -l ~{input_vcf} | awk '{print $0, "/mnt/dataSeq/DATA_REPOSITORY/GVCFS_HG38/"$0"/"$0".cram", "/mnt/dataSeq/DATA_REPOSITORY/GVCFS_HG38/"$0"/"$0".cram.crai"}' OFS="\t" > samples_locations.txt
        python3 /srv/data/bravo_data_prep/data_prep/py_tools/prepare_sequences2.py cram -i ~{chromosomeVCF} -c samples_locations.txt -w 100 -r ~{referenceFasta} -o ~{chromosome}.cram
        samtools index ~{chromosome}.cram
        samtools sort ~{chromosome}.cram -o ~{chromosome}.bam
        samtools index ~{chromosome}.bam
        rm ~{chromosome}.cram*
        samtools view -T ~{referenceFasta} -O CRAM ~{chromosome}.bam -o ~{chromosome}.cram
        samtools index ~{chromosome}.cram
        rm ~{chromosome}.bam*
    >>>
    output {
        File combined_cram = "~{chromosome}.cram"
        File combined_cram_index = "~{chromosome}.cram.crai"
    }
    runtime {
        docker: "alesmaver/bravo-pipeline-sgp:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
