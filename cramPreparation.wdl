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
    String sampleLocationPath
    #File sampleLocationFile
    #File sampleIndexLocationFile

    Int threads = 40
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
            sampleLocationPath = sampleLocationPath,
            #sampleLocationFile = sampleLocationFile,
            #sampleIndexLocationFile = sampleIndexLocationFile,
            referenceFasta = referenceFasta,
            threads = threads
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
        docker: "alesmaver/bravo-pipeline-sgp:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
        runtime_minutes: 10
    }
}

task prepareSequences {
    input {
        String chromosome
        File input_vcf
        File chromosomeVCF
        File chromosomeVCFIndex
        String sampleLocationPath
        #File sampleLocationFile
        #File sampleIndexLocationFile
        File referenceFasta
        Int threads
        }

    String chromosome_filename = sub(sub(chromosome, "-", "_"), ":", "__")

    # Need to read index and bam files into WDL task scope
    # Array[File] sampleFiles = read_lines(sampleLocationFile)
    # Array[File] sampleFilesIndex = read_lines(sampleIndexLocationFile)

    command <<<
        non_comment_lines=$(zcat "~{chromosomeVCF}" | grep -v '^#' | wc -l)
        
        mkdir -p cram
        if [ "$non_comment_lines" -ne 0 ]; then
            echo "OK, we have variants that need to be processed..."

            bcftools query -l ~{input_vcf} | awk '{print $0, "~{sampleLocationPath}/"$0".cram", "~{sampleLocationPath}/"$0".cram.crai"}' OFS="\t" > samples_locations.txt
            python3 /srv/data/bravo_data_prep/data_prep/py_tools/prepare_sequences2.py cram -i ~{chromosomeVCF} -c samples_locations.txt -w 100 -r ~{referenceFasta} -o ~{chromosome_filename}.cram
            samtools index ~{chromosome_filename}.cram
            samtools sort -@ ~{threads} ~{chromosome_filename}.cram -o ~{chromosome_filename}.bam
            samtools index ~{chromosome_filename}.bam
            rm ~{chromosome_filename}.cram*
            samtools view -T ~{referenceFasta} -O CRAM ~{chromosome_filename}.bam -o cram/~{chromosome_filename}.cram
            samtools index cram/~{chromosome_filename}.cram
            rm ~{chromosome_filename}.bam*
        else
            echo "Sorry no variants to process, creating empty cram and crai to allow continuation of the WF..."

            touch cram/~{chromosome_filename}.cram
            touch cram/~{chromosome_filename}.cram.crai
        fi
    >>>
    output {
        File combined_cram = "cram/~{chromosome_filename}.cram"
        File combined_cram_index = "cram/~{chromosome_filename}.cram.crai"
    }
    runtime {
        docker: "alesmaver/bravo-pipeline-sgp:latest"
        #cpu: "4"
        cpu: threads
        bootDiskSizeGb: "50"
        #continueOnReturnCode: true
        runtime_minutes: 720
    }
}
