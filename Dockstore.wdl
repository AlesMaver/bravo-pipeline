import "https://raw.githubusercontent.com/statgen/bravo-pipeline/master/vcfPercentilesPreparation.wdl" as vcfPercentilesPreparation
import "https://raw.githubusercontent.com/statgen/bravo-pipeline/master/cramPreparation.wdl" as cramPreparation
import "https://raw.githubusercontent.com/statgen/bravo-pipeline/master/coveragePreparation.wdl" as coveragePreparation

workflow mainWorkflow {
    ### Shared Inputs ###
    # Array of input chromosome VCF files read from .txt file
    File inputChromosomesList
    Array[File] inputChromosomes = read_lines(inputChromosomesList)
    File referenceFasta
    File referenceDir
    File samplesFile

    ### Prepare VCF and Percentiles Inputs ###
    Int bufferSize
    String assembly
    String lofOptions
    File cadScores
    File cadScoresIndex
    Array[String] infoFields
    Int threads
    Int numberPercentiles
    String description

    ### Prepare Cram Inputs ###
    File sampleLocationFile
    File sampleIndexLocationFile

    ### Prepare Coverage Inputs ###
    Array[File] inputCramFiles
    Int chromosomeNum

    scatter (chromosome in inputChromosomes) {
        call vcfPercentilesPreparation.prepareVCFPercentiles {
            input: chromosomeVCF = chromosome,
                samplesFile = samplesFile,
                bufferSize = bufferSize,
                assembly = assembly,
                lofOptions = lofOptions,
                referenceDir = referenceDir,
                referenceFasta = referenceFasta,
                cadScores = cadScores,
                cadScoresIndex = cadScoresIndex,
                infoFields = infoFields,
                threads = threads,
                numberPercentiles = numberPercentiles,
                description = description
        }
        call cramPreparation.prepareCram {
            input: chromosomeVCF = chromosome,
                samplesFile = samplesFile,
                sampleLocationFile = sampleLocationFile,
                sampleIndexLocationFile = sampleIndexLocationFile
        }
        call coveragePreparation.prepareCoverage {
            input: inputCramFiles = inputCramFiles,
                #TODO - get chromosome number from VCF Filename
                chromosome = chromosomeNum,
                referenceFasta = referenceFasta
        }
    }

    meta {
        author: "Jacob Pleiness"
        email: "pleiness@umich.edu"
        description: "Data preperation pipeline for BRAVO https://github.com/statgen/bravo-pipeline"
    }
}
