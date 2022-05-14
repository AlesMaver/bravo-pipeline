# BRAVO Data Pipeline
Prepare data for [BRAVO](https://github.com/statgen/bravo)

## Setup
This pipeline is written in [WDL](https://software.broadinstitute.org/wdl/) using the Cromwell execution engine.

# Running the BRAVO data preparation pipeline
### Data preparation 
1. Merge multiple vcf.gz files (if Joint Genotyping emits multiple GVCF files rather than a single one)
`bcftools concat cohort1.filtered.{0..9999}.vcf.gz -o ../VCF_MERGED/merged.vcf.gz -Oz`
Index the merged gvcf file
`tabix ../VCF_MERGED/merged.vcf.gz`
(replace "cohort1" with the name of you cohort)

2. Run the data preparation workflow
   - Prepare the inputs json file
     ``` 
     {
     "BravoDataPreparation.input_vcf": "merged.vcf.gz",
     "BravoDataPreparation.input_vcf_index": "merged.vcf.gz.tbi",
     "BravoDataPreparation.samplesFile": "cohort.samples.txt",
     "BravoDataPreparation.sampleLocationFile": "cohort.samples.locations.txt",
     "BravoDataPreparation.referenceFasta": "Homo_sapiens_assembly38.fasta",
     "BravoDataPreparation.cadScores": "/cadd_scores/whole_genome_SNVs.tsv.gz",
     "BravoDataPreparation.cadScoresIndex": "/cadd_scores/whole_genome_SNVs.tsv.gz.tbi",
     "BravoDataPreparation.infoFields": [
       "AVGDP",
       "VQSLOD",
       "SOR",
       "AS_InbreedingCoeff",
       "AS_QD",
       "ExcessHet",
       "FS",
       "InbreedingCoeff",
       "MLEAF",
       "MQ",
       "QD",
       "RAW_MQ",
       "BaseQRankSum",
       "MQRankSum",
       "ClippingRankSum"
        ]
      } 
      ```

      The samplesFile should have the following format:
      ```
      SAMPLE1
      SAMPLE2
      ```

      The sampleLocationFile should have the following format:
      ```
      SAMPLE1 /path_to_crams/SAMPLE1.cram
      SAMPLE2 /path_to_crams/SAMPLE2.cram
      ```

      The CADD scores are available from the following URLs
      ```
      https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
      https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
      ```
