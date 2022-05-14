# BRAVO Data Pipeline - SGP manual
Prepare data for [BRAVO](https://github.com/statgen/bravo) and is based on workflows in WDL. Data preparation for BRAVO consists of two workflows, outlined below. 

# Workflow 1. Annotate variants, compute metrics and prepare CRAMs
This step will create two sets of data:
1. An annotated VCF file with percentiles, VEP, CADD and metrics, and
2. Precomputed CRAM files for plotting raw data in the browser

Perform the following steps:
1. Merge scattered vcf.gz files and index the created file 
Perform this step if Joint Genotyping emits multiple GVCF files rather than a single one:

`bcftools concat cohort1.filtered.{0..9999}.vcf.gz -o merged.vcf.gz -Oz` 

Replace "cohort1.filtered" with the prefix of your gvcfs in the command above.

2. Index the merged gvcf file:

`tabix merged.vcf.gz`

This step can take several hours. You may consider using GatherVcfs or bcftools if they require speeding up.

3. Prepare the inputs json file as follows:
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
SAMPLE3
...
```

The sampleLocationFile should have the following tab-separated format:
```
SAMPLE1 /path_to_crams/SAMPLE1.cram
SAMPLE2 /path_to_crams/SAMPLE2.cram
SAMPLE2 /path_to_crams/SAMPLE3.cram
...
```

The CADD scores are available from the following URLs
```
https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

4. Run the following workflow with the prepared inputs.json file: `https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/BravoDataPreparation.wdl`


# WORKFLOW 2: Prepare coverage histograms
This step will create JSONs with data that will be plotted in the coverage histograms. 

1. Prepare an inputs json file:

```
{
  "prepareCoverage.inputCramFiles": ["/path_to_crams/SAMPLE1.cram", "/path_to_crams/SAMPLE2.cram", "/path_to_crams/SAMPLE3.cram"],
  "prepareCoverage.inputCraiFiles": ["/path_to_crams/SAMPLE1.cram.crai", "/path_to_crams/SAMPLE2.cram.crai", "/path_to_crams/SAMPLE3.cram.crai"],
  "prepareCoverage.chromosome": "chrY",
  "prepareCoverage.referenceFasta": "Homo_sapiens_assembly38.fasta",
  "prepareCoverage.referenceFastaCache": "Homo_sapiens_assembly38.ref_cache.tar.gz"
} 
```

List paths to input cram and cram.crai file in the `inputCramFiles` and `inputCraiFiles` parameters.
*Note: Initially, do not analyse more than 500 cram files (the workflow scatters quite widely and can take a lot of time with hundreds of samples) - the goal for the coverage presentation is not to include all the samples, but to make an average coverage estimation across representative samples.*

Provide `chromosome` input value with the relevant chromosome from the following list: `["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]`

Get the referenceFastaCache using: `wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz`

2. Run the following workflow for each chromosome: `https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/coveragePreparation.wdl`

**WORKFLOWS 1 and 2 can be run concurrently**

# Expected outputs
The two pipelines will generate the following files:
- An indexed VCF file, containing VEP, CADD, histogram and percentiles information
- A metrics.json file containing calculated metrics for percentile presentation
- One CRAM file per chromosome containing pre-computed reads for IGV.js display
- One JSON file per chromosome containing histogram data with coverage information