library("CMG")
# Calculate coverage over cohort1
samplesDF<-read.table("/home/ales/JOINT_COHORTS/cohort1/cohort1.samples.txt")
crams<-paste0("/mnt/dataSeq/DATA_REPOSITORY/GVCFS_HG38/", samplesDF$V1,"/", samplesDF$V1, ".cram")

# Loop over all chrs, except for chrM
chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chromosomes = c("chr3", "chr4", "chr5", "chr6")


# Loop across all chromosomes
for (chromosome in chromosomes) {
  print(paste0("Submitting workflow for: ", chromosome))

  WORKFLOW="prepareCoverage"
  WDLParamList <- list()
  WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                              WORKFLOW = WORKFLOW, INPUT_NAME = "inputCramFiles", INPUT = crams)
  WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                              WORKFLOW = WORKFLOW, INPUT_NAME = "inputCraiFiles", INPUT = gsub(".cram", ".cram.crai", crams))
  WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                              WORKFLOW = WORKFLOW, INPUT_NAME = "chromosome", INPUT = chromosome)
  WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                              WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFasta", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.fasta")
  WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                              WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFastaCache", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.ref_cache.tar.gz")
  
  
  CromwellOptionList<-list()
  CromwellOptionList <- makeDefaultOptionList(CromwellOptionList = CromwellOptionList,
                                              CromwellOutputsDir = "/mnt/dataSeq/BRAVO/coverage/",
                                              CromwellWfLogDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_WF_LOG_DIR"]),
                                              CromwellCallLogsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_CALL_LOGS_DIR"]),
                                              CromwellServer=CromwellServ,
                                              CromwellHogGroup =  "coverage_quick")
  
  POST_RESPONSE = submitWF(EXOME = "PX0000",
                           CromwellHost = CONFIG_LIST["CROMWELL_HOST"],
                           CromwellPort = CONFIG_LIST[["CROMWELL_PORT_CONIFER"]][CromwellServ],
                           WFUrl = "http://10.3.248.96:8080/coveragePreparation.wdl",
                           InputsJsonList = WDLParamList,
                           OptionsJsonList = CromwellOptionList,
                           LabelsJsonList = list(sample=EXOME))
}
