library("CMG")

crams = list.files("/mnt/dataSeq/DATA_REPOSITORY/GVCFS_HG38/", recursive = T, pattern = ".cram$", full.names = T)[1:10]
crais = list.files("/mnt/dataSeq/DATA_REPOSITORY/GVCFS_HG38/", recursive = T, pattern = ".crai$", full.names = T)[1:10]

WORKFLOW="prepareCoverage"
WDLParamList <- list()
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "inputCramFiles", INPUT = crams)
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "inputCraiFiles", INPUT = gsub(".cram", ".cram.crai", crams))
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "chromosome", INPUT = "chr1")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFasta", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.fasta")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFastaCache", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.ref_cache.tar.gz")


CromwellOptionList<-list()
CromwellOptionList <- makeDefaultOptionList(CromwellOptionList = CromwellOptionList,
                                            CromwellOutputsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_OUT_DIR"]),
                                            CromwellWfLogDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_WF_LOG_DIR"]),
                                            CromwellCallLogsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_CALL_LOGS_DIR"]),
                                            CromwellServer=CromwellServ,
                                            CromwellHogGroup =  "coverage")

POST_RESPONSE = submitWF(EXOME = "PX0000",
                         CromwellHost = CONFIG_LIST["CROMWELL_HOST"],
                         CromwellPort = CONFIG_LIST[["CROMWELL_PORT_CONIFER"]][CromwellServ],
                         WFUrl = "https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/coveragePreparation.wdl",
                         InputsJsonList = WDLParamList,
                         OptionsJsonList = CromwellOptionList,
                         LabelsJsonList = list(sample=EXOME))
