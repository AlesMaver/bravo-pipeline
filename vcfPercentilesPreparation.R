library("CMG")

WORKFLOW="prepareVCFPercentiles"
WDLParamList <- list()
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "chromosomeVCF", INPUT = "/home/cmg/CMGpipeline_devel/bravo-pipeline/DataPrep/test/input.vcf.gz")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "samplesFile", INPUT = "/home/cmg/CMGpipeline_devel/bravo-pipeline/DataPrep/test/samples.txt")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "assembly", INPUT = "GRCh38")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "bufferSize", INPUT = 100000)
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "lofOptions", INPUT = "")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "lofteeDir", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/vep_data/loftee/")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "referenceDir", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/vep_data/")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFasta", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.fasta")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "cadScores", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/cadd_scores/whole_genome_SNVs.tsv.gz")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "cadScoresIndex", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/cadd_scores/whole_genome_SNVs.tsv.gz.tbi")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "infoFields", INPUT = list("AVGDP"))
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "threads", INPUT = 10)
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "numberPercentiles", INPUT = 10)
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "description", INPUT = "Description")


CromwellOptionList<-list()
CromwellOptionList <- makeDefaultOptionList(CromwellOptionList = CromwellOptionList,
                                            CromwellOutputsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_OUT_DIR"]),
                                            CromwellWfLogDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_WF_LOG_DIR"]),
                                            CromwellCallLogsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_CALL_LOGS_DIR"]),
                                            CromwellServer=CromwellServ,
                                            CromwellHogGroup =  WORKFLOW)

POST_RESPONSE = submitWF(EXOME = "PX0000",
                         CromwellHost = CONFIG_LIST["CROMWELL_HOST"],
                         CromwellPort = CONFIG_LIST[["CROMWELL_PORT_CONIFER"]][CromwellServ],
                         WFUrl = "https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/vcfPercentilesPreparation.wdl",
                         InputsJsonList = WDLParamList,
                         OptionsJsonList = CromwellOptionList,
                         LabelsJsonList = list(sample=EXOME))
