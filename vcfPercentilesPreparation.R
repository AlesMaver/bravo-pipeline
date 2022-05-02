library("CMG")

WORKFLOW="BravoDataPreparation"
WDLParamList <- list()
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "input_vcf", INPUT = "/mnt/dataSeq/BRAVO/joint_vcf/merged.vcf.gz")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "input_vcf_index", INPUT = "/mnt/dataSeq/BRAVO/joint_vcf/merged.vcf.gz.tbi")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "samplesFile", INPUT = "/home/ales/JOINT_COHORTS/cohort1/cohort1.samples.txt")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "sampleLocationFile", INPUT = "/home/ales/JOINT_COHORTS/cohort1/cohort1.samples.locations.txt")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList,
                            WORKFLOW = WORKFLOW, INPUT_NAME = "referenceFasta", INPUT = "/mnt/dataSeq/CROMWELL/wgs_reference/Homo_sapiens_assembly38.fasta")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "cadScores", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/cadd_scores/whole_genome_SNVs.tsv.gz")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "cadScoresIndex", INPUT = "/mnt/dataSeq/CROMWELL/vep_reference/cadd_scores/whole_genome_SNVs.tsv.gz.tbi")
WDLParamList <- addWDLInput(WDLParamList = WDLParamList, 
                            WORKFLOW = WORKFLOW, INPUT_NAME = "infoFields", INPUT = list("AVGDP", "VQSLOD", "SOR", "AS_InbreedingCoeff", "AS_QD", "ExcessHet", "FS", "InbreedingCoeff", "MLEAF", "MQ", "QD", "RAW_MQ", "BaseQRankSum", "MQRankSum", "ClippingRankSum"))
                                                                                          #AS_InbreedingCoeff AS_QD BaseQRankSum ClippingRankSum DP ExcessHet FS InbreedingCoeff MLEAC MLEAF MQ MQRankSum QD RAW_MQ ReadPosRankSum SOR VQSLOD

CromwellOptionList<-list()
CromwellOptionList <- makeDefaultOptionList(CromwellOptionList = CromwellOptionList,
                                            CromwellOutputsDir = "/mnt/dataSeq/BRAVO/cromwell_outputs/",
                                            CromwellWfLogDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_WF_LOG_DIR"]),
                                            CromwellCallLogsDir = paste0(CONFIG_LIST[["REFERENCE_DIR"]][CromwellServ], CONFIG_LIST["CROMWELL_CALL_LOGS_DIR"]),
                                            CromwellServer=CromwellServ,
                                            CromwellHogGroup = WORKFLOW)

POST_RESPONSE = submitWF(EXOME = "PX0000",
                         CromwellHost = CONFIG_LIST["CROMWELL_HOST"],
                         CromwellPort = CONFIG_LIST[["CROMWELL_PORT_CONIFER"]][CromwellServ],
                         #WFUrl = "https://raw.githubusercontent.com/AlesMaver/bravo-pipeline/master/BravoDataPreparation.wdl",
                         WFUrl = "http://10.3.248.96:8080/BravoDataPreparation.wdl",
                         InputsJsonList = WDLParamList,
                         OptionsJsonList = CromwellOptionList,
                         LabelsJsonList = list(sample=EXOME))
