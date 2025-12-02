library(biomaRt)
library(dplyr)

## set some parameter
# DataDir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/'
# FileName = 'microarray/'
DataDir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/'
FileName = 'expression/'

DataPath = paste(DataDir, FileName, sep = "")
# DataName = 'GeneExpression_glasser360.csv'
DataName = 'hcp_3d_gxrm.csv'

## read data
datExpr0 = read.csv(paste(DataPath, DataName, sep = ""))

datExpr0[1:4,1:4]
dim(datExpr0)

# 创建一个转换字典
gene_replace_dict <- c(
  '01-Mar'='MARCH1', '02-Mar'='MARCH2', '03-Mar'='MARCH3', '04-Mar'='MARCH4',
  '05-Mar'='MARCH5', '06-Mar'='MARCH6', '07-Mar'='MARCH7', '08-Mar'='MARCH8', '09-Mar'='MARCH9', '10-Mar'='MARCH10', '11-Mar'='MARCH11',
  '01-Sep'='SEPT1', '02-Sep'='SEPT2', '03-Sep'='SEPT3', '04-Sep'='SEPT4', '05-Sep'='SEPT5', '06-Sep'='SEPT6', '07-Sep'='SEPT7', '08-Sep'='SEPT8',
  '09-Sep'='SEPT9', '10-Sep'='SEPT10', '11-Sep'='SEPT11', '12-Sep'='SEPT12', '13-Sep'='SEPT13', '14-Sep'='SEPT14', '15-Sep'='SEPT15',
  '01-Dec'='DECR1', '02-Dec'='DECR2'
)

# 为防止开头数字没有补0，我们进行统一格式化
datExpr0 <- datExpr0 %>%
  mutate(Gene = if_else(
    str_detect(Gene, "^[0-9]+-Mar|^[0-9]+-Sep|^[0-9]+-Dec"),
    str_replace_all(Gene, "^([0-9]{1})-(...)", "0\\1-\\2"),  # 补0
    Gene
  )) %>%
  mutate(Gene = recode(Gene, !!!gene_replace_dict))



mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 查询基因的 GO 注释
genes <- datExpr0[,1]

go_annotations <- getBM(attributes = c("external_gene_name", "go_id", "name_1006", "namespace_1003"),
                        filters = "external_gene_name",
                        values = genes,
                        mart = mart)

go_annotations$Index <- match(go_annotations$external_gene_name, datExpr0$Gene)
go_annotations <- go_annotations[order(go_annotations$go_id), ]


go_annotations_clean <- go_annotations %>%
  filter(!is.na(go_id) & go_id != "")

write.csv(go_annotations_clean, paste(DataPath, "GO_label.csv", sep = ""))







