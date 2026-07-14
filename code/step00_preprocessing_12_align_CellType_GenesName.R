library(biomaRt)
library(dplyr)
library(stringr)

## set some parameter
MarkerDataPath = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2021NN_PFC_LayerGene/'
# MarkerSaveName = 'GeneExp_NN.csv'
# MarkerSaveName = 'GeneExp_Mine.csv'
DataName = 'hcp_3d_gxrm.csv'

GEDataPath = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/expression/'
GEDataName = 'hcp_3d_gxrm.csv'
# GEDataPath = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/microarray/'
# GEDataName = 'GeneExpression_glasser360.csv'

## read data
datExpr0 = read.csv(paste(GEDataPath, GEDataName, sep = ""))

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



# 查询基因的ensembl
genes <- datExpr0[,1]

GE_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "external_gene_name",
                        values = genes,
                        mart = mart)
GE_annotations$Index <- match(GE_annotations$external_gene_name, datExpr0$Gene)

write.csv(GE_annotations, paste(MarkerDataPath, MarkerSaveName, sep = ""))






# 
# 
# go_annotations$Index <- match(go_annotations$external_gene_name, datExpr0$Gene)
# go_annotations <- go_annotations[order(go_annotations$go_id), ]
# 
# 
# go_annotations_clean <- go_annotations %>%
#   filter(!is.na(go_id) & go_id != "")
# 
# write.csv(go_annotations_clean, paste(DataPath, "GO_label.csv", sep = ""))
# 
# 
# 
# 
# 
# 
# library(biomaRt)
# library(dplyr)
# 
# # 初始化 biomart
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # 假设你的两个基因向量如下：
# gene_list1 <- datExpr0[,1]
# gene_list2 <- CTdata[,2]
# 
# # 获取注释信息（symbol + Ensembl ID）
# get_gene_info <- function(genes) {
#   getBM(
#     attributes = c("ensembl_gene_id", "external_gene_name", "external_synonym"),
#     filters = "external_gene_name",
#     values = gene_list1,
#     mart = mart
#   ) %>%
#     distinct()
# }
# 
# info1 <- get_gene_info(gene_list1)
# info2 <- get_gene_info(gene_list2)
# 
# # 给来源添加标记
# info1 <- info1 %>% mutate(source = "list1", original_name = external_gene_name)
# info2 <- info2 %>% mutate(source = "list2", original_name = external_gene_name)
# 
# # 合并两个注释表
# combined_info <- bind_rows(info1, info2)
# 
# # 找出相同的 Ensembl ID 或 symbol 或 synonym 来配对
# matched_pairs <- combined_info %>%
#   group_by(ensembl_gene_id) %>%
#   filter(n() > 1) %>%
#   summarise(
#     from_list1 = paste(unique(original_name[source == "list1"]), collapse = "; "),
#     from_list2 = paste(unique(original_name[source == "list2"]), collapse = "; "),
#     .groups = "drop"
#   ) %>%
#   filter(from_list1 != "" & from_list2 != "")
# 
# # 结果展示
# print(matched_pairs)
# 
# 
# 
# 
# 
# # 获取注释信息
# get_gene_info <- function(genes, list_name) {
#   annotations <- getBM(
#     attributes = c("ensembl_gene_id", "external_gene_name", "external_synonym"),
#     filters = "external_gene_name",
#     values = gene_list1,
#     mart = mart
#   ) %>%
#     distinct()
#   
#   # 添加来源和原始位置
#   annotations <- annotations %>%
#     rowwise() %>%
#     mutate(
#       original_name = external_gene_name,
#       index = which(genes == external_gene_name)[1],  # 取第一个匹配的位置
#       source = list_name
#     ) %>%
#     ungroup()
#   
#   return(annotations)
# }
# 
# # 获取两个列表的注释信息
# info1 <- get_gene_info(gene_list1, "list1")
# info2 <- get_gene_info(gene_list2, "list2")
# 
# # 合并注释
# combined_info <- bind_rows(info1, info2)
# 
# # 匹配 Ensembl ID 相同的基因
# matched_pairs <- combined_info %>%
#   group_by(ensembl_gene_id) %>%
#   filter(n() > 1) %>%
#   summarise(
#     from_list1 = paste(unique(original_name[source == "list1"]), collapse = "; "),
#     index_list1 = paste(unique(index[source == "list1"]), collapse = "; "),
#     from_list2 = paste(unique(original_name[source == "list2"]), collapse = "; "),
#     index_list2 = paste(unique(index[source == "list2"]), collapse = "; "),
#     .groups = "drop"
#   ) %>%
#   filter(from_list1 != "" & from_list2 != "")


