library(GO.db)
library(AnnotationDbi)

# ============================================================
# direct_roots：只保留"化学物质转化"类过程
# 移除：GO:0043170 (macromolecule metabolic，太宽，含转录/翻译)
#        GO:0006139 (nucleobase-containing，含核酸代谢/转录)
#        GO:0006810 (transport，存在摇摆，移至可选)
# ============================================================
direct_roots <- c(
  "GO:0008152", # metabolic process
  "GO:0044238", # primary metabolic process
  "GO:0019748", # secondary metabolic process
  "GO:0009058", # biosynthetic process
  "GO:0009056", # catabolic process
  "GO:0055114", # oxidation-reduction process
  "GO:0044281", # small molecule metabolic process
  "GO:0006629", # lipid metabolic process
  "GO:0005975", # carbohydrate metabolic process
  "GO:0006091", # generation of precursor metabolites and energy
  "GO:0006519", # amino acid metabolic process
  "GO:0006082", # organic acid metabolic process
  "GO:0006732", # coenzyme metabolic process
  "GO:0006766", # vitamin metabolic process
  "GO:0006790", # sulfur compound metabolic process
  "GO:0009308", # amine metabolic process
  "GO:0006810", # transport（加回：代谢物转运是代谢流的一部分）
  "GO:0006869", # lipid transport（保留显式列出）
  "GO:0006865", # amino acid transport
  "GO:0015749", # monosaccharide transmembrane transport
  "GO:0015711", # organic anion transport
  "GO:0006820",  # anion transport
  
  # 新增：明确列出多糖/大分子代谢中属于直接代谢的子类
  "GO:0005976", # polysaccharide metabolic process（多糖代谢）
  "GO:0005977", # glycogen metabolic process（糖原代谢）
  "GO:0006112", # energy reserve metabolic process（能量储备代谢）
  "GO:0044042", # glucan metabolic process（葡聚糖代谢）
  "GO:0006073"  # cellular glucan metabolic process
)

# ============================================================
# exclude_from_all：封锁 direct + indirect（纯信息流过程）
# ============================================================
exclude_from_all <- c(
  "GO:0032774", # RNA biosynthetic process
  "GO:0016070", # RNA metabolic process
  "GO:0006351", # transcription, DNA-templated
  "GO:0006366", # transcription by RNA polymerase II
  "GO:0006396", # RNA processing
  "GO:0010467", # gene expression
  "GO:0006412", # translation
  "GO:0006281", # DNA repair
  "GO:0006260", # DNA replication
  "GO:0006259", # DNA metabolic process
  "GO:0051276", # chromosome organization
  "GO:0006468", # protein phosphorylation（信号修饰）
  "GO:0016567", # protein ubiquitination（信号修饰）
  "GO:0006457", # protein folding
  "GO:0006461", # protein complex assembly
  "GO:0008104", # protein localization
  "GO:0045184", # establishment of protein localization
  "GO:0006886", # intracellular protein transport
  "GO:0006605"  # protein targeting
)

# ============================================================
# exclude_from_direct：只封锁 direct，不影响 indirect 路径
# 用于"本身不是直接代谢，但可能通过祖先归为间接"的情况
# ============================================================
exclude_from_direct <- c(
  "GO:0019538", # protein metabolic process（蛋白质代谢，非小分子代谢）
  "GO:0006996", # organelle organization（结构组织，非化学转化）
  "GO:0022607", # cellular component assembly（结构组装）
  "GO:0070925"  # organelle assembly
)

# ============================================================
# indirect_roots：只保留与代谢有明确机制联系的过程
# 原则：能说清楚"为什么这个过程影响代谢"才保留
# 移除了：signaling、response to stimulus、development、
#          cell cycle、gene expression、biological regulation 等
# ============================================================
indirect_roots <- c(
  # 自噬：直接回收代谢底物
  "GO:0006914", # autophagy
  
  # 囊泡运输：代谢物区室间转运
  "GO:0016192", # vesicle-mediated transport
  
  # 线粒体/过氧化物酶体：能量代谢和脂质氧化的关键区室
  "GO:0007005", # mitochondrion organization
  "GO:0007007", # inner mitochondrial membrane organization
  "GO:0007031", # peroxisome organization
  "GO:0006839", # mitochondrial transport
  
  # 蛋白质降解：释放氨基酸和能量回代谢池
  "GO:0030163", # protein catabolic process
  
  # 细胞死亡：代谢底物释放与回收
  "GO:0006915", # apoptotic process
  "GO:0012501", # programmed cell death
  
  # 代谢直接调控
  "GO:0019222", # regulation of metabolic process
  
  # === 新增：离子/金属辅因子稳态 ===
  # 锌、铁、铜等是数百种代谢酶的必需辅因子，其稳态直接影响代谢
  "GO:0006875", # cellular metal ion homeostasis
  "GO:0030003", # cellular cation homeostasis
  "GO:0006873", # cellular ion homeostasis
  "GO:0006826", # iron ion transport（铁是呼吸链关键辅因子）
  "GO:0006824", # cobalt ion transport（维生素B12前体）
  
  # === 新增：营养/饥饿响应 ===
  # 这类过程是细胞感知代谢状态并调整代谢流的直接机制
  "GO:0009267", # cellular response to starvation（饥饿响应）
  "GO:0031667", # response to nutrient levels（营养水平响应）
  "GO:0014070", # response to organic cyclic compound（代谢物信号）
  "GO:0010033", # response to organic substance（有机物响应）
  
  # === 新增：活性氧与氧化还原稳态 ===
  # ROS 是代谢副产物，其清除直接关联线粒体代谢
  "GO:0000302", # response to reactive oxygen species
  "GO:0045454", # cell redox homeostasis
  
  # === 新增：脂滴/能量储存 ===
  "GO:0019915", # lipid storage（脂质储存）
  "GO:0006001",  # fructose catabolic process（能量储备分解，示例性扩展）
  
  # === 新增：代谢性器官/组织的发育分化 ===
  # 这些过程直接编程代谢酶的表达和代谢表型
  "GO:0045444", # fat cell differentiation（脂肪细胞分化，直接建立脂质代谢表型）
  "GO:0006986", # response to unfolded protein（UPR，ER 代谢应激响应）
  "GO:0006096", # glycolytic process（糖酵解，属于 carbohydrate metabolic 子节点，应已覆盖）
  
  # 线粒体生物发生（比 mitochondrion organization 更明确）
  "GO:0007006", # mitochondrial membrane organization
  "GO:0072655", # establishment of protein localization to mitochondrion
  # （线粒体蛋白输入，代谢酶定位的关键步骤）
  
  # 脂滴生物发生（脂质代谢区室）
  "GO:0097084"  # lipid droplet organization
)

# ============================================================
# 核心函数
# ============================================================
classify_go_metabolism <- function(go_id) {
  
  if (!grepl("^GO:\\d{7}$", go_id)) {
    warning(paste("无效的 GO ID 格式:", go_id))
    return(NA_character_)
  }
  
  ancestors <- tryCatch({
    anc <- get(go_id, GOBPANCESTOR)
    c(go_id, anc)
  }, error = function(e) {
    tryCatch({
      anc <- get(go_id, GOMFANCESTOR)
      c(go_id, anc)
    }, error = function(e2) {
      tryCatch({
        anc <- get(go_id, GOCCANCESTOR)
        c(go_id, anc)
      }, error = function(e3) NULL)
    })
  })
  
  if (is.null(ancestors) || length(ancestors) == 0) {
    warning(paste("未找到 GO term:", go_id))
    return(NA_character_)
  }
  
  ancestors <- ancestors[ancestors != "all"]
  ancestors_only <- ancestors[ancestors != go_id]
  
  is_self_direct     <- go_id %in% direct_roots
  is_ancestor_direct <- any(ancestors_only %in% direct_roots)
  
  # direct 路径：被两层 exclude 任一拦截都不算 direct
  is_excluded_direct <- any(ancestors %in% c(exclude_from_all, exclude_from_direct))
  
  is_self_indirect     <- go_id %in% indirect_roots
  is_ancestor_indirect <- any(ancestors_only %in% indirect_roots)
  
  # indirect 路径：只被 exclude_from_all 拦截（纯信息流才封死）
  is_excluded_all <- any(ancestors %in% exclude_from_all)
  
  if (is_self_direct) {
    return("direct_metabolism")
  } else if (is_ancestor_direct && !is_excluded_direct) {
    return("direct_metabolism")
  } else if ((is_self_indirect || is_ancestor_indirect) && !is_excluded_all) {
    return("indirect_metabolism")
  } else {
    return("non_metabolism")
  }
}

# 向量化批量版本
classify_go_metabolism_batch <- function(go_ids) {
  results <- vapply(go_ids, classify_go_metabolism, character(1))
  data.frame(go_id = go_ids, category = results, stringsAsFactors = FALSE)
}



# ============================================================
# real
# ============================================================
goid_ccsi = read.csv('/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data/group/GO_sigCCSI_goid.csv')
out <- classify_go_metabolism_batch(goid_ccsi$go_id)
write.csv(out, '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data/group/GO_sigCCSI_goid_wCate.csv')





# ============================================================
# 使用示例
# ============================================================

# 单个查询
classify_go_metabolism("GO:0006006")  # glucose metabolic process → direct_metabolism
classify_go_metabolism("GO:0006915")  # apoptotic process         → indirect_metabolism
classify_go_metabolism("GO:0007049")  # cell cycle                → indirect_metabolism
classify_go_metabolism("GO:0060429")  # epithelium development    → non_metabolism

# 批量查询
test_ids <- c(
'GO:0015701',
'GO:0006869',
'GO:0051592',
'GO:0007616',
'GO:0050673',
'GO:0051289',
'GO:0015914',
'GO:0007163',
'GO:0000045',
'GO:0001662',
'GO:0006897',
'GO:0042789',
'GO:0051056',
'GO:0007173',
'GO:0001568',
'GO:0033209',
'GO:0034198',
'GO:0045893',
'GO:0070534',
'GO:0120163',
'GO:0030155',
'GO:1902459',
'GO:0046488',
'GO:0030326',
'GO:0034446',
'GO:0006508',
'GO:0030282',
'GO:1904262',
'GO:0009117',
'GO:0042733',
'GO:0006915',
'GO:0051966',
'GO:0006661',
'GO:0005977',
'GO:0032092',
'GO:0065003',
'GO:0010508',
'GO:1990830',
'GO:0030148',
'GO:0045944',
'GO:0030216',
'GO:0071805',
'GO:0007517',
'GO:0045669',
'GO:0006821',
'GO:0010507',
'GO:0043124',
'GO:0043491',
'GO:0048538',
'GO:1902476',
'GO:0032956',
'GO:0032508',
'GO:0007059',
'GO:0008360',
'GO:0009791',
'GO:0000122',
'GO:0031398',
'GO:0006665',
'GO:0030163',
'GO:0007005',
'GO:0030514',
'GO:0031333',
'GO:0060291'
)
out <- classify_go_metabolism_batch(test_ids)







