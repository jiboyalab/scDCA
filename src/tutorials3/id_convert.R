#' @TODO 基因ID转换
#' @param input 输入需要转换的基因，字符型向量
# > head(input)
# [1] "RRM2"   "HMMR"   "NUSAP1" "UBE2C"  "TOP2A"  "AURKA" 
#' @param from 输入的类型
#' @param to 转换的类型
#' @param database 使用的注释数据库
#' @returnType data.frame
#' @return 转换后的矩阵，两列，第一列为from，第二列为to
# > res
#      SYMBOL ENTREZID
# 1      RRM2     6241
# 2      HMMR     3161
# 3    NUSAP1    51203
# 4     UBE2C    11065
# 5     TOP2A     7153
#
#' @author ZGX
#
id_convert <- function(input=NULL, from="SYMBOL", to="ENTREZID", database="org.Hs.eg.db", homo_gene_file="~/Data/mouse.human.homology.genes.csv"){
  library(clusterProfiler)
  # fromType和toType的取值范围如下
  #   > columns(org.Hs.eg.db)
  #  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  #  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  # [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  # [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
  # [26] "UNIPROT" 
  res <- NA
  if(database=="org.Mm.eg.db"){
    print("Running for mouse genes ...")
    library(org.Mm.eg.db)
    res <- unique(bitr(input, fromType=from, toType=to, OrgDb=database))
    if(!is.null(homo_gene_file)){
      print("Converting homology genes ... The file of homology genes is")
      print(homo_gene_file)
      homo_mouse2human <- read.csv(homo_gene_file)
      homo_mouse2human <- unique(homo_mouse2human[,c("mgi_symbol","hgnc_symbol")])
      res2 <- merge(res, homo_mouse2human, by.x="SYMBOL", by.y="mgi_symbol", all.x=TRUE)
      # res2 <- unique(res2[, c("ENSEMBL", "hgnc_symbol")])
      res2 <- unique(res2[which(!is.na(res2[,"hgnc_symbol"])), ])
      # colnames(res2) <- c("ENSEMBL", "SYMBOL_HOMO")
      return(list(res=res, res_homo=res2))
    }
  }
  if(database=="org.Hs.eg.db"){
    print("Running for human genes ...")
    library(org.Hs.eg.db)
    res <- unique(bitr(input, fromType=from, toType=to, OrgDb=database))
  }
  return(res)
  # 注意，也可以使用 bioMart 对基因进行ID转换，见 https://www.jianshu.com/p/3a0e1e3e41d0
}


exprs_matrix_rowname2gene <- function(input=NULL, gene_name_col=NULL){
  dup_idx <- which(duplicated(input[, gene_name_col]))
  if(length(dup_idx)>0){
    dup_genes <- unique(as.character(input[dup_idx, gene_name_col]))
    rm_idx <- which(input[, gene_name_col] %in% dup_genes)
    input <- input[-rm_idx, ]
  }
  rownames(input) <- input[, gene_name_col]
  input <- input[, setdiff(colnames(input), gene_name_col)]
  return(input)
}

# load("~/Projects/Psoriasis_Mouse_A2A/Res/0.prepare_data/RData/mouse_A2A_exprs_fpkm.RData")
# rownames(mouse_A2A_exprs_fpkm) <- gsub("\\..*", "", rownames(mouse_A2A_exprs_fpkm))
# mouse_A2A_exprs_fpkm <- data.frame(ENSEMBL=rownames(mouse_A2A_exprs_fpkm), mouse_A2A_exprs_fpkm)

# ENSEMBL2SYMBOL_list <- id_convert(input=rownames(mouse_A2A_exprs_fpkm), from="ENSEMBL", to="SYMBOL", database="org.Mm.eg.db")

# mouse_A2A_exprs_fpkm_homo <- merge(ENSEMBL2SYMBOL_list$res_homo, mouse_A2A_exprs_fpkm, by="ENSEMBL")
# mouse_A2A_exprs_fpkm <- merge(ENSEMBL2SYMBOL_list$res, mouse_A2A_exprs_fpkm, by="ENSEMBL")

# mouse_A2A_exprs_fpkm_homo <- exprs_matrix_rowname2gene(input=mouse_A2A_exprs_fpkm_homo, gene_name_col="SYMBOL_HOMO")
# mouse_A2A_exprs_fpkm_homo <- mouse_A2A_exprs_fpkm_homo[, -1]

# mouse_A2A_exprs_fpkm <- exprs_matrix_rowname2gene(input=mouse_A2A_exprs_fpkm, gene_name_col="SYMBOL")
# mouse_A2A_exprs_fpkm <- mouse_A2A_exprs_fpkm[, -1]
# save(mouse_A2A_exprs_fpkm, file="~/Projects/Psoriasis_Mouse_A2A/Res/0.prepare_data/RData/mouse_A2A_exprs_fpkm.RData")
# save(mouse_A2A_exprs_fpkm_homo, file="~/Projects/Psoriasis_Mouse_A2A/Res/0.prepare_data/RData/mouse_A2A_exprs_fpkm_homo.RData")


