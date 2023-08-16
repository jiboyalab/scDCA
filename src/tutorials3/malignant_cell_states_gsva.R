
###############################
#' TODO: calculate functional states of malignant cells
#' 
#' Author:  Boya Ji && Liwen Xu  date(2023-06)
#'###############################################################################
# R package needs to be installed in advance: org.Hs.eg.db, clusterProfiler, GSVA)
print('############ ------------- TODO: calculate functional states of malignant cells --------------- ############')
print('############ ------------- Author:  Boya Ji && Liwen Xu  date(2023-06) --------------- ############')
print('############ ------------- R package needs to be installed in advance: org.Hs.eg.db, clusterProfiler, GSVA) --------------- ############')
library(optparse)
options(stringsAsFactors = FALSE)


option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count matrix / normalized count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltypes annotation) path"),
  make_option(c("-o", "--output_file_name"), type="character",
              help="the output dir")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);
count_path <- opts$count
meta_path <- opts$meta
output_file_name <- opts$output_file_name


print(paste0('>>> CancerSEA signature gene list ID convert <<< [', Sys.time(),']'))
###
# 1. CancerSEA signature gene list ID convert --------------------
#' 
#'###############################################################################)
source("id_convert.R")
CancerSEA_signature <- get(load("CancerSEA_signature.RData"))
ENSEMBL2SYMBOL <- id_convert(input = unique(unlist(CancerSEA_signature)), to = "SYMBOL", from = "ENSEMBL", database = "org.Hs.eg.db") # 1247 2

CancerSEA_signature_symbol <- lapply(CancerSEA_signature, function(x){
    idx <- match(x, ENSEMBL2SYMBOL$ENSEMBL)
    symbol <- ENSEMBL2SYMBOL$SYMBOL[idx]
    return(symbol)
})
# save(CancerSEA_signature_symbol, file = "./Data/Signature/CancerSEA_signature_symbol.RData")
# CancerSEA_signature_symbol <- get(load(file = "./Data/Signature/CancerSEA_signature_symbol.RData"))
#    Angiogenesis       Apoptosis       CellCycle Differentiation       DNAdamage 
#              73              66             137             201             110 
#       DNArepair             EMT         Hypoxia    Inflammation        Invasion 
#             119              90              83             112              97 
#      Metastasis   Proliferation      Quiescence        Stemness 
#             166              88              66             166 

print(paste0('>>> Extract malignant cell subsets <<< [', Sys.time(),']'))
###
# 2.Extract malignant cell subsets --------------------
#' 
#'###############################################################################)
cell_counts <- read.table(count_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1, check.names = F)
cell_anno <- read.table(meta_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1, check.names = F)

Malignant_cell_counts <- cell_counts[,rownames(cell_anno)[which(cell_anno$cell_type == "Malignant")]] # 32718  2573

# Gene filtering, filtering out genes that are expressed as 0 in all cells
gene_exp_counts <- rowSums(Malignant_cell_counts)
length(which(gene_exp_counts == 0)) # 9641

Malignant_cell_counts <- Malignant_cell_counts[-which(gene_exp_counts == 0),] # 23077  2573


print(paste0('>>> Calculation of malignant cell state activity based on GSVA <<< [', Sys.time(),']'))
###
# 3.Calculation of malignant cell state activity based on GSVA --------------------
#' 
#'###############################################################################)
library(GSVA)
tumor_cell_states_gsva_mat <- gsva(expr = as.matrix(Malignant_cell_counts), 
                gset.idx.list = CancerSEA_signature_symbol, 
                kcdf = "Gaussian" # "Gaussian" for logCPM, logRPKM, logTPM, "Poisson" for counts
                ) # 

write.table(tumor_cell_states_gsva_mat, file = output_file_name, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


