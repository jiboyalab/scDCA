import numpy as np
import scanpy as sc
import torch
from sklearn.metrics import mean_squared_error
from function import *
from training import *
from interpretation import *
from performance import *
import argparse
if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(usage="it's usage tip.", description="Prioritize the dominant cell communication assmebly that regulates the target gene expression pattern.")
    
    parser.add_argument("--count", required=True, help="Count matrix / normalized count matrix path")
    parser.add_argument("--meta", required=True, help="Meta data (celltypes annotation) path")
    parser.add_argument("--lr_file", required=True, help="LR pair file path")
    parser.add_argument("--gene", required=True, help="Specific target gene name")
    parser.add_argument("--dca_rank_result", required=True, help="The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern")
    parser.add_argument("--ccc_ratio_result", required=True, help="The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern")
    
    args = parser.parse_args()
    
    countfile=args.count
    metafile=args.meta
    LR_pair_database_path=args.lr_file
    gene_name=args.gene
    adata=construct_anndata(countfile,metafile) 
    LR_df = load_lr_df(LR_pair_database_path)
    expressed_LR_df = get_expressed_lr_df(LR_df)
    CE_tensor = compute_ce_tensor(adata, lr_df=expressed_LR_df)
    CE_tensor_filtered = filter_ce_tensor(CE_tensor, adata, lr_df=expressed_LR_df,n_pairs=200)
    Sum_CE_tensor_filtered=torch.sum(CE_tensor_filtered, dim=0)
    cell_type_adj,cell_type_pair_df=construct_cell_type_adj(adata,Sum_CE_tensor_filtered)
    # The target gene should be a highly variable gene
    sc.pp.highly_variable_genes(adata)
    target_all_gene_expr, used_gene_list = get_gene_expr(adata, expressed_LR_df)
    for i in range(len(used_gene_list)):
        if used_gene_list[i]==gene_name:
            target = get_one_case_expr(target_all_gene_expr, cases_list=used_gene_list,used_case_name=used_gene_list[i])
            X, cell_type_names = get_one_hot_cell_type_tensor(adata, categorical_cell_type_col = 'cell_type',)
            adj = cell_type_adj.to(torch.float32)
            trained_MGC_model_list,device = mgc_repeat_training(X, adj, target, device='cuda:1',repeat_num=50,max_epoch=200)
            predict_result = plot_mgc_result(trained_MGC_model_list, adata, X, adj)
            mse = mean_squared_error(predict_result.T, target.T)
            print("The mean squared error of original and predicted gene expression profiles:",mse)
            print("The Pearson correlation of original and predicted gene expression profiles:",np.corrcoef(predict_result_MMP11.T, target.T)[0,1])
            fname=args.dca_rank_result
            ranked_LR_df_for_MMP11 = lr_rank_in_mgc(trained_MGC_model_list, cell_type_pair_df,repeat_attention_scale=True) 
            ranked_LR_df_for_MMP11.to_csv(fname)
            fname=args.ccc_ratio_result
            device="cpu"
            delta_e = delta_e_proportion(trained_MGC_model_list, X, adj,cell_type_names,fname=fname,device=device) 
        else:
            print("Please ensure that the gene is highly variable or adjust the screening parameters.",end="\n")
            exit(0)
               
    
