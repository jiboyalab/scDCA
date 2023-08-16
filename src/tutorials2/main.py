import numpy as np
import scanpy as sc
import torch
from sklearn.metrics import mean_squared_error
from function import *
from training import *
from interpretation import *
from performance import *
import argparse,time
import warnings

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser(usage="it's usage tip.", description="Prioritize the dominant cell communication assmebly that regulates the target gene expression pattern.")
    parser.add_argument("--count", required=True, help="Count matrix / normalized count matrix path")
    parser.add_argument("--meta", required=True, help="Meta data (celltypes annotation) path")
    parser.add_argument("--lr_file", required=True, help="LR pair file path")
    parser.add_argument("--gene", required=True, help="Specific target gene name")
    parser.add_argument("--cell_type", required=True, help="Specific cell type")
    parser.add_argument("--facked_LR", required=True, help="Faked ligand and receptor genes",default=200)
    parser.add_argument("--device", required=True, help="The device for model training",default="cpu")
    parser.add_argument("--repeat_num", required=True, help="The repeat number for model training",default=50)
    parser.add_argument("--max_epoch", required=True, help="The max epoch for model training",default=200)
    parser.add_argument("--learning_rate", required=True, help="The learning rate for model training",default=0.01)
    parser.add_argument("--display_loss", required=True, help="display training loss for model training",default=True)
    parser.add_argument("--dca_rank_result", required=True, help="The result filename of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern")
    parser.add_argument("--ccc_ratio_result", required=True, help="The result filename of the ratio of different cell types affected by cellular communication")
    args = parser.parse_args()
    countfile=args.count
    metafile=args.meta
    LR_pair_database_path=args.lr_file
    gene_name=args.gene
    cell_type=args.cell_type
    facked_LR=int(args.facked_LR)
    device=args.device
    repeat_num=int(args.repeat_num)
    max_epoch=int(args.max_epoch)
    learning_rate=float(args.learning_rate)
    display_loss=args.display_loss

    print('############ ------------- scDecipher (the key factors in specific cell type)--------------- ############')
    print('>>> arguments <<< ', args)
    print('>>> loading library and data <<< ', time.ctime())
    print('>>> construct an AnnData object from a count file and a metadata file <<< ', time.ctime())
    adata=construct_anndata(countfile,metafile)
    print('>>> load the provided dataframe with the information on ligands and receptors <<< ', time.ctime())
    LR_df = load_lr_df(LR_pair_database_path)
    expressed_LR_df = get_expressed_lr_df(LR_df)
    print('>>> calculate the CE tensor by considering the expression levels of ligand and receptor genes <<< ', time.ctime())
    CE_tensor = compute_ce_tensor(adata, lr_df=expressed_LR_df)
    print('>>> filter the edge in calculated CE tensor, removing the edges with low specificities <<< ', time.ctime())
    CE_tensor_filtered = filter_ce_tensor(CE_tensor, adata, lr_df=expressed_LR_df,n_pairs=facked_LR)
    
    Sum_CE_tensor_filtered=torch.sum(CE_tensor_filtered, dim=0)

    print('>>> construct a cell type adjacency tensor based on the specific cell type and the summed LR-CE tensor. <<< ', time.ctime())
    cell_type_adj,cell_type_pair_df=construct_cell_type_adj(adata,Sum_CE_tensor_filtered,cell_type)

    print('>>> detect the highly variable genes <<< ', time.ctime())
    sc.pp.highly_variable_genes(adata) # The target gene should be a highly variable gene
    target_all_gene_expr, used_gene_list = get_gene_expr(adata, expressed_LR_df)
    
    if gene_name in used_gene_list:
        target = get_one_case_expr(target_all_gene_expr, cases_list=used_gene_list,used_case_name=gene_name)
        X, cell_type_names = get_one_hot_cell_type_tensor(adata, categorical_cell_type_col = 'cell_type')
        adj = cell_type_adj.to(torch.float32)
        print('>>> start training the multi-view graph convolutional neural network <<< ', time.ctime())
        trained_MGC_model_list = mgc_repeat_training(X, adj, target, device=device,repeat_num=repeat_num,max_epoch=max_epoch)
        print('>>> calculate the generated expression profile of the target gene. <<< ', time.ctime())
        predict_result = plot_mgc_result(trained_MGC_model_list, adata, X, adj)
        mse = mean_squared_error(predict_result.T, target.T)
        print("The mean squared error of original and predicted gene expression profiles:",mse)
        print("The Pearson correlation of original and predicted gene expression profiles:",np.corrcoef(predict_result.T, target.T)[0,1])
        
        fname=args.dca_rank_result
        print('>>> the dominant cell communication assmebly that regulates the target gene expression pattern is stored at: <<< ',fname, time.ctime())
        ranked_df = dca_rank_in_mgc(trained_MGC_model_list, cell_type_pair_df,repeat_attention_scale=True) 
        ranked_df.to_csv(fname)

        fname=args.ccc_ratio_result
        print('>>> the ratio of different cell types affected by cellular communication is stored at: <<< ',fname, time.ctime())
        delta_e = delta_e_proportion(trained_MGC_model_list, X, adj,cell_type_names) 
        delta_e.to_csv(fname)
        
    else:
        print("Please ensure that the gene is highly variable or adjust the screening parameters.",end="\n")
        print("We detect the highly variable genes by running sc.pp.highly_variable_genes with default parameters",end="\n")
        print(used_gene_list)
        
            
    
