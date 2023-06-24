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
    parser.add_argument("--cell_type", required=True, help="Specific cell type")
    parser.add_argument("--cell_state", required=True, help="Functional state of malignant cells")
    parser.add_argument("--dca_rank_result", required=True, help="The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern")
    parser.add_argument("--ccc_ratio_result", required=True, help="The result of the ratio of different cell types affected by cellular communication")
    
    args = parser.parse_args()
    
    countfile=args.count
    metafile=args.meta
    LR_pair_database_path=args.lr_file
    gene_name=args.gene
    cell_type=args.cell_type
    cell_state=args.cell_state
    adata=construct_anndata(countfile,metafile) 
    LR_df = load_lr_df(LR_pair_database_path)
    expressed_LR_df = get_expressed_lr_df(LR_df)
    CE_tensor = compute_ce_tensor(adata, lr_df=expressed_LR_df)
    CE_tensor_filtered = filter_ce_tensor(CE_tensor, adata, lr_df=expressed_LR_df,n_pairs=200)
    Sum_CE_tensor_filtered=torch.sum(CE_tensor_filtered, dim=0)
    cell_type_adj,cell_type_pair_df=construct_cell_type_adj(adata,Sum_CE_tensor_filtered,cell_type)
    
    target = torch.zeros((len(adata.obs_names.values), 1)).float()
    tumor_cell_states = pd.read_table('./data/tumor_cell_states_gsva_mat.txt', sep='\t', header=0, index_col=0)
    a=tumor_cell_states.loc[str(cell_state)].index.values
    b=adata.obs_names.values
    indexes = np.where(np.isin(b, a))
    indexes_mask=indexes[0].tolist()
    tumor_cell_states_values=tumor_cell_states.loc[str(cell_state)].values
    tumor_cell_states_values = (tumor_cell_states_values - np.min(tumor_cell_states_values)) / (np.max(tumor_cell_states_values) - np.min(tumor_cell_states_values))
    j=0
    for i in range(len(indexes_mask)):
        target[indexes_mask[i]][0]=tumor_cell_states_values[j]
        j+=1
    
    X, cell_type_names = get_one_hot_cell_type_tensor(adata, categorical_cell_type_col = 'cell_type',)
    torch.manual_seed(1)
    X = torch.randn(adata.n_obs, 130)
    adj = cell_type_adj.to(torch.float32)
    trained_MGC_model_list,device = mgc_repeat_training(X, adj, target, device='cuda:1',repeat_num=50,max_epoch=200)
    predict_result = plot_mgc_result(trained_MGC_model_list, adata, X, adj)
    mse = mean_squared_error(predict_result.T, target.T)
    print("The mean squared error of original and predicted gene expression profiles:",mse)
    print("The Pearson correlation of original and predicted gene expression profiles:",np.corrcoef(predict_result_MMP11.T, target.T)[0,1])
    fname=args.dca_rank_result
    ranked_LR_df_for_MMP11 = lr_rank_in_mgc(trained_MGC_model_list, cell_type_pair_df,repeat_attention_scale=True) 
    ranked_LR_df_for_MMP11.to_csv(fname)
    
       
               
    
