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
    parser.add_argument("--cell_type", required=True, help="Specific cell type")
    parser.add_argument("--cell_state_file_path", required=True, help="File path of functional state of malignant cells")
    parser.add_argument("--cell_state", required=True, help="Functional state of malignant cells")
    parser.add_argument("--facked_LR", required=True, help="Faked ligand and receptor genes",default=200)
    parser.add_argument("--device", required=True, help="The device for model training",default="cpu")
    parser.add_argument("--repeat_num", required=True, help="The repeat number for model training",default=50)
    parser.add_argument("--max_epoch", required=True, help="The max epoch for model training",default=200)
    parser.add_argument("--learning_rate", required=True, help="The learning rate for model training",default=0.1)
    parser.add_argument("--display_loss", required=True, help="display training loss for model training",default=True)
    parser.add_argument("--dca_rank_result", required=True, help="The result filename of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern")
    
    args = parser.parse_args()
    countfile=args.count
    metafile=args.meta
    LR_pair_database_path=args.lr_file
    cell_state_file_path=args.cell_state_file_path
    cell_state=args.cell_state
    cell_type=args.cell_type
    facked_LR=int(args.facked_LR)
    device=args.device
    repeat_num=int(args.repeat_num)
    max_epoch=int(args.max_epoch)
    learning_rate=float(args.learning_rate)
    display_loss=args.display_loss

    print('############ ------------- scDecipher (functional states of malignant cells)--------------- ############')
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
    CE_tensor_filtered=CE_tensor
    Sum_CE_tensor_filtered=torch.sum(CE_tensor_filtered, dim=0)

    print('>>> construct a cell type adjacency tensor based on the specific cell type and the summed LR-CE tensor. <<< ', time.ctime())
    cell_type_adj,cell_type_pair_df=construct_cell_type_adj(adata,Sum_CE_tensor_filtered,cell_type)

    print('>>> get the functional states of malignant cells in the dataset <<< ', time.ctime())
    target = torch.zeros((len(adata.obs_names.values), 1)).float()
    tumor_cell_states = pd.read_table(cell_state_file_path, sep='\t', header=0, index_col=0)
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
    print('>>> start training the multi-view graph convolutional neural network <<< ', time.ctime())
    trained_MGC_model_list = mgc_repeat_training(X, adj, target, device=device,repeat_num=repeat_num,max_epoch=max_epoch)
    print('>>> calculate the functional states of malignant cells. <<< ', time.ctime())
    predict_result = plot_mgc_result(trained_MGC_model_list, adata, X, adj)
    mse = mean_squared_error(predict_result.T, target.T)
    print("The mean squared error of original and predicted the functional states of malignant cells:",mse)
    print("The Pearson correlation of original and predicted the functional states of malignant cells:",np.corrcoef(predict_result.T, target.T)[0,1])
    
    fname=args.dca_rank_result
    print('>>> the dominant cell communication assmebly that affected the functional states of malignant cells is stored at: <<< ',fname, time.ctime())
    ranked_df = dca_rank_in_mgc(trained_MGC_model_list, cell_type_pair_df,repeat_attention_scale=True) 
    ranked_df.to_csv(fname)

    
        
    
        
            
    
