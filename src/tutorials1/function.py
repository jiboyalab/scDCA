import matplotlib.pyplot as plt
import csv,os
from anndata._core.anndata import AnnData
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import torch
from tqdm import tqdm
from copy import deepcopy
from scipy.stats import gmean
from typing import Tuple, List, Optional
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import LabelBinarizer,OneHotEncoder


def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return
def ReadMyTsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName),delimiter="\t")
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return
def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return
def set_figure_params(tex_fonts=False, ):
    if tex_fonts:
        # tex configuration
        tex_fonts = {
            "text.usetex": True,
            "font.family": "serif",
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8
        }
        plt.rcParams.update(tex_fonts)
        plt.rcParams['text.latex.preamble'] = [r'\usepackage{underscore}']

def construct_anndata(countfile,metafile):
    """
    This function constructs an AnnData object from a count file and a metadata file.
    
    Parameters:
        countfile (str): The path to the count file. The file should be a tab-separated text file with the first row as gene names and the first column as cell names.
        metafile (str): The path to the metadata file. The file should be a tab-separated text file with one row per cell and one column per metadata field.
    
    Returns:
        adata (AnnData): An AnnData object containing the count data and metadata.
    """
    counts=pd.read_table(countfile,sep="\t",header=0,index_col=0)
    counts_values=np.array(counts.values)
    counts_values=np.transpose(counts_values) # scanpy handles single-cell data with rows as cell names, columns as gene names, and therefore transposes the gene expression matrix
    counts_values=np.array(counts_values,dtype=np.float32)
    counts_sparse = csr_matrix(counts_values)
    adata = ad.AnnData(counts_sparse)
    cell_meta=pd.read_table(metafile,sep="\t")
    adata.obs=cell_meta
    adata.obs_names = counts.columns.values.tolist()  
    adata.var_names = counts.index.values.tolist()
    return adata
    


def load_lr_df(LR_pair_database_path) -> pd.DataFrame:
    """
    Load the provided dataframe with the information on ligands and receptors.
    
    Returns
    -------
    The LR-gene dataframe.
    
    """
    
    if os.path.exists(LR_pair_database_path):
        
        connectomeDB = pd.read_csv(LR_pair_database_path)
    
    used_connectomeDB = connectomeDB.loc[:,['ligand','receptor']]
    used_connectomeDB.columns = ['Ligand_gene_symbol','Receptor_gene_symbol']
    
    return used_connectomeDB

def get_expressed_lr_df(lr_df: pd.DataFrame) -> pd.DataFrame:
    
    expressed_lr_df=lr_df
    expressed_lr_df['LR_Pair'] = [expressed_lr_df.Ligand_gene_symbol[i] + ':' +
                                  expressed_lr_df.Receptor_gene_symbol[i]
                                  for i in range(len(expressed_lr_df))]

    return expressed_lr_df



def compute_ce_tensor(adata: AnnData,lr_df: pd.DataFrame,) -> torch.Tensor:
    
    """
    Compute the Cell-Cell Communication Strength (CE) tensor based on gene expression data.

    Args:
        adata (AnnData): Annotated data containing gene expression information.
        lr_df (pd.DataFrame): DataFrame containing ligand-receptor pairs' gene symbols.

    Returns:
        torch.Tensor: Cell-Cell Communication Strength (CE) tensor.

    Note:
        This function calculates the CE tensor by considering the expression levels of ligand and receptor genes.
        If the data is large, it may require substantial memory for computation.
        We're working on improving this piece of code.
    """
    print("""
    Notes:
        This function calculates the CE tensor by considering the expression levels of ligand and receptor genes.
        If the data is large, it may require substantial memory for computation.
        We're working on improving this piece of code.
    """)
    # Normalize gene expression data
    adata.X=adata.X/adata.X.max()
    # Extract ligand and receptor gene symbols
    expressed_ligand = lr_df.loc[:, 'Ligand_gene_symbol'].tolist()
    expressed_receptor = lr_df.loc[:, 'Receptor_gene_symbol'].tolist()
    # Get gene expression tensor for expressed ligands and receptors
    expressed_ligand_tensor = get_gene_expr_tensor(adata, expressed_ligand)
    expressed_receptor_tensor = get_gene_expr_tensor(adata, expressed_receptor).permute(0, 2, 1)
    # Calculate Cell-Cell Communication Efficacy (CE) tensor
    ce_tensor = expressed_ligand_tensor.mul(expressed_receptor_tensor)
    # Apply a scaling factor to CE tensor
    ce_tensor=ce_tensor / (ce_tensor+0.5)
    return ce_tensor.to(torch.float32)


def dist_factor_calculate(adata: AnnData,
                          w_best: float,
                          obsm_spatial_slot: str = 'spatial',
                          ) -> torch.Tensor:
    """

    Parameters
    ----------
    adata :
        Annotated data matrix.
    w_best :
        A distance parameter in edge weighting function controlling the covering region of ligands.
        'default_w_visium' function provides a recommended value of w_best.
    obsm_spatial_slot :
        The slot name storing the spatial position information of each spot。

    Returns
    -------
    A tensor describing the distance factor.

    """
    position_mat = adata.obsm[obsm_spatial_slot]
    dist_mat = squareform(pdist(position_mat, metric='euclidean'))
    dist_factor = dist_mat / w_best

    dist_factor = np.exp((-1) * dist_factor * dist_factor)
    dist_factor = torch.tensor(dist_factor).to(torch.float32)

    return dist_factor


def distinguish_dist_factor_calculate(adata: AnnData,
                                      lr_df: pd.DataFrame,
                                      w_best: float,
                                      distinguish=False,
                                      ) -> torch.Tensor:
    if distinguish:
        w_best2 = w_best * 2
        dist_factor1 = dist_factor_calculate(adata=adata, w_best=w_best, )
        dist_factor2 = dist_factor_calculate(adata=adata, w_best=w_best2, )
        dist_factor_tensor = dist_factor1.repeat(lr_df.shape[0], 1, 1)
        secreted_index = lr_df[lr_df.Ligand_location == 'Secreted Signaling'].index
        dist_factor_tensor[secreted_index, :, :] = dist_factor2
    else:
        dist_factor_tensor = dist_factor_calculate(adata, w_best=w_best, )

    return dist_factor_tensor


def get_gene_expr_tensor(adata: AnnData,
                         gene_name: List[str],
                         ) -> torch.Tensor:
    """
    Extract gene expression data for a given list of gene names and convert it into a tensor.

    Args:
        adata (AnnData): Annotated data containing gene expression information.
        gene_name (List[str]): List of gene names for which to extract expression data.

    Returns:
        torch.Tensor: Tensor containing gene expression data. The tensor shape is (genes, 1, cells).

    Note:
        This function processes gene expression data from the provided AnnData object.
        For gene names containing '+', it calculates the geometric mean of expression values for the pairs.
        The function then stacks the gene expression data column-wise to create the tensor.
    """
    # gene_expr_mat = adata[:, gene_name].X.toarray().astype(np.float32)
    gene_expr_mat=np.empty(shape=(adata.n_obs,1))# Create an empty array with random values, and remove the first column of values from the back.
    for i in range(len(gene_name)):
        
        if "+" in gene_name[i]:
            pair=[]
            temname=gene_name[i].split("+")
            for i in range(len(temname)):
                data1=adata[:, temname[i]].X.toarray().astype(np.float32).T.tolist()
                pair.append(data1)
            gene_expr_mat=np.column_stack((gene_expr_mat,gmean(pair).T))
        else:
            data=adata[:, gene_name[i]].X.toarray().astype(np.float32)
            gene_expr_mat=np.column_stack((gene_expr_mat,data))
    gene_expr_mat=gene_expr_mat[:,1:]
    gene_expr_tensor = torch.tensor(np.expand_dims(gene_expr_mat, 2)).permute(1, 2, 0)
    
    return gene_expr_tensor
def filter_ce_tensor(ce_tensor: torch.Tensor,
                     adata: AnnData,
                     lr_df: pd.DataFrame,
                     n_pairs: int = 200,
                     thres: float = 0.05,
                     copy: bool = True,
                     ) -> torch.Tensor:
    """
    Filter the edge in calculated CE tensor, removing the edges with low specificities.

    For each LR pair, select faked ligand and receptor genes, which have similar expression levels
    with the ligand and receptor gene in the dataset. Then calculate the background CE tensor using faked LR genes,

    Using permutation tests, require filtered edges with communication event strength larger than
    a proportion of background strengthes.

    Parameters
    ----------
    ce_tensor :
        Calculated CE tensor (LR_pair_num * cell_num * cell_num) by "compute_sc_CE_tensor" function
    adata :
        Annotated data matrix.
    lr_df :
        A preprocessed LR-gene dataframe.
    w_best :
        A distance parameter in edge weighting function controlling the covering region of ligands.
        'default_w_visium' function provides a recommended value of w_best.
    n_pairs :
        The number of faked ligand and receptor genes.
    thres :
        We require filtered edges with communicatin event strength larger than a proportion of background strengthes.
        The parameter is the proportion.
    
    copy :
        If False, change the input ce_tensor and save memory consumption

    Returns
    -------
    A CE tensor which removed the edges with low specificities (LR_pair_num * cell_num * cell_num)

    """
    print("""
    Notes:
        This process will take a long time, if you want to reduce the calculation time, please reduce the facked_LR number, the default value is 200
    """)
    if copy:
        if ce_tensor.is_sparse:
            ce_tensor = deepcopy(ce_tensor.to_dense())
        else:
            ce_tensor = deepcopy(ce_tensor)
    all_genes = [item for item in adata.var_names.tolist()]
               
    means = adata.to_df()[all_genes].mean().sort_values()

    for i in tqdm(range(lr_df.shape[0])):
        

        lr1 = lr_df.Ligand_gene_symbol[i]
        lr2 = lr_df.Receptor_gene_symbol[i]
        
        lr1list=lr1.split("+")
        lr1expression_pair=[]
        for m in range(len(lr1list)):
            i2=means.index.get_loc(lr1list[m])
            lr1expression_temp=means.iloc[i2]
            lr1expression_pair.append(lr1expression_temp)
        lr1expression=gmean(lr1expression_pair)

        lr2list=lr2.split("+")
        lr2expression_pair=[]
        for m in range(len(lr2list)):
            i2=means.index.get_loc(lr2list[m])
            lr2expression_temp=means.iloc[i2]
            lr2expression_pair.append(lr2expression_temp)
        lr2expression=gmean(lr2expression_pair)
        if lr1expression<=lr2expression:
            min=lr1expression
            max=lr2expression
        else:
            min=lr2expression
            max=lr1expression
        pair=[]
        for n in range(len(means)):
            if means.iloc[n]>=min and means.iloc[n]<max:
                pair.append(means.iloc[n])
        
        im=np.argmin(abs(means.values - np.median(pair)))
        if lr1.split("+") == lr2.split("+"):
            selected = (
            abs(means - means.iloc[im]).sort_values().drop(lr1.split("+"))[: n_pairs * 2].index.tolist()
        )
        else:
            selected = (
                abs(means - means.iloc[im]).sort_values().drop(lr1.split("+")).drop(lr2.split("+"))[: n_pairs * 2].index.tolist()
            )
        
        faked_ligand = selected[-n_pairs:]
        faked_receptor = selected[:n_pairs]

        faked_expressed_ligand_tensor = get_gene_expr_tensor(adata, faked_ligand)
        faked_expressed_receptor_tensor = get_gene_expr_tensor(adata, faked_receptor).permute(0, 2, 1)
        
        faked_ce_tensor = faked_expressed_ligand_tensor.mul(faked_expressed_receptor_tensor)
        listlr1=[]
        listlr1.append(lr1)
        listlr2=[]
        listlr2.append(lr2)
        expressed_ligand_tensor = get_gene_expr_tensor(adata, listlr1)
        expressed_receptor_tensor = get_gene_expr_tensor(adata, listlr2).permute(0, 2, 1)
        
        true_ce_tensor = expressed_ligand_tensor.mul(expressed_receptor_tensor)
        tmp = (true_ce_tensor > faked_ce_tensor).sum(0) > n_pairs * (1 - thres)
        ce_tensor[i, :, :] = ce_tensor[i, :, :].mul(tmp)

    return ce_tensor

def construct_cell_type_adj(adata,Sum_LR_CE_tensor):
    """
    Construct a cell type adjacency tensor based on the cell type and the summed LR-CE tensor.

    Args:
        adata (AnnData): Annotated data containing cell type information.
        Sum_LR_CE_tensor (torch.Tensor): Summed LR-CE tensor containing cell-cell communication efficacies.

    Returns:
        torch.Tensor: Cell type adjacency tensor.
        pd.DataFrame: DataFrame with cell type pairs.

    """
    cell_type=np.unique(adata.obs["cell_type"].values)
    columns1=['LR_Pair']
    num=-1
    index1=[]
    cell_type_name=[]
    for i in range(len(cell_type)):
        for j in range(i,len(cell_type)):
            num+=1
            index1.append(num)
            cell_type_name.append(str(cell_type[i])+":"+str(cell_type[j]))
    kong=pd.DataFrame(columns=columns1,index=index1)
    for i in range(len(cell_type_name)):
        kong["LR_Pair"][i]=cell_type_name[i]
    
    cell_type_num=len(cell_type)
    cell_pair_num=int((cell_type_num*(cell_type_num+1))/2)
    cell_type_adj=torch.zeros(cell_pair_num,adata.n_obs,adata.n_obs)
    multi_view_matrix_num=0
    for i in range(cell_type_num):
        print("cell type:",cell_type[i])
        for j in tqdm(range(i,cell_type_num)):
            for n in range(Sum_LR_CE_tensor.shape[0]):
                if adata.obs["cell_type"][n] == cell_type[i]:
                    for m in range(Sum_LR_CE_tensor.shape[1]):
                        if adata.obs["cell_type"][m] == cell_type[j]:
                            cell_type_adj[multi_view_matrix_num][n][m]=float(Sum_LR_CE_tensor[n][m])
                            
                        
                if adata.obs["cell_type"][n] == cell_type[j]:
                    
                    for m in range(Sum_LR_CE_tensor.shape[1]):
                        if adata.obs["cell_type"][m] == cell_type[i]:
                            cell_type_adj[multi_view_matrix_num][n][m]=float(Sum_LR_CE_tensor[n][m])
                            
                        
            multi_view_matrix_num+=1 
    return cell_type_adj,kong

def get_gene_expr(adata: AnnData,
                  lr_df: pd.DataFrame,
                  gene_list_on_interest: Optional[List[str]] = None,
                  min_mean: float = 0.05,
                  min_disp: float = 0.5,
                  max_sparse: float = 0.1,
                  remove_lr_gene: bool = True,
                  ) -> Tuple[torch.Tensor, List[str]]:
    """\
    Filter out the genes with too low expression levels or too low dispersions.
    Filter out the genes expressed in too fewer cells (or spots).
    Filter out the mitochondria genes and ligand (or receptor) genes.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    lr_df
        A pandas dataframe, must contain three columns: 'Ligand_gene_symbol', 'Receptor_gene_symbol' and 'LR_pair'.
    gene_list_on_interest:
        User provided gene list for the predefined genes on interest.
    min_mean:
        The minimum value of the mean expression level of filtered genes. (provided in adata.var.means)
    min_disp:
        The minimum value of the dispersions of filtered gene expressions. (provided in adata.var.dispersions_norm)
    max_sparse:
        The percentage of cells required to express the filtered genes.
    remove_lr_gene
        If True, filter out the ligand (or receptor) genes, to avoid data Leakage in predicting target gene expression
        using ligand–receptor pairs.
    
    Returns
    -------
    target_all_gene_expr
        Generated target gene expression matrix (cell by gene, scaled to make the maximum expression level of
        each target gene become 1)   
    used_gene_list
        The filtered target gene list for the following workflow.
        
    """

    if ('dispersions_norm' not in adata.var.columns) or ('means' not in adata.var.columns):
        raise Exception(
            'Not detect the highly variable genes info in adata.var, please run sc.pp.highly_variable_genes before.')

    hvgenes = adata.var[adata.var.dispersions_norm > min_disp][adata.var.means > min_mean].index
    mt_genes = [adata.var.index[i] for i in range(adata.shape[1]) if adata.var.index[i].split('-')[0] == 'MT']

    non_zero_exp_genes = (np.array(adata.to_df()) != 0).sum(0)
    adata.var['non_zero_exp'] = non_zero_exp_genes / adata.shape[0]
    non_zero_exp_genes = adata.var.index[np.where(non_zero_exp_genes > adata.shape[0] * max_sparse)[0]]

    if gene_list_on_interest is not None:
        used_gene_list = gene_list_on_interest
    else:
        if remove_lr_gene:
            used_gene_list = list(set(non_zero_exp_genes).intersection(set(hvgenes)) - set(mt_genes)
                                  - set(lr_df.Ligand_gene_symbol) - set(lr_df.Receptor_gene_symbol))
        else:
            used_gene_list = list(set(non_zero_exp_genes).intersection(set(hvgenes)) - set(mt_genes))
    
    target_all_gene_expr = adata[:, np.array(used_gene_list)].X.toarray()
    target_all_gene_expr = torch.Tensor(target_all_gene_expr / target_all_gene_expr.max(0)).T

    return target_all_gene_expr, used_gene_list
def get_one_case_expr(target_all_expr: torch.Tensor,
                      cases_list: List[str],
                      used_case_name: str,
                      ) -> torch.Tensor:
    """

    Get a cell_num*1 tensor representing the scaled expression profile of one gene, using as the target of GNN.

    Parameters
    ----------
    target_all_expr :
        Expression matrix of all target genes (from 'get_gene_expr' function).
    cases_list :
        All target gene names (from 'get_gene_expr' function).
    used_case_name :
        The gene name for the output expression vector.

    Returns
    -------
    cell_num*1 tensor representing the scaled expression profile of one gene

    """
    target = target_all_expr[cases_list.index(used_case_name)]
    target = target / target.max()
    target = target.to(torch.float32)
    target = target.unsqueeze(1)

    return target
def find_neighbors(adata):
    position_mat = adata.obsm['spatial']
    dist_mat = squareform(pdist(position_mat, metric='euclidean'))

    boundry = 1.4 * dist_mat[dist_mat > 0].min()
    dist_mat[dist_mat < boundry] = 1
    dist_mat[dist_mat >= boundry] = 0
    dist_mat[range(dist_mat.shape[0]), range(dist_mat.shape[0])] = 0
    return dist_mat
def get_one_hot_cell_type_tensor(adata: AnnData,
                                 categorical_cell_type_col: str = 'predictions',
                                 ) -> Tuple[torch.Tensor, List[str]]:
    """

    Get categorical cell-type labels and one-hot encoded into a matrix, used as the feature matrix of GNN.
    The categorical cell-type labels are stored at one slot in adata.obs.

    Parameters
    ----------
    adata :
        Annotated data matrix with cell-type information.
    categorical_cell_type_col :
        The column name of categorical cell-type labels in in adata.obs.

    Returns
    -------
        cell_type_tensor: cell_num * cell_type_num
        cell_type_names: used cell-type names

    """
    cell_type = list(adata.obs[categorical_cell_type_col])

    one_hot = LabelBinarizer()
    one_hot_cell_type = one_hot.fit_transform(cell_type)

    cell_type_index = list(one_hot.classes_)

    return torch.FloatTensor(one_hot_cell_type), cell_type_index

def train_test_mask(cell_num: int,
                    train_set_ratio: float = 0.6,
                    val_set_ratio: float = 0.2,
                    ) -> Tuple[List[int], List[int], List[int]]:
    """

    Get the index of cells using as the training set, the validation set, or the testing set.
    Set the train_set_ratio and val_set_ratio, the last part is the testing set.

    Parameters
    ----------
    cell_num :
        The number of cell.
    train_set_ratio :
        A value from 0-1. The ratio of cells using as the training set.
    val_set_ratio :
        A value from 0-1. The ratio of cells using as the validation set.

    Returns
    -------
    Three list of cell index, respectively for the training set, the validation set, and the testing set.

    """

    train_cell_num = round(cell_num * train_set_ratio)
    val_cell_num = round(cell_num * val_set_ratio)

    train_mask = np.random.choice(list(range(cell_num)), train_cell_num + val_cell_num, replace=False)
    val_mask = list(np.random.choice(train_mask, val_cell_num, replace=False))

    test_mask = list(set(range(cell_num)) - set(train_mask))
    train_mask = list(set(train_mask) - set(val_mask))

    return train_mask, test_mask, val_mask
def FindLR():
    """
    https://github.com/wanglabtongji/CCI/tree/main/scripts
    #CellChat -> LR  https://www.jianshu.com/p/03f9c2f3ee4f
    Rscript run_cc.R
    
    
    #CellPhoneDB -> LR
    run_cellphonedb.py

   
    #NicheNet -> LR
    
    # 使用这个代码
    # https://github.com/wanglabtongji/CCI
    run_nichenet.R


    #iCellNet
    run_icellnet_refine.R


    接着在每种方法中对lr结果进行处理，运行各个文件夹下的process_lr.py
    """

   
def construct_cell_type(filename,get_filename,meta_file):
    meta=pd.read_table(meta_file,sep='\t',header=0,index_col=None)
    cell_name_list=meta["cell_name"].values.tolist()
    data=pd.read_table(filename,sep='\t',header=0,index_col=0)
    
    index_all=[]
    for index, row in data.iteritems():
        index_all.append(index)
    diff_list=[y for y in index_all if y not in cell_name_list]
    data=data.drop(columns=diff_list)
    data.to_csv(get_filename,sep="\t")
    