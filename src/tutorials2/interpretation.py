from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import seaborn as sns
import torch
from matplotlib import pyplot as plt
from torch.nn import functional as F
from tqdm import tqdm

from model import MGC_Model


def lr_rank_in_mgc(trained_MGC_model_list: List[MGC_Model],
                   lr_df: pd.DataFrame,
                   repeat_filter_num: int = 0,
                   repeat_attention_scale: bool = True,
                   ) -> pd.DataFrame:
    """

    Analyze the MGC attention value corresponding to each LR pair, and plot the LR pairs serving as the core mediators
    of FCEs.

    Parameters
    ----------
    trained_MGC_model_list :
        A list of trained MGC model for generating the expression of one target gene.
    lr_df :
        The used preprocessed LR-gene dataframe, must contain the 'LR_pair' column.
    repeat_filter_num :
        The number of repetitions to be filtered out. If the attention obtained from a certain training
        has too low variance, delete the training result. If 0, not filter any training.
    repeat_attention_scale :
        If True, scale the attention of each training to 0-1.
    

    Returns
    -------
    lr_df added column 'MGC_layer_attention', which contains the attention value of each LR pair.

    """
    
    repeat_num = len(trained_MGC_model_list)
    
    layer_attention_final = [trained_MGC_model_list[i].mgc.att.data.view(-1) for i in range(repeat_num)]
    layer_attention_torch = torch.stack(layer_attention_final)

    used_repeat = np.sort(np.argsort(abs(layer_attention_torch).var(1))[repeat_filter_num:])

    df1 = pd.DataFrame(np.array(abs(layer_attention_torch)).T[:, used_repeat],
                       index=lr_df.LR_Pair)
    if repeat_attention_scale:
        df1 = (df1 - df1.min(0)) / (df1.max(0) - df1.min(0))
    
    tmp = pd.DataFrame(df1.mean(1), columns=['weight_sum']).sort_values(by='weight_sum', ascending=False)
    tmp['LR_Pair'] = tmp.index

    related_LR_df_result = lr_df.copy()
    related_LR_df_result.index = lr_df.LR_Pair
    related_LR_df_result['MGC_layer_attention'] = pd.DataFrame(df1.mean(1), columns=['weight_sum']).weight_sum
    related_LR_df_result.index = lr_df.index
    related_LR_df_result = related_LR_df_result.sort_values(by='MGC_layer_attention', ascending=False)


    return related_LR_df_result




def delta_e_proportion(trained_MGC_model_list: List[MGC_Model],
                       X: torch.Tensor,
                       adj: torch.Tensor,
                       cell_type_names: List[str],
                       fname: Optional[Union[str, Path]] = None,
                       device:str="",
                       **kwargs,
                       ) -> pd.DataFrame:
    """

    Plotting the proportion of delta_e in the sum of delta_e and e_0 in each cell-type.

    Parameters
    ----------
    trained_MGC_model_list :
        A list of trained MGC model for generating the expression of one target gene.
    X :
        The feature matrix used as the input of the trained_MGC_model_list.
    adj :
        The adjacency matrix used as the input of the trained_MGC_model_list.
    cell_type_names :
        The list of cell-type names.
    
    fname :
        The output file name. If None, not save the figure.
    
    kwargs :
        Other parameters in seaborn.barplot

    Returns
    -------
    The dataframe containing the proportion of delta_e in the sum of delta_e and e_0 in each cell-type.

    """

    ce_list = []
    b_list = []
      
    for i in tqdm(range(len(trained_MGC_model_list))):
        model = trained_MGC_model_list[i]
        x = model.mgc(X,adj,device)
        x = F.relu(x)

        ce = model.linear_ce(x)
        b = model.linear_b(X)

        ce_list.append(ce)
        b_list.append(b)

    tmp_df = pd.DataFrame()
    b_result = abs(X.T.mm(torch.stack(b_list).mean(0)).sum(1).detach().numpy())
    ce_result = abs(X.T.mm(torch.stack(ce_list).mean(0)).sum(1).detach().numpy())
    # b_result = abs(X.T.mul(torch.stack(b_list).mean(0)).sum(1).detach().numpy())
    # ce_result = abs(X.T.mul(torch.stack(ce_list).mean(0)).sum(1).detach().numpy())
    tmp_df['delta_e'] = ce_result
    tmp_df['e0'] = b_result
    tmp_df['delta_e_proportion'] = ce_result / (b_result + ce_result)
    tmp_df['cell_type'] = cell_type_names
    tmp_df.to_csv(fname)
    

    return tmp_df
