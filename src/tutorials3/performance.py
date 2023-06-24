
from typing import List, Optional, Union
import torch
from anndata._core.anndata import AnnData
from model import MGC_Model
from training import get_mgc_result


def plot_mgc_result(trained_MGC_model_list: List[MGC_Model],
                    adata: AnnData,
                    X: torch.Tensor,
                    adj: torch.Tensor,
                    hide_repeat_tqdm: bool = False,
                    **kwarg,
                    ) -> torch.Tensor:
    """

    Plot the generated expression profile of the target gene.

    Parameters
    ----------
    trained_MGC_model_list :
        A list of trained MGC model for generating the expression of one target gene.
    adata :
        Annotated data matrix.
    X :
        The feature matrix used as the input of the trained_MGC_model_list.
    adj :
        The adjacency matrix used as the input of the trained_MGC_model_list.
    hide_repeat_tqdm :
        If true, hide the tqdm for getting the result of repeated training.
    kwarg :
        See in 'feature_plot' function.

    Returns
    -------
    The generated expression profile of the target gene (cell_num * 1).

    """
    result = get_mgc_result(trained_MGC_model_list, X, adj, hide_repeat_tqdm=hide_repeat_tqdm, )
    

    return result



