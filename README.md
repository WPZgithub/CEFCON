# CEFCON
This is the PyTorch implementation of the paper:

Peizhuo Wang, Xiao Wen, Peng Lang, Hantao Shu, Lin Gao, Dan Zhao and Jianyang Zeng, 
[Deciphering Driver Regulators of Cell Fate Decisions by Constructing and Controlling Cell-Lineage-Specific Gene Regulatory Networks](https://) 

![Overview.png](https://github.com/WPZgithub/CEFCON/Overview.png)

CEFCON is a computational tool to characterize gene regulatory dynamics from a perspective of network control theory  and identifies putative regulators that
drive cell fate decisions. CEFCON takes a prior gene interaction network and expression profiles from scRNA-seq data as inputs and consists of three main components, including cell-lineage-specific GRN construction, driver regulator
identification and regulon-like gene module identification.

## Installation
### Requirements
- python>=3.7
- numpy, scipy, pandas, scikit-learn, tqdm
- [pytorch>=1.8.0](https://pytorch.org/get-started/locally/) 
- [torch-geometric>=2.1.0](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html)
- anndata>=0.7.0
- networkx>=2.8.0
- cvxpy>=1.2.0
- gurobipy>=9.5.0
- pyscenic>=0.12.0
- Recommended: An NVIDIA GPU with CUDA support for GPU acceleration
### Optional
- matplotlib-venn
- rpy2
### Install using pip
```
pip install git+https://github.com/WPZgithub/CEFCON
```

### Using GUROBI

We recommend using [GRUOBI](https://www.gurobi.com/) to solve the integer linear programming (ILP) problem for identifying driver genes.
GUROBI is a commercial solver that requires licenses to run. Thankfully, it provides free licenses in academia, as well as trial
licenses outside academia. If there is no problem about the licenses, you need to install the
`gurobipy` package.

If you have difficulty using GUROBI, a non-commercial solver, [SCIP](https://www.scipopt.org/), will be used by default.

### Using GPU

We recommend using GPU. If so, you will need to install the GPU version of PyTorch.

## TODO
- [ ] Add support for visualization
- [ ] Add notebook file for analyses

## Input data
The pre-processed data in the paper can be downloaded from [here](https://www.dropbox.com/s/48oe7shjq0ih151/data.tar.gz?dl=0). 
- `scRNA-seq data`: a 'csv' file in which rows represent cells and columns represent genes
- `Prior gene interaction network`: a edgelist formatted network file
```
   From      To      Directed
   Gene1     Gene1   True
   Gene1     Gene2   True
   ...
   GeneX     GeneY   False
```
We provide prior gene interaction networks for human and mouse respectively, located in `/prior_data`.
- `Differential expression level`: 
```
   GeneName    absLogFC
   Gene1       3.25
   Gene2       1.06
   ...
   GeneX       0.71
```

## Usage example
### Command line usage
```
python CEFCON.py [-h] --input_expData PATH --input_priorNet PATH [--input_genesDE PATH] [--TFs PATH] [--additional_edges_pct ADDITIONAL_EDGES_PCT] [--cuda CUDA] [--seed SEED] [--hidden_dim HIDDEN_DIM]
                 [--output_dim OUTPUT_DIM] [--heads HEADS] [--attention {COS,AD,SD}] [--miu MIU] [--epochs EPOCHS] [--repeats REPEATS] [--edge_threshold_param EDGE_THRESHOLD_PARAM] [--remove_self_loops]
                 [--topK_drivers TOPK_DRIVERS] --out_dir OUT_DIR
```
You can run the `run_CEFCON.sh` bash file for an usage example.

- Output (in the output folder `${OUT_DIR}/`):
    - A cell-lineage-specific GRN with default name "cl_GRN.csv"
    - Gene embeddings with default name "gene_embs.csv"
    - A list of driver regulators with default name "driver_regulators.csv"
    - A list of RGMs with default name "RGMs.csv"


## Citation
Please cite the following paper, if you find the repository or the paper useful.

Peizhuo Wang, Xiao Wen, Peng Lang, Hantao Shu, Lin Gao, Dan Zhao and Jianyang Zeng, [Deciphering Driver Regulators of Cell Fate Decisions by Constructing and Controlling Cell-Lineage-Specific Gene Regulatory Networks](https://arxiv.org/abs/2102.07810), Preprint, 2023 

```
@article{wang2022cefcon,
  title={Deciphering Driver Regulators of Cell Fate Decisions by Constructing and Controlling Cell-Lineage-Specific Gene Regulatory Networks},
  author={Wang, peizhuo and Wen, Xiao and Lang, Peng and Shu, Hantao and Gao, Lin and Zhao, Dan and Zeng, Jianyang},
  journal={Preprint},
  year={2023}
}
```

## Bugs & Suggestions
Please contact `wangpeizhuo_37@163.com` or raise an issue in the github repo with any questions.
