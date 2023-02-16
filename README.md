# CEFCON

CEFCON is a computational tool for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data.
CEFCON takes a prior gene interaction network and expression profiles from scRNA-seq data associated with a given 
developmental trajectory as inputs, and consists of three main components, including cell-lineage-specific gene 
regulatory network (GRN) construction, driver regulator identification and regulon-like gene module (RGM) identification.

![Overview.png](https://github.com/WPZgithub/CEFCON/blob/main/Overview.png)

CEFCON first uses the graph attention neural networks under a contrastive learning framework to construct reliable GRNs 
for individual developmental cell lineages. Then, CEFCON characterizes the gene regulatory dynamics from a perspective 
of network control theory and identifies the driver regulators that steer cell fate decisions. 
CEFCON also detects the gene regulatory modules (i.e., RGMs) involving the identified driver regulators and measure 
their activities based on the [AUCell](https://github.com/aertslab/AUCell) method. 

## Installation
This code was originally run on a Linux x86_64 machine with a GTX3090 NVIDIA GPU.
### Requirements
Please ensure that the following packages are installed in order to run the codes.
- python>=3.8
- [pytorch>=1.8.0](https://pytorch.org/get-started/locally/) 
- [torch-geometric>=2.1.0](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html)
- [scanpy>=1.8.2](https://scanpy.readthedocs.io/en/stable/installation.html)
- networkx>=2.8.0
- cvxpy>=1.2.0
- gurobipy>=9.5.0
- [pyscenic>=0.12.0](https://pyscenic.readthedocs.io/en/latest/installation.html)
- numpy, scipy, pandas, scikit-learn, tqdm
- Recommended: An NVIDIA GPU with CUDA support for GPU acceleration
### Optional (for performance evaluation and other analyses)
- rpy2>=3.4.1
- R>=3.6
  - PRROC (R package)
- matplotlib>=3.5.3
- matplotlib-venn>=0.11.7
- seaborn>=0.12.1
- [palantir==1.0.1](https://github.com/dpeerlab/palantir)
### Install using pip
```
pip install git+https://github.com/WPZgithub/CEFCON.git
```
It may take about 10-20 minutes to install these dependencies.

### Using GUROBI

We recommend using [GRUOBI](https://www.gurobi.com/) to solve the integer linear programming (ILP) problem for identifying driver genes.
GUROBI is a commercial solver that requires licenses to run. Thankfully, it provides free licenses in academia, as well as trial
licenses outside academia. If there is no problem about the licenses, you need to install the
`gurobipy` package.

If you have difficulty using GUROBI, a non-commercial solver, [SCIP](https://www.scipopt.org/), will be used. But it does not ensure a successful solution.

### Using GPU

We recommend using GPU. If so, you will need to install the GPU version of PyTorch.

## Input data

- `Prior gene interaction network`: an edgelist formatted network file.\
&emsp;We provide prior gene interaction networks for human and mouse respectively, located in `/prior_data`.
- `scRNA-seq data`: a '.csv' file in which rows represent cells and columns represent genes, or a '.h5ad' formatted file with AnnData objects. 
- `Differential expression level`: a 'csv' file contains the log fold change of each gene.

An example of input data (i.e., the hESC dataset with 1,000 highly variable genes) are located in `/example_data`
All the input data in the paper can be downloaded from [here](https://zenodo.org/record/7564872). 

## Usage example
### Command line usage
```
cefcon [-h] --input_expData PATH --input_priorNet PATH [--input_genesDE PATH] [--TFs PATH] \
           [--additional_edges_pct ADDITIONAL_EDGES_PCT] [--cuda CUDA] [--seed SEED] \
           [--hidden_dim HIDDEN_DIM] [--output_dim OUTPUT_DIM] [--heads HEADS] [--attention {COS,AD,SD}] \
           [--miu MIU] [--epochs EPOCHS] [--repeats REPEATS] [--edge_threshold_param EDGE_THRESHOLD_PARAM] \
           [--remove_self_loops] [--topK_drivers TOPK_DRIVERS] --out_dir OUT_DIR
```
Please use `cefcon.py -h` to view parameters information. \
You can run the `run_CEFCON.sh` bash file for a usage example.

- The output results can be found in the folder `${OUT_DIR}/`:
    - "cl_GRN.csv": the constructed cell-lineage-specific GRN;
    - "gene_embs.csv": the obtained gene embeddings;
    - "driver_regulators.csv": a list of identified driver regulators;
    - "RGMs.csv": a list of obtained RGMs;
    - "AUCell_mtx.csv": the AUCell activity matrix of the obtained RGMs.

It may take about 2-5 minutes to run on the example data.

## Citation
Please cite the following paper, if you find the repository or the paper useful.

Peizhuo Wang, Xiao Wen, Peng Lang, Han Li, Hantao Shu, Lin Gao, Dan Zhao and Jianyang Zeng, [A network-based framework for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data](https://github.com/WPZgithub/CEFCON), Preprint, 2023 

```
@article{wang2023cefcon,
  title={A network-based framework for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data},
  author={Wang, peizhuo and Wen, Xiao and Lang, Peng and Li, Han and Shu, Hantao and Gao, Lin and Zhao, Dan and Zeng, Jianyang},
  journal={Preprint},
  year={2023}
}
```

## Bugs & Suggestions
Please contact `wangpeizhuo_37@163.com` or raise an issue in the github repo with any questions.
