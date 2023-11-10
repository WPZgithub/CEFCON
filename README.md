# CEFCON

[![DOI](https://zenodo.org/badge/423767675.svg)](https://zenodo.org/doi/10.5281/zenodo.10101434)

CEFCON is a computational tool for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data.
It takes a prior gene interaction network and expression profiles from scRNA-seq data associated with a given 
developmental trajectory as inputs, and consists of three main components, including cell-lineage-specific gene 
regulatory network (GRN) construction, driver regulator identification and regulon-like gene module (RGM) identification.

![Overview.png](https://github.com/WPZgithub/CEFCON/blob/main/Overview.png)

### About method
CEFCON initially employs the graph attention neural networks under a contrastive learning framework to construct reliable GRNs 
for specific developmental cell lineages (Fig. b). Subsequently, CEFCON characterizes gene regulatory dynamics from the perspective 
of network control theory and identifies the driver regulators that steer cell fate decisions (Fig. c). 
Moreover, CEFCON detects gene regulatory modules (i.e., RGMs) involving the identified driver regulators and measure 
their activities using [AUCell](https://github.com/aertslab/AUCell) (Fig. d). 

## Installation
CEFCON was originally tested on Ubuntu 20.04 with Python (3.8~3.10). 
We recommend running CEFCON on CUDA if possible. 
The following packages are required to be able to run this code:

### Requirements
- python(>=3.8,<=3.10)
- [pytorch(>=1.13.0,<2.0)](https://pytorch.org/get-started/locally/) 
- [torch-geometric(>=2.1.0)](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html)
- [scanpy(>=1.8.2)](https://scanpy.readthedocs.io/en/stable/installation.html)
- networkx(>=3.0)
- cvxpy(>=1.2.0)
- [gurobipy(>=10.0.0)](https://pypi.org/project/gurobipy/)
- [pyscenic(>=0.12.0)](https://pyscenic.readthedocs.io/en/latest/installation.html)
- numpy, scipy, pandas, scikit-learn, tqdm
- Recommended: An NVIDIA GPU with CUDA support for GPU acceleration
### Optional (for performance evaluation, visualization and other analyses)
- matplotlib(>=3.5.3)
- matplotlib-venn(>=0.11.7)
- seaborn(>=0.12.1)
- [palantir(==1.0.1)](https://github.com/dpeerlab/palantir)
- rpy2(>=3.4.1)
- R(>=4.0)
  - PRROC (R package)
  - slingshot (R package)
  - MAST (R package)
### Setup a [conda](https://docs.conda.io/projects/miniconda/en/latest/) environment
```
conda create -y --name CEFCON python=3.10
conda activate CEFCON
```
### Install R and the required packages
```
conda install -y -c conda-forge r
R --no-save -q < ./r_env.R
```
### Install using pip
```
pip install git+https://github.com/WPZgithub/CEFCON.git
```

### Using GUROBI

We recommend using [GRUOBI](https://www.gurobi.com/) to solve the integer linear programming (ILP) problem when identifying driver genes.
GUROBI is a commercial solver that requires licenses to run. Thankfully, it provides free licenses in academia, as well as trial
licenses outside academia. If there is no problem about the licenses, you need to install the
`gurobipy` package.

If difficulties arise while using GUROBI, the non-commercial solver, [SCIP](https://www.scipopt.org/), will be employed as an alternative. But the use of SCIP does not come with a guarantee of achieving a successful solution.

### Using GPU

We recommend using GPU. If you choose to do so, you will need to install the GPU version of PyTorch.

## Usage example
### Command line usage
```bash
cefcon [-h] --input_expData PATH --input_priorNet PATH [--input_genesDE PATH] \
           [--additional_edges_pct ADDITIONAL_EDGES_PCT] [--cuda CUDA] [--seed SEED] \
           [--hidden_dim HIDDEN_DIM] [--output_dim OUTPUT_DIM] [--heads HEADS] [--attention {COS,AD,SD}] \
           [--miu MIU] [--epochs EPOCHS] [--repeats REPEATS] [--edge_threshold_param EDGE_THRESHOLD_PARAM] \
           [--remove_self_loops] [--topK_drivers TOPK_DRIVERS] --out_dir OUT_DIR
```
Please use `cefcon.py -h` to view parameters information. \
Please run the `run_CEFCON.sh` bash file for a usage example.

#### Input data

- `scRNA-seq data`: a '.csv' file in which rows represent cells and columns represent genes, or a '.h5ad' formatted file with AnnData objects.
- `Prior gene interaction network`: an edgelist formatted network file.\
&emsp;We provide prior gene interaction networks for human and mouse respectively, located in `/prior_data`.
- `Gene differential expression level`: a 'csv' file contains the log fold change of each gene.

An example of input data (i.e., the hESC dataset with 1,000 highly variable genes) can be found in `/example_data`.
All the input data mentioned in the paper can be downloaded from [here](https://zenodo.org/record/7564872). 


#### The output results can be found in the folder `${OUT_DIR}/`:
    - "cell_lineage_GRN.csv": the constructed cell-lineage-specific GRN;
    - "gene_embs.csv": the gene embeddings;
    - "driver_regulators.csv": a list of identified driver regulators with their influence scores;
    - "RGMs.csv": a list of obtained RGMs;
    - "AUCell_mtx.csv": the AUCell activity matrix of the obtained RGMs.

### Package usage
**Quick start by an example ([Jupyter Notebook](https://github.com/WPZgithub/CEFCON/blob/main/notebooks/run_CEFCON_nestorowa16_data.ipynb)).** \
**Please check this [Notebook](https://github.com/WPZgithub/CEFCON/blob/main/notebooks/preprocessing_nestorowa16_data.ipynb) for scRNA-seq preprocessing.**
```python
import cefcon as cf

# We assume you have an AnnData object containing scRNA-seq data, cell lineages information,
# and gene differential expression levels (optional).
# We also assume you have a pandas dataframe containing the prior gene interaction network
# in edgelist format.

# Data preparation
data = cf.data_preparation(adata, prior_network)

for lineage, data_li in data.items():
    # Construct cell-lineage-specific GRN
    cefcon_GRN_model = cf.NetModel(epochs=350, repeats=3, cuda='0')
    cefcon_GRN_model.run(data_li)
    
    cefcon_results = cefcon_GRN_model.get_cefcon_results(edge_threshold_avgDegree=8)
    
    # Identify dirver regulators
    cefcon_results.gene_influence_score()
    cefcon_results.driver_regulators()

    # Identify regulon-like gene modules
    cefcon_results.RGM_activity()
```
**Please check this [Notebook](https://github.com/WPZgithub/CEFCON/blob/main/notebooks/run_CEFCON_nestorowa16_data.ipynb) for results visualization and analyses.**


## Citation
Please cite the following paper, if you find the repository or the paper useful.

Peizhuo Wang, Xiao Wen, Han Li, Peng Lang, Shuya Li, Yipin Lei, Hantao Shu, Lin Gao, Dan Zhao and Jianyang Zeng, [A network-based framework for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data](https://github.com/WPZgithub/CEFCON), --, 2023 

```
@article{wang2023cefcon,
  title={A network-based framework for deciphering driver regulators of cell fate decisions from single-cell RNA-seq data},
  author={Wang, peizhuo and Wen, Xiao and Li, Han and Lang, Peng and Li, Shuya and Yipin, Lei and Shu, Hantao and Gao, Lin and Zhao, Dan and Zeng, Jianyang},
  journal={-},
  year={2023}
}
```

## Bugs & Suggestions
Please contact `wangpeizhuo_37@163.com` or raise an issue in the github repo with any questions.
