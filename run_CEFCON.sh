#! /usr/bin/env bash

exp="./data/example/exp_hESC_hvg1000.csv"
prior="./prior_data/network_human.csv.zip"
DE="./data/example/DE_hESC_hvg1000.csv"
TFs="./prior_data/hs_hgnc_tfs_lambert2018.txt"

cefcon --input_expData $exp --input_priorNet $prior --input_genesDE $DE --cuda 3 --repeats 3 --edge_threshold_param 8 --out_dir ./out_test --TFs $TFs
