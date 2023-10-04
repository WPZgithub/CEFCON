#! /usr/bin/env bash

input_exp="./example_data/exp_hESC_hvg1000.csv"
input_prior="./prior_data/network_human.zip"
input_genesDE="./example_data/DE_hESC_hvg1000.csv"

python CEFCON.py --input_expData $input_exp --input_priorNet $input_prior --input_genesDE $input_genesDE --cuda 3 --repeats 3 --edge_threshold_param 8 --out_dir ./out_test
