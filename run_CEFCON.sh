#! /usr/bin/env bash
echo "Please run this script in the folder of $0"

exp="./data/example/exp_hESC_hvg1000.csv"
prior="./prior_data/network_human.csv"
DE="./data/example/DE_hESC_hvg1000.csv"

echo python CEFCON.py --input_expData $exp --input_priorNet $prior --input_genesDE $DE --cuda 0 --repeats 3 --edge_threshold_param 8 --out_dir ./out_test
