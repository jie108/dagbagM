# Ths folder contains the R scripts used for simulation in the manuscript


## data_generation.R

This file has the scripts for generating data for simulations with mixture of continuous and binary nodes. 

## simulation_continuous: 
This file has the sripts for running dagbagM with and without bootstrap and MMHC and PC-alg for simulations with the continuous nodes only

## simulation_cont+bin: 
This file has the sripts for running dagbagM, dagbagC, bnlearnD and mDAG for simulations with the mixture of continuous and bnary nodes.

## Simulation_results_cont+bin: 
This file has the scripts for summarizing the results obatained using simulation_cont+bin.R scripts.

## Run_time_comparison
This file has the scipt for comparing the run time of one DAG learning (no bootstrap) by hc implemented in DAGBagM with that in bnlearn and mDAG.

## Run_time_comp_adj_matrices
This folder contains all the true adjacency matrices for different number of nodes and edges used for the run time comparison among DAGBagM, bnlearn and mDAG.
