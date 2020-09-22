# Ths folder contains the R scripts used for simulation in the manuscript


## data_generation.R

This file has the scripts for generating data for simulation. 

## simulation_continuous: 
This file has the sripts for running dagbagM with and without bootstrap and MMHC and PC-alg for simulations with the continuous nodes only

## simulation_cont+bin: 
This file has the sripts for running dagbagM, dagbagC and bnlearnD (results reported in Tables 3 and 4 in Simultion section of the manuscript) simulations with the mixture of continuous and bnary nodes.

## Simulation_results_cont+bin: 
This file has the scripts for summarizing the results obatained using simulation_cont+bin.R scripts.

## Run_time_comparison
This file has the scipt for comparing the run time (Fig. 3) of one DAG learning (no bootstrap) by DAGBagM and bnlearn R package (Scutari, 2009).

## Run_time_comp_adj_matrices
This folder contains all the true adjacency matrices for different number of nods and edges used for the run time comparison between DAGBagM and bnlearn.
