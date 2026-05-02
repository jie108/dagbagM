# DAGBagM: Learning Directed Acyclic Graphs via Bootstrap Aggregation for Mixture of Continuous and Binary Variables

<img src="fig1.png" width="700" align="center">

## Cotents
- [Reference](#Reference)
- [Overview](#Overview)
- [Installation](#Installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Examples](#Examples)
- [Contributions](#contributions)

## Reference

If you use DAGBagM in your research please cite:

Chowdhury, S., Wang, R., Yu, Q. et al. DAGBagM: learning directed acyclic graphs of mixed variables with an application to identify protein biomarkers for treatment response in ovarian cancer. BMC Bioinformatics 23, 321 (2022). https://doi.org/10.1186/s12859-022-04864-y.


## Overview
```
This repository contains two folders. 

dagbagM: 
contains the R package "dagbagM" for learning directed acycic graphs for mixture of continuous and binary variables

Simulation_scripts: 
contains the R scripts for replicating the simulation results in the manuscript. (Notes: dagbag::score_shd is now replaced by dagbagM::score_shd. See Usage and Examples.)
```


## Installation
```
#install.packages("devtools")
library(devtools)
install_github("jie108/dagbagM",subdir="dagbagMv2")
```
or alternatively
```
#install.packages("remotes")
library(remotes)
remotes::install_github("jie108/dagbagM",subdir="dagbagMv2")
```
or install to a local library
```
.libPaths("~/R_libs")
install_github("jie108/dagbagM",subdir="dagbagMv2", lib = "~/R_libs")
```

## Usage

```
hc: Learn a DAG model by hill climbing for mixture of continuous and binary variables.
dagbagMv2::hc(Y, nodeType, whiteList, blackList, tol, standardize, maxStep, restart, seed,
              verbose, debug, addDeleteOnly)

hc_boot: Learn a DAG model for every bootstrap resample (sequential).
dagbagMv2::hc_boot(Y, n.boot, nodeType, whiteList, blackList, standardize, tol, maxStep,
                   restart, seed, nodeShuffle, verbose, debug, addDeleteOnly, return)

hc_boot_parallel: Learn a DAG model for every bootstrap resample (parallel).
dagbagMv2::hc_boot_parallel(Y, n.boot, nodeType, whiteList, blackList, standardize, tol,
                            maxStep, restart, seed, nodeShuffle, numThread, verbose, debug,
                            addDeleteOnly, return)

score_shd: Aggregate an ensemble of DAGs using generalized structural hamming distance.
dagbagMv2::score_shd(boot.adj, alpha, freqCutoff, whiteList, blackList, maxStep, verbose)

score_shd_freq: Aggregate from an edge frequency matrix using generalized structural hamming distance.
dagbagMv2::score_shd_freq(freq, alpha, freqCutoff, whiteList, blackList, maxStep, verbose)
```


## Arguments

### Arguments for dagbagMv2::hc, dagbagMv2::hc_boot, and dagbagMv2::hc_boot_parallel
  
| Parameter | Default | Description |
| :--- | :---: | :--- |
| Y | | an n by p data matrix: n -- sample size, p -- number of variables |
| n.boot (hc_boot, hc_boot_parallel) | 1 | an integer: the number of bootstrap resamples of the data matrix Y |
| nodeType | NULL | a vector of length p specifying the type of each variable: "c" for continuous and "b" for binary |
| whiteList | NULL | a p by p logical or 0-1 matrix: if the (i,j)th-entry is TRUE/1, the edge i-->j will always be included in the DAG |
| blackList | NULL | a p by p logical or 0-1 matrix: if the (i,j)th-entry is TRUE/1, the edge i-->j will be excluded from the DAG |
| standardize | TRUE | logical: whether to standardize continuous nodes to have mean zero and sd one |
| tol | 1e-06 | a scalar: minimum BIC improvement threshold for accepting a step |
| maxStep | 2000 | an integer: the maximum number of search steps of the hill climbing algorithm |
| restart | 1 | an integer: number of random restarts to search for a global optimum |
| seed | 1 | an integer: seed for bootstrap resampling, node shuffling, and tie-breaking |
| nodeShuffle (hc_boot, hc_boot_parallel) | TRUE | logical: whether to shuffle the order of the variables before DAG learning |
| numThread (hc_boot_parallel only) | 2 | an integer: number of parallel workers for bootstrap fitting |
| verbose | FALSE | logical: whether to print step information |
| debug | FALSE | logical: whether to run additional cache and acyclicity checks during HC |
| addDeleteOnly | FALSE | logical: whether to skip edge reversal candidates in HC |
| return (hc_boot, hc_boot_parallel) | "array" | bootstrap output mode: "array" (p x p x B), "freq" (p x p frequencies), or "both" |


### Arguments for dagbagMv2::score_shd and dagbagMv2::score_shd_freq
  
| Parameter | Default | Description |
| :--- | :---: | :--- |
| boot.adj (score_shd only) | | a p by p by B array of bootstrap DAG adjacency matrices |
| freq (score_shd_freq only) | | a p by p edge selection frequency matrix |
| alpha | 1 | a positive scalar: controls the gSHD family member used for aggregation; larger alpha retains fewer edges (lower FDR, less power) |
| freqCutoff | 0.5 | a scalar in [0, 1]: edge selection frequency cutoff; edges with generalized selection frequency above this value are candidates for inclusion. Default 0.5 corresponds to majority-vote selection |
| whiteList | NULL | a p by p logical or 0-1 matrix: if the (i,j)th-entry is TRUE/1, the edge i-->j will always be included |
| blackList | NULL | a p by p logical or 0-1 matrix: if the (i,j)th-entry is TRUE/1, the edge i-->j will be excluded |
| maxStep | NULL | deprecated legacy parameter with no effect |
| verbose | FALSE | logical: whether to print step information |



## Value

### Value for dagbagMv2::hc

a list of five components

| Object | Description |
| :--- | :--- |
| adjacency | adjacency matrix of the learned DAG |
| score | BIC score for each node under the final parent configuration |
| operations | a list of 3-element vectors recording the selected operation (1=add, 2=delete, 3=reverse) at each step |
| deltaMin | score change (delta) at every accepted step |
| steps | total number of accepted steps |

### Value for dagbagMv2::hc_boot and dagbagMv2::hc_boot_parallel

depends on the `return` argument

| return | Value |
| :--- | :--- |
| "array" | a p by p by B array of bootstrap DAG adjacency matrices |
| "freq" | a p by p matrix of edge selection frequencies |
| "both" | a list with components `adjacency` (array) and `freq` (matrix) |


### Value for dagbagMv2::score_shd and dagbagMv2::score_shd_freq

| Object | Description |
| :--- | :--- |
| (matrix) | a p by p integer adjacency matrix of the aggregated DAG |


  
## Examples
```
rm(list=ls())
library(dagbagMv2)
data(example)
Y.n=example$Y # data matrix
p<- dim(Y.n)[2] # no. of nodes: p=102
n<-dim(Y.n)[1] # sample size: n=102

true.dir=example$true.dir  # adjacency matrix of the data generating DAG
true.moral=moral_graph(true.dir) ## moral graph of the data generating DAG
true.ske=skeleton(true.dir)  # skeleton graph of the data generating DAG
true.vstr=vstructures(true.dir) ## vstructures of the data generating DAG

sum(true.dir) #number of edges: |E|= 109

#(i) DAG learning by hill climbing: no bootstrap resample

temp<- dagbagMv2::hc(Y=Y.n,nodeType=rep("c",p), whiteList=NULL, blackList=NULL, tol = 1e-6, standardize=TRUE, maxStep = 1000, restart=1, seed = 1,  verbose = FALSE)
adj.temp=temp$adjacency

# Evaluations
# results on DAG estimation
sum(adj.temp==1&true.dir==0)/sum(adj.temp==1) ## FDR: 0.8803681
sum(adj.temp==1&true.dir==1)/sum(true.dir==1) ## Power: 0.3577982

# results on skeleton graph estimation
adj.temp.ske=skeleton(adj.temp)
sum(adj.temp.ske==1&true.ske==0)/sum(adj.temp.ske==1) ## FDR: 0.7055215
sum(adj.temp.ske==1&true.ske==1)/sum(true.ske==1) ## Power: 0.8807339


#(ii) DAG learning by hill climbing: for bootstrap resamples

library(foreach)
library(doFuture)
boot.adj<- dagbagMv2::hc_boot_parallel(Y=Y.n, n.boot=50, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=1, seed = 1,  nodeShuffle=TRUE, numThread = 2, return="array", verbose = FALSE)

#(iii) Bootstrap aggregation of DAGs learnt from bootstrap resamples: freqCutoff=0.5 corresponds to 50% selection freq. cutoff
adj.bag=dagbagMv2::score_shd(boot.adj, alpha = 1, freqCutoff=0.5) 

#(iv) Evaluations
## results on DAG estimation
sum(adj.bag==1&true.dir==0)/sum(adj.bag==1) ## FDR:  0.3636364
sum(adj.bag==1&true.dir==1)/sum(true.dir==1) ## Power:  0.6422018

## results on skeleton graph estimation
adj.bag.ske=skeleton(adj.bag)
sum(adj.bag.ske==1&true.ske==0)/sum(adj.bag.ske==1) ## FDR: 0.1818182
sum(adj.bag.ske==1&true.ske==1)/sum(true.ske==1) ## Power: 0.8256881

## results on moral graph estimation
adj.bag.moral=moral_graph(adj.bag)
sum(adj.bag.moral==1&true.moral==0)/sum(adj.bag.moral==1) ## FDR: 0.2369942
sum(adj.bag.moral==1&true.moral==1)/sum(true.moral==1) ## Power:0.7173913


## results on vstructures estimation
adj.bag.vstr=vstructures(adj.bag)
vstr.corr=compare.vstructures(target.vstructures=adj.bag.vstr, true.vstructures=true.vstr)
1-nrow(vstr.corr)/nrow(adj.bag.vstr) ## FDR: 0.421875
nrow(vstr.corr)/nrow(true.vstr) ## Power: 0.4805195


```

## Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <jiepeng108@gmail.com>. Contributions (via pull requests or otherwise) are welcome.
