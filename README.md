# Learning directed acyclic graphs via bootstrap aggregation
- [Overview](#Overview)
- [Installation](#installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Examples](#Examples)
- [Contributions](#contributions)


## Overview

This repository contains 3 folders. 

## dagbag: 
contains the R package "dagbag" for learning directed acycic graphs via bootstrap aggregation for continuous variables

## dagbagM: 
contains the R package "dagbagM" for learning directed acycic graphs via bootstrap aggregation for mixture of continuous and binary variables

## Simulation_scripts: 
contains the R scripts for replicating the simulation results in the manuscript: "Dagbag: an algorithm for learning directed acyclic graphs to 
identify prognostic protein biomarkers in ovarian cancer".



## Installation

```

### Install dagbag
```
install_github("jie108/dagbag/dagbag")
```

### Install dagbagM
```
install_github("jie108/dagbag/dagbagM")
```
```

## Usage
```
dagbag
score: A function to learn a DAG model by the hill climbing algorithm. It can be used to build an ensemble of DAGs (in form of adjacency matrices) based on bootstrap resamples of the data

score(Y, n.boot=0, score.type="BIC", threshold=0, max.step=500,  ini.adj.matrix=NULL, blacklist=NULL, whitelist=NULL, standardize=TRUE,  standardize.boot=TRUE, 
random.forest=FALSE, random.step.length=NULL, nrestart=0, perturb=0, shuffle=FALSE, print=FALSE, EPS=1e-06)

score_shd: A function to use structural hamming distance to aggregate DAGs. It aggregates an ensemble of DAGs to get a DAG that minimizes the overall distance to the ensemble.

score_shd(boot.adj, alpha = 1, threshold=0, max.step = 500, blacklist = NULL, whitelist = NULL, print = FALSE)
```
```
dagbagM

hc: A function to learn a DAG model by the hill climbing algorithm for mixture of continuous and binary variables

hc(Y, node.type, whiteList, blackList, maxStep = 5, verbose = FALSE)
```


## Arguments for score
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| Y	       |	           |an n by p data matrix: n – sample size, p – number of variables
| n.boot         | 0           |an integer: the number of bootstrap resamples of the data matrix Y, default = 0, meaning no bootstrapping
| score.type 	       |	BIC	            |a string: "BIC" or "likelihood"
| threshold  		       | 0	           | a nonnegative scalar: the cutoff value for the change of the score to decide whether to stop the search
| max.step		           | 500             |an integer: the maximum number of search steps of the hill climbing algorithm
| ini.adj.matrix	 	 | NULL           | a p by p 0-1 matrix: the initial graph, default = NULL, meaning the empty graph
| blacklist	         | NULL             | a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will be excluded from the DAG during the search
| whitelist          | NULL           |  a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will always be included in the DAG during the search
| standardize       | TRUE  | logical: whether to standardize the data to have mean zero and sd one
| standardize.boot   | TRUE         | logical: whether to standardize the bootstrap resamples
| random.forest		   | FALSE	    | logical: whether to use the "random forest" idea for further variance reduction
| random.step.length	  |          | a vector: specify “random forest" steps
| nrestart		  | 0    	     | an integer: number of times to restart the search algorithm after a local optimal is achieved. The purpose is to search for global optimal
| perturb		    | 0     	     | an integer: how many random addition/deletion/reversal operations should be used in each random restart
| shuffle		  |   FALSE 	   | logic: whether to shuffle the order of variables before DAG learning. The purpose is to avoid potential systematic biases in simulation studies
| print	  |     FALSE     | logical: whether print the step information
| EPS     |     1e-06     | a scalar: a number to indicate a threshold below which values will be treated as zero


## Arguments for score_shd
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| boot.adj	       |	           | A p by p by B array, where B is the number of DAGs to be aggregated. It records the adjacency matrices. It may be the output of the "score" function.
| alpha         | 1          |a positive scalar: alpha defines which member of the gSHD family should be used to aggregate the DAGs. In general, the larger the alpha, the more aggressive of the aggregation, in that less edges are retained leading to smaller FDR and less power
| threshold 	       |	0	     |a scalar: it defines the frequency cut-off value, "0" corresponds to cut-off 0.5
| max.step		           | 500             |an integer: the maximum number of search steps 
| blacklist	         | NULL             | a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will be excluded from the DAG during the search
| whitelist          | NULL           |  a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will always be included in the DAG during the search
| print		     |     FALSE     | logical: whether print the step information


## Arguments for hc
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| Y	       |	           |an n by p data matrix: n – sample size, p – number of variables
| node.type  		       |         | a vector of length equal to the number of variables specifying the type of variable/node type: "c" for continuous and "b" for binary
| maxStep		           | 500    |an integer: the maximum number of search steps of the hill climbing algorithm
| blacklist	         | NULL    | a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will be excluded from the DAG during the search
| whitelist          | NULL   |  a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will always be included in the DAG during the search
| verbose		     | FALSE   | logical: whether print the step information
| tol     |     1e-06     | a scalar: a number to indicate a threshold below which values will be treated as zero

## Value for score and score_shd

a list of three components

| Object       | Description   |
| :------------------------ | :-------------|
| adj.matrix	  | adjacency matrix of the learned DAG
| final.step    | a number recording how many search steps are conducted before the procedure stops
| movement	    | a matrix recording the selected operation, addition, deletion or reversal of an edge, at each search step

## Value for score and hc

a list of three components

| Object       | Description   |
| :------------------------ | :-------------|
| adjacency	  | adjacency matrix of the learned DAG
| score       | BIC score at each search step
| operations  | a matrix recording the selected operation, addition, deletion or reversal of an edge, at each search step
| deltaMin    |   

## Examples
```
(i) **DAG learning by hill climbing for continuous nodes: no aggregation**

data(example)
Y.n=example$Y # data matrix 
true.dir=example$true.dir  #adjacency matrix of the data generating DAG
true.ske=example$true.ske  # skeleton graph of the data generating DAG

temp=score(Y=Y.n, n.boot=0, score.type="BIC") 
adj=temp$adj.matrix


(ii) DAG learning by bootstrap aggregation for continuous nodes 

set.seed(1)

### generating DAGs for bootstrap resamples

temp.boot=score(Y.n, n.boot=10, score.type="BIC") 
boot.adj=temp.boot$adj.matrix

### aggregating DAGs to lern an ensemble

temp.bag=score_shd(boot.adj, alpha = 1) 
adj.bag=temp.bag$adj.matrix

(iii) DAG learning by bootstrap aggregation for mixture of continuous and binary nodes

temp<- hc(Y, rep("c",p), whiteList, blackList, verbose=F)
adj=temp$adjacency

```


## Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <jiepeng108@gmail.com>. Contributions (via pull requests or otherwise) are welcome.
