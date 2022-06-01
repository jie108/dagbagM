############### load the following libraries ##############
rm(list=ls())
library(dagbagM)
library(dagbag)
library(doParallel)
library(foreach)
library(bnlearn)
library(mDAG)

########### source the data generating functions ##################
source("data_genertion.R")

### generate data for different n and p ###
n<-100
p<-11 ###### change to 21, 61, 121, 201, 501 and 210 
nrep<-100 ######## number of independent replicates of the data ######## 

data<- data.gen(n=n, p=p)

######## running DAG methods and aggregating bootstrap DAGs using dagbag::score_shd ########

adj.res.m=array(NA, c(p,p,nrep))
adj.res.c=array(NA, c(p,p,nrep))
adj.res.d=array(NA, c(p,p,nrep))
adj.res.md<- array(NA, c(p,p,nrep))

for(rep in 1:nrep){
  
  print(paste0("data_",rep))
  
  load(paste0("data_n",n,"_p",p,"_rep",rep, ".RData"))
  
  ########### dagbagM::cont ############
  
  res.hc.c<- dagbagM::hc_boot_parallel(Y=Y.n, n.boot=100, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)
  
  net.hc.c=score_shd(res.hc.c,threshold=0)$adj.matrix
  
  adj.res.c[,,rep]=net.hc.c
  
  ########### dagbagM::bin ############
  
  res.hc.m<- dagbagM::hc_boot_parallel(Y=Y.n, n.boot=100, nodeType=c(rep("c",p-1),"b"),  whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)
  
  net.hc.m=score_shd(res.hc.m,threshold=0)$adj.matrix
  
  adj.res.m[,,rep]=net.hc.m
  
  ################## bnlearnD ######################
  Y.n.disc<-matrix(0,n,(dim(Y.n)[2]-1))
  for(k in 1:(dim(Y.n)[2]-1))
  {
    temp1<-paste("a_",k,sep="")
    temp2<-paste("b_",k,sep="")
    
    Y.n.disc[,k]<-ifelse((Y.n[,k]>median(Y.n[,k])),temp1,temp2)
  }
  Y.n.disc<-cbind.data.frame(Y.n.disc,Y.n[,p])
  Y.n.disc[,p]<-as.factor(Y.n.disc[,p])
  colnames(Y.n.disc)<- colnames(Y.n)
  
  seed.u<-100
  seed.e<-0
  
  res.hc.d<- array(0,c(p,p,100))
  res.hc.md<- array(0,c(p,p,100))
  for(i in 1:100){
    
    s.pick=sample(1:n, n, replace=TRUE)
    node.rand<-sample(1:p,replace=F) ### randomize node order
    
    Y.n.disc.boot=Y.n.disc[s.pick,node.rand]
    
    s<-sort(node.rand,decreasing=F,ind=T)
    
    Y.n.disc.boot<-as.data.frame(Y.n.disc.boot)
    
    bn.ad<-bnlearn::hc(Y.n.disc.boot,restart = 10, max.iter = 1000)
    res.hc.d[,,i]<-amat(bn.ad)[s$ix,s$ix]
    
    ############## mDAG ##################
    Y.n.boot<- Y.n[s.pick,node.rand]
    
    l<- rep(1,p)
    l[which(colnames(Y.n.boot)=="y")]<- 2
    
    t<- rep("g",p)
    t[which(colnames(Y.n.boot)=="y")]<- "c"
    
    mdag.ad<- mDAG(data=Y.n.boot, type=t, level=l, nperm=10000)
    res.hc.md[,,i]<- mdag.ad$skeleton[s$ix,s$ix]
  }
  
  net.hc.d=score_shd(res.hc.d,threshold=0)$adj.matrix 
  
  adj.res.d[,,rep]=net.hc.d 
  
  net.hc.md=score_shd(res.hc.md,threshold=0)$adj.matrix 
  
  adj.res.md[,,rep]=net.hc.md
}
