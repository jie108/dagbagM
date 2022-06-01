library(dagbagM)
library(dagbag)
library(bnlearn)
library(pcalg)
library(doParallel)
library(foreach)
############ functions to obtain v-structures and skeleton with number of total and correct edges #############
vstructures<-function(adj.matrix){
  ##para: adj.matrix: adjacency matrix
  ##return: 3-column matrix: (par1, child , par2): note, par 1< par 2
  
  p=ncol(adj.matrix)
  res=NULL
  for(i in 1:p){
    parent.node=which(adj.matrix[,i]!=0)
    n.par=length(parent.node)
    
    if(n.par>1){
      for(j in 1:(n.par-1)){
        for(k in (j+1):n.par){
          res=rbind(res,c(parent.node[j],i,parent.node[k]))
        }
      }##for loop
    }##end if
  }## end i loop
  
  if(!is.null(res)){
    colnames(res)=c("par1","child","par2")
  }
  return(res)
}

result_skeleton<-function(adj.m, network.skel){
  ##para: adj.m: adjacency matrix, true.ske -- true skeleton (note: symmetric)
  
  diag(adj.m)=0
  tt=adj.m+t(adj.m) 
  correct.c=sum((tt>0)&(true.ske>0))/2    
  total.c=sum(tt>0)/2
  return(c(total.c, correct.c))
}

######### obtain the correct and total number of edges #############
##############

n<-250
k=500
p<- 1000 ####### Simulation I
p<- 504 ######## Simulation II

##############################################
true.adj<- matrix(0,p,p) ######## for simulation I

load("true_adj_p_504.RData")
true.adj<- adjn ############## for simulation II

network.skel=(true.adj+t(true.adj))>0+0

true.ske=network.skel
sum(true.ske)
######################################################
############### data generation ###############
nrep<- 100

for(rep in 1:nrep){

   Y.n<- NULL
  for(i in 1:p){
    temp<-rnorm(n,0,1)
    temp.s<- (temp-mean(temp))/sd(temp) 

  Y.n<- cbind(Y.n, temp.s)
  }
   file.name=paste0("data_n250_","p_1000_","rep",rep,".RData")
   save(Y.n,file=file.name)
}

########## running DAGBagM hc, hc with no boot, mmhc and pcalg ##############
adj.res.c<- array(NA, c(p,p,nrep))
adj.res.c.noboot<- array(NA, c(p,p,nrep))
adj.res.c.mmhc<- array(NA, c(p,p,nrep))
adj.res.c.pcalg<- array(NA, c(p,p,nrep))

dagbagm<- NULL
hc<- NULL
mmhc<- NULL
pcalg<- NULL
for(rep in 1:100)
{
  load(paste0("data_n250_","p_1000_","rep",rep,".RData"))
  
  print(paste0("running ",rep))

  ########### running dagbagM with bootstrap ################  
  res.hc.c<- dagbagM::hc_boot_parallel(Y=Y.n, n.boot=100, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)

 net.hc.c=dagbag::score_shd(res.hc.c,threshold=0)$adj.matrix
 
 adj.res.c[,,rep]=net.hc.c
 
 ############## running dagbagM without bootstrap ##############
 
 net.hc.c.noboot<- dagbagM::hc(Y=as.matrix(Y.n), nodeType=rep("c",p), blackList=NULL, whiteList=NULL, standardize=TRUE, restart=10,maxStep=1000, verbose=FALSE, tol=1e-06)
 
 adj.res.c.noboot[,,rep]<- net.hc.c.noboot$adjacency
 
 ####################### runnning MMHC #########################
 Y.n.d<- as.data.frame(Y.n)
 colnames(Y.n.d)<- paste0("Y_",1:1000)
   
 bn<- bnlearn::mmhc(as.data.frame(Y.n.d), whitelist = NULL, blacklist = NULL, debug = FALSE)
 
 adj.res.c.mmhc[,,rep]<- amat(bn)
 
####################### runnning PC-alg #########################
 
 pc<- pcAlgo(dm = as.matrix(Y.n), alpha=0.005, corMethod = "standard", verbose=FALSE, directed=FALSE, G=NULL, datatype = "continuous",
             NAdelete=TRUE, m.max=Inf, u2pd = "rand", psepset=FALSE)
 adj.res.c.pcalg[,,rep]<- showAmat(pc@graph)
 
 dagbagm<- c(dagbagm, result_skeleton(adj.res.c[,,rep],network.skel)[1])
 hc<- c(hc, result_skeleton(adj.res.c.noboot[,,rep],network.skel)[1])
 mmhc<- c(mmhc, result_skeleton(adj.res.c.mmhc[,,rep],network.skel)[1])
 pcalg<- c(pcalg, result_skeleton(adj.res.c.pcalg[,,rep],network.skel)[1])
}
mean(dagbagm); mean(hc); mean(mmhc); mean(pcalg)
sd(dagbagm); sd(hc); sd(mmhc); sd(pcalg)

############ Simulation II ###################
n<-250 ## also change it to 100
k=500
p<- 504
load("true_adj_p_504.RData")

gene.pancr.label<-function(N, pancr.adj,beta.simu, SN=runif(nrow(pancr.adj), min=0.5, max=1.5))
{
  ### pancr.adj[i,j] indicates the existence of edge i--->j
  ### para: N -- sample size, pancr.adj: adjacency matrix, bmin -- lower bound for the coefficients, bmax -- upper bound for the coefficients 
  ### SN: a p-vector , signal to noise ratio, used to determine the variance of the residuals 
  ### return: N by p data matrix 
  
  size=nrow(pancr.adj)
  
  result=matrix(0, N, size)
  
  
  tt=which(colSums(pancr.adj)==0)  ### find out nodes without parents 
  for(j in tt){
    result[,j]=rnorm(N, mean=0, sd=1)  ### generate them from N(0,1)
  }
  
  esd.v=rep(0, size)
  esd.v[tt]=1
  
  tt2=(1:size)[-tt]   ##nodes with at least one parents 
  current=tt2
  
  while(length(current)>0){
    
    i=current[1]                      ##look at node i  
    temp=beta.simu[,i]
    temp2=which(temp!=0)     ##parents for node i
    check=numeric(length(temp2))
    
    for(m in 1:length(temp2)){
      check[m]=all(result[,temp2[m]]==0)   ##check whether node i's parents have been generated: check =0 yes, check =1, not yet 
    }
    
    if(any(check!=0)){  ##if some parents of node i have not been generated 
      current=current[-1]
      current=c(current,i)
    }
    
    if(all(check==0)){  ## if all parenets of node i have been generated 
      
      yhat=as.matrix(result[,temp2],nrow=N)%*%as.matrix(temp[temp2])
      ysd=sd(as.vector(yhat))
      esd.v[i]=ysd/SN[i]
      #esd.v[i]<-sd.vec[i]
      error=rnorm(N, mean=0, sd=esd.v[i])
      result[,i]=yhat+error
      current=current[-1]
    } 
    
  }##end while 
  
  return(list(data=result, beta=beta.simu, esd=esd.v))
}

beta.sim<-function(pancr.adj,bmin=0.3, bmax=0.5,flip=FALSE)
{
  size=nrow(pancr.adj)
  
  beta.simu=matrix(0, nrow=size, ncol=size)
  beta.simu[pancr.adj==1]=runif(sum(pancr.adj), min=bmin, max=bmax)
  if(flip==TRUE){
    sign.temp=numeric(sum(pancr.adj==1))
    sign.temp=sample(c(-1,1),length(sign.temp),replace=TRUE)
    beta.simu[pancr.adj==1]=beta.simu[pancr.adj==1]*sign.temp
  }
  return(beta.simu)
}

beta.simu<-beta.sim(true.adj)

data_generation<-function(n,p,pancr.adj,beta.simu)
{
  Y.n=NULL
  panc.simu.Data=gene.pancr.label(n, pancr.adj,beta.simu,SN=runif(nrow(pancr.adj), min=0.5, max=1.5))  ##SNR , 0.5, 1.5 
  
  Y=panc.simu.Data$data
  Y.gm=apply(Y, 2, mean)
  Y.gsd=apply(Y, 2, sd)
  Y.n=(Y-matrix(Y.gm, n, p, byrow=T))/matrix(Y.gsd, n, p, byrow=T)
  
  return(Y.n=Y.n)
}

nrep=100
for(j in 1:nrep){
  print(j)
  Y.n=data_generation(n=n,p=p,pancr.adj=adjn,beta.simu)
  file.name=paste0("data_n_250_p_",504,"_rep",j,".RData")
  save(Y.n,file=file.name)
}

adj.res.c<- array(NA, c(p,p,nrep))
adj.res.c.noboot<- array(NA, c(p,p,nrep))
adj.res.c.mmhc<- array(NA, c(p,p,nrep))
adj.res.c.pcalg<- array(NA, c(p,p,nrep))

dagbagm.power<- NULL
hc.power<- NULL
mmhc.power<- NULL
pcalg.power<- NULL

dagbagm.fdr<- NULL
hc.fdr<- NULL
mmhc.fdr<- NULL
pcalg.fdr<- NULL
for(rep in 1:100)
{
  load(paste0("data_n_250_","p_504_","rep",rep,".RData"))
  
  print(paste0("running ",rep))
  
  ########### running dagbagM with bootstrap ################  
  res.hc.c<- dagbagM::hc_boot_parallel(Y=Y.n, n.boot=100, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)
  
  net.hc.c=dagbag::score_shd(res.hc.c,threshold=0)$adj.matrix
  
  adj.res.c[,,rep]=net.hc.c
  
  ############## running dagbagM without bootstrap ##############
  
 net.hc.c.noboot <- dagbagM::hc(Y=as.matrix(Y.n), nodeType=rep("c",p), blackList=NULL, whiteList=NULL, standardize=TRUE, restart=10,maxStep=1000, verbose=FALSE, tol=1e-06)
  
  adj.res.c.noboot[,,rep]<- net.hc.c.noboot$adjacency
    
  ####################### runnning MMHC #########################
  bn<- bnlearn::mmhc(as.data.frame(Y.n), whitelist = NULL, blacklist = NULL, debug = FALSE)
  
  adj.res.c.mmhc[,,rep]<- amat(bn)
  
  ####################### runnning PC-alg #########################
  
  pc<- pcAlgo(dm = as.matrix(Y.n), alpha=0.005, corMethod = "standard", verbose=FALSE, directed=FALSE, G=NULL, datatype = "continuous",
              NAdelete=TRUE, m.max=Inf, u2pd = "rand", psepset=FALSE)
  adj.res.c.pcalg[,,rep]<- showAmat(pc@graph)
  
  dagbagm.power<- c(dagbagm.power, (result_skeleton(adj.res.c[,,rep],network.skel)[2])/sum(true.adj))
  hc.power<- c(hc.power, (result_skeleton(adj.res.c.noboot[,,rep],network.skel)[2])/sum(true.adj))
  mmhc.power<- c(mmhc.power, (result_skeleton(adj.res.c.mmhc[,,rep],network.skel)[2])/sum(true.adj))
  pcalg.power<- c(pcalg.power, (result_skeleton(adj.res.c.pcalg[,,rep],network.skel)[2])/sum(true.adj))
  
  dagbagm.fdr<- c(dagbagm.fdr, (result_skeleton(adj.res.c[,,rep],network.skel)[1] - result_skeleton(adj.res.c[,,rep],network.skel)[2])/result_skeleton(adj.res.c[,,rep],network.skel)[1])
  
  hc.fdr<- c(hc.fdr, (result_skeleton(adj.res.c.noboot[,,rep],network.skel)[1] - result_skeleton(adj.res.c.noboot[,,rep],network.skel)[2])/result_skeleton(adj.res.c.noboot[,,rep],network.skel)[1])
  
  mmhc.fdr<- c(mmhc.fdr, (result_skeleton(adj.res.c.mmhc[,,rep],network.skel)[1] - result_skeleton(adj.res.c.mmhc[,,rep],network.skel)[2])/result_skeleton(adj.res.c.mmhc[,,rep],network.skel)[1])
  
  pcalg.fdr<- c(pcalg.fdr, (result_skeleton(adj.res.c.pcalg[,,rep],network.skel)[1] - result_skeleton(adj.res.c.pcalg[,,rep],network.skel)[2])/result_skeleton(adj.res.c.pcalg[,,rep],network.skel)[1])
}

mean(dagbagm.power); mean(hc.power); mean(mmhc.power); mean(pcalg.power)
mean(dagbagm.fdr); mean(hc.fdr); mean(mmhc.fdr); mean(pcalg.fdr)

