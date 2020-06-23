source("network_gene_function_10172012.R")
load("Y_parent_child.RData") #### p=12

load("Y_parent_child_22.RData") #### p=22

# adjn<-cbind(true.adj[,1:20],c(rep(0,10),1,rep(0,10)),true.adj[,21])
# adjn1<-rbind(adjn[1:20,],rep(0,22),adjn[21,])
# adjn1[22,21]<-1
# colnames(adjn1)[21:22]<-rownames(adjn1)[21:22]<-c("X21","Y")
#   
# dim(adjn1)
# sum(adjn1)
# obj<- graph_from_adjacency_matrix(adjn1, mode = "directed", weighted = TRUE)
# is.dag(obj)
# plot(obj)

#true.adj<-adjn1
#save(true.adj,file="Y_parent_child_22.RData")

n<- 102

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

true.dir<-true.adj

beta.simu<-beta.sim(true.dir)

beta0.y<-runif(1,1,3)
beta2.y<-runif(1,1,3)
beta3.y<-runif(1,1,3)

beta.y<-c(beta0.y,beta2.y,beta3.y)

data_generation<-function(n,p,pancr.adj,beta.simu,beta.y)
{
  Y.n=NULL
  panc.simu.Data=gene.pancr.label(n, pancr.adj,beta.simu,SN=runif(nrow(pancr.adj), min=0.5, max=1.5))  ##SNR , 0.5, 1.5 
  
  Y=panc.simu.Data$data[,-22]
  Y.gm=apply(Y, 2, mean)
  Y.gsd=apply(Y, 2, sd)
  Y.n=(Y-matrix(Y.gm, n, p, byrow=T))/matrix(Y.gsd, n, p, byrow=T)

  ##### generate discrete node #####
  p.y<-exp(beta0.y + beta2.y*Y.n[,1] + beta3.y*Y.n[,2])/(1+exp(beta.y[1] + beta.y[2]*Y.n[,1] + beta.y[3]*Y.n[,2]))
  Y.disc<-vapply(p.y, rbinom, FUN.VALUE = 0,n=1, size=1)
  
  
  Y.n<-cbind(Y.n,Y.disc)
  return(list(Y.n=Y.n,SN=panc.simu.Data$SN))
}

SN<-NULL
nrep=100
set.seed(2000)
seeds=sample(2:10000,nrep)
for(rep in 1:nrep){
  print(rep)
  set.seed(seeds[rep])
  Y.n=data_generation(n=102,p=21,pancr.adj=true.dir,beta.simu,beta.y)
  SN<-c(SN,Y.n$SN[21])
  file.name=paste("bin_par_child_22/data_n102_","p_22_","rep",rep,".RData",sep='')
  t=Y.n$Y.n
  save(t,file=file.name)
}

##################### dagbag/bnlearnC ####################
iter<-100
n=102

for(j in 1:iter)
{
  load(paste("bin_par_child_22/data_n102_","p_22_","rep",j,".RData",sep=''))
  
  print(paste("iter",j,"is running"),sep="")
  x1<-rnorm(n,0,0.5)
  x2<-rnorm(n,0,0.5)
  
  x1.s<-(x1-mean(x1))/sd(x1) ### standardize x1
  x2.s<-(x2-mean(x2))/sd(x2) ### standardize x2
  
  beta0.y<-runif(1,1,3)
  beta1.y<-runif(1,1,3)
  beta2.y<-runif(1,1,3)
  
  p.y<-exp(beta0.y + beta1.y*x1.s + beta2.y*x2.s)/(1+exp(beta0.y + beta1.y*x1.s + beta2.y*x2.s))
  Y.disc<-vapply(p.y, rbinom, FUN.VALUE = 0,n=1, size=1)
  sum(Y.disc)
  
  Y.n<-t
  
  yhat=as.matrix(Y.n[,c(11,22)],nrow=102)%*%as.matrix(beta.simu[c(11,22),21])
  ysd=sd(as.vector(yhat))
  esd.v=ysd/SN[j]
  error=rnorm(102, mean=0, sd=esd.v)
  x21=yhat+error
  x21.s<-(x21-mean(x21))/sd(x21)
  
  Y<-cbind("x1"=x1.s,"x2"=x2.s,Y.n[,3:20],"x21"=x21.s,"y"=Y.disc)
  
  file.name=paste("bin_par_child_22/data_n102_","p_22_","rep",j,".RData",sep='')
  save(Y,file=file.name)
}


################## for bnlearnD #######################
for(k in 1:100)
{
  load(paste("bin_par_child_22/data_n102_","p_22_","rep",k,".RData",sep=''))
  Y.disc<-matrix(0,102,21)
  for(i in 1:21)
  {
    temp1<-paste("a_",i,sep="")
    temp2<-paste("b_",i,sep="")
    
    Y.disc[,i]<-ifelse((Y[,i]>median(Y[,i])),temp1,temp2)
  }
  
  Y.disc<-cbind.data.frame(Y.disc,Y[,22])
  Y.disc[,22]<-as.factor(Y.disc[,22])
  
  file.name=paste("bin_par_child_22_disc/datad_n102_","p_22_","rep",k,".Rdata",sep='')
  
  save(Y.disc,file=file.name)
}

################### test run ######################

source("C:\\Users\\chowds14\\Desktop\\MSSM\\MSSM\\PTRC\\dag\\dagbag_functions_scoreBIC.R")
library(dagbag)

seed.u<-100
seed.e<-0
blacklist<-NULL
n.B<-100 ## no. of bootstrap samples
iter<-1 ### number of simulations
n<-102
p<-12

dagbag.list<-list()
hc.list<-list()

for(j in 1:iter)
{
  print(paste("iter",j,"is running"),sep="")
  
  require(bnlearn)
  
  n=nrow(Y)
  p=ncol(Y)
  
  hc.adj=array(0,c(p,p,n.B))
  dagbag_bin.adj=array(0,c(p,p,n.B))
  
  Y.B=array(0,c(n,p,n.B))
  
  i=1
  while(i < (n.B+1)){
    
    set.seed(i*seed.u+seed.e)
    
    s.pick=sample(1:n, n, replace=TRUE)
    node.rand<-sample(1:p,replace=F) ### randomize node order
    
    Y.B[,,i]=Y[s.pick,node.rand]
    
    x11<-Y.B[,which(node.rand==1),i]
    x22<-Y.B[,which(node.rand==2),i]
    y<-Y.B[,which(node.rand==p),i]
    
    yy<-y
    for(k in 1:n)
    {
      if(x11[k]< -0.5)
      {
        yy[k]<-0
      }else if(x22[k]> 0.5)
      {
        yy[k]<-1
      }
    }
    
    Y.B[,which(node.rand==p),i]<- yy
    
    temp<-c(rep("c",(p-1)),"b") ## change "b" to "c" to treat binary as continuous in dagbag
    
    node.type<-temp[node.rand]

    result<-try(score_BIC_opt(Y=Y.B[,,i], node.type = node.type, blacklist=blacklist, threshold = 0, nstart=10, seed=Sys.time(), print=FALSE),silent=T)
    if(class(result)=="try-error"){
      next;
    }else{
      print(i)
      s<-sort(node.rand,decreasing=F,ind=T)
      
      dagbag_bin.adj[,,i]<-result$adj.matrix[s$ix,s$ix]
      
      hc.B=hc(as.data.frame(Y.B[,,i]),restart=10)
      hc.adj[,,i]=amat(hc.B)[s$ix,s$ix]
      
      i=i+1
    }
  }
  
  hc.list[[j]]<-score_shd(hc.adj,threshold=0)$adj.matrix
  remove(hc.adj)
  
  dagbag.list[[j]]<-score_shd(dagbag_bin.adj,threshold=0)$adj.matrix
  remove(dagbag_bin.adj)
}

