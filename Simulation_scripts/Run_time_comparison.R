library(igraph)
library(dagbag)
library(dagbagM)
library(bnlearn)

################################## data generation #########################

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

for(i in p)
{
  load(paste0("adjn",i,".RData"))
  print(paste0("running ",i))
  n<-500
  p1<-nrow(adjn)
  
  beta.simu<-beta.sim(adjn)
  
  nrep=100
  set.seed(2000)
  seeds=sample(2:10000,nrep)
  for(j in 1:nrep){
    print(j)
    set.seed(seeds[j])
    Y.n=data_generation(n=n,p=p1,pancr.adj=adjn,beta.simu)
    file.name=paste0("data_n_500_p_",i,"/data_rep",j,".RData")
    save(Y.n,file=file.name)
  }
}

######################## different p ########################

p<-c(504,
     1008,
     1512,
     2016,
     2520,
     3024,
     3528,
     4032,
     4536,
     5040
)

########################## running dagbagM ####################  
time.dagbag.final<-NULL
for(i in p)
{
  print(paste0("running ",i))
  
  time.dag.vec<-NULL
 for(j in 1:100)
 {
   load(paste0("data_n_500_p_",i,"/data_rep",j,".RData"))
   print(paste0(c(j,dim(Y.n))))
   
   init<-proc.time()
   
   dag<-dagbagM::hc(Y=as.matrix(Y.n), nodeType=rep("c",p), blackList=NULL, whiteList=NULL, standardize=TRUE, restart=0,maxStep=1000, verbose=FALSE, tol=1e-06)

   time.dag<-(proc.time() - init)
   time.dag.vec<-c(time.dag.vec,time.dag[3])
 }

  time.dagbag.final<-c(time.dagbag.final,mean(time.dag.vec))
} 


####################### running bnlearn ####################  

time.bn.final<-NULL

for(i in p)
{
  print(paste0("running ",i))
  
  time.bn.vec<-NULL
  
  for(j in 1:100)
  {
  
    load(paste0("data_n_500_p_",i,"/data_rep",j,".RData"))
    
  data<-as.data.frame(Y.n)
  
   init<-proc.time()
  
  bn<-hc(data,restart = 0, perturb = 0,max.iter = 1000)
  adj<-amat(bn)

  time.bn<-(proc.time() - init)[3] 
  
  time.bn.vec<-c(time.bn.vec,time.bn)
  }
  time.bn.final<-c(time.bn.final,mean(time.bn.vec))
}

