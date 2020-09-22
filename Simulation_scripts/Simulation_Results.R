###################### summary #################################

########### True adjacency matrix for DAG with p=11 ###############
p<- 11
true.adj<- matrix(0,p,p)
true.adj[1,4]<- true.adj[2,3]<- true.adj[p,5]<- true.adj[1,7]<- true.adj[2,8]<- true.adj[7,p]<- true.adj[8,p]<- true.adj[5,6]<- 1

########### True adjacency matrix for DAG with p=21 ###############
p<- 21
true.adj<- matrix(0,p,p)

true.adj[1,4]<- true.adj[2,3]<- true.adj[p,5]<- true.adj[1,7]<- true.adj[2,8]<- true.adj[7,p]<- true.adj[8,p]<- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

########### True adjacency matrix for DAG with p=61, 121, 201 ###############

########## p=61 ##############
p<- 61
true.adj<- matrix(0,20,20)
true.adj1<- matrix(0,20,20)

true.adj[1,4]<- true.adj[2,3]<- true.adj[1,7]<- true.adj[2,8] <- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

true.adj1[1,4]<- true.adj1[2,3]<- true.adj1[1,7]<- true.adj1[2,8] <- true.adj1[7,5]<-true.adj1[8,5]<- true.adj1[5,6]<- true.adj1[11,14]<- true.adj1[12,13]<- true.adj1[11,17]<- true.adj1[12,18]<- true.adj1[17,15]<-true.adj1[18,15]<- true.adj1[15,16]<-true.adj1[4,14]<-true.adj1[3,13]<- 1

A<- rbind(cbind(true.adj,matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), true.adj1,matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), true.adj1))

A<- rbind(cbind(A, rep(0,60)), rep(0,61))

A[7,p]<- A[8,p]<- A[p,5]<- 1

true.adj<- A

############# p=121 #############
p<- 121

A<- rbind(cbind(true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj1, matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), true.adj1,matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),true.adj1))

A<- rbind(cbind(A, rep(0,120)), rep(0,121))

A[7,p]<- A[8,p]<- A[p,5]<- 1

true.adj<- A

######### p=201  ############
p<- 201
A<- rbind(cbind(true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj1, matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj1,matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),true.adj1,matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),true.adj1))

A<- rbind(cbind(A, rep(0,200)), rep(0,201))

A[7,p]<- A[8,p]<- A[p,5]<- 1

true.adj<- A
############ Calculate power and FDR #################################

true.ske=(true.adj==1)|(t(true.adj)==1)
sum(true.adj)
sum(true.ske)/2

nrep=dim(adj.res.m)[3]
power.ske=numeric(nrep)
fdr.ske=numeric(nrep)

for (rep in 1:nrep){
  adj.cur=adj.res.c[,,rep]
  ske.cur=(adj.cur==T)|(t(adj.cur)==T)
  power.ske[rep]=sum(ske.cur==T&true.ske==T)/sum(true.ske==T) 
  fdr.ske[rep]=sum(ske.cur==T&true.ske==F)/sum(ske.cur==T) 
}

c(mean(power.ske), mean(fdr.ske)) ### power and FDR

########### number of times true edge between binary and continuous nodes are detected with correct direction #################

s<-NULL
for(k in 1:100)
{
  s<-c(s,(adj.res.c[7,p,k]+adj.res.c[8,p,k]))  
}
sum(s==2) 
sum(s>0) 