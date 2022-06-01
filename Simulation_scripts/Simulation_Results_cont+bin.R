###################### summary #################################

########### True adjacency matrix for DAG with p=11 ###############
p<- 11
true.adj<- matrix(0,p,p)
true.adj[1,4]<- true.adj[2,3]<- true.adj[p,5]<- true.adj[1,7]<- true.adj[2,8]<- true.adj[7,p]<- true.adj[8,p]<- true.adj[5,6]<- 1

########### True adjacency matrix for DAG with p=21 ###############
p<- 21
true.adj<- matrix(0,p,p)

true.adj[1,4]<- true.adj[2,3]<- true.adj[p,5]<- true.adj[1,7]<- true.adj[2,8]<- true.adj[7,p]<- true.adj[8,p]<- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

########### True adjacency matrix for DAG with p=61, 121, 201, 501 ###############

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

true.adj<- matrix(0,20,20)
true.adj1<- matrix(0,20,20)

true.adj[1,4]<- true.adj[2,3]<- true.adj[1,7]<- true.adj[2,8] <- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

true.adj1[1,4]<- true.adj1[2,3]<- true.adj1[1,7]<- true.adj1[2,8] <- true.adj1[7,5]<-true.adj1[8,5]<- true.adj1[5,6]<- true.adj1[11,14]<- true.adj1[12,13]<- true.adj1[11,17]<- true.adj1[12,18]<- true.adj1[17,15]<-true.adj1[18,15]<- true.adj1[15,16]<-true.adj1[4,14]<-true.adj1[3,13]<- 1

A<- rbind(cbind(true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), true.adj1,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj1, matrix(0,20,20),matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), true.adj1,matrix(0,20,20)), cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),true.adj1))

A<- rbind(cbind(A, rep(0,120)), rep(0,121))

A[7,p]<- A[8,p]<- A[p,5]<- 1

true.adj<- A

######### p=201  ############
p<- 201

true.adj<- matrix(0,20,20)
true.adj1<- matrix(0,20,20)

true.adj[1,4]<- true.adj[2,3]<- true.adj[1,7]<- true.adj[2,8] <- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

true.adj1[1,4]<- true.adj1[2,3]<- true.adj1[1,7]<- true.adj1[2,8] <- true.adj1[7,5]<-true.adj1[8,5]<- true.adj1[5,6]<- true.adj1[11,14]<- true.adj1[12,13]<- true.adj1[11,17]<- true.adj1[12,18]<- true.adj1[17,15]<-true.adj1[18,15]<- true.adj1[15,16]<-true.adj1[4,14]<-true.adj1[3,13]<- 1

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
######### p=501 ################
p<- 501

true.adj<- matrix(0,20,20)
true.adj1<- matrix(0,20,20)

true.adj[1,4]<- true.adj[2,3]<- true.adj[1,7]<- true.adj[2,8] <- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

true.adj1[1,4]<- true.adj1[2,3]<- true.adj1[1,7]<- true.adj1[2,8] <- true.adj1[7,5]<-true.adj1[8,5]<- true.adj1[5,6]<- true.adj1[11,14]<- true.adj1[12,13]<- true.adj1[11,17]<- true.adj1[12,18]<- true.adj1[17,15]<-true.adj1[18,15]<- true.adj1[15,16]<-true.adj1[4,14]<-true.adj1[3,13]<- 1

A<- rbind(cbind(true.adj,matrix(0,20,480)),cbind(matrix(0,20,20),true.adj1,matrix(0,20,460)),cbind(matrix(0,20,40),true.adj1,matrix(0,20,440)),cbind(matrix(0,20,60),true.adj1,matrix(0,20,420)),cbind(matrix(0,20,80),true.adj1,matrix(0,20,400)),cbind(matrix(0,20,100),true.adj1,matrix(0,20,380)),cbind(matrix(0,20,120),true.adj1,matrix(0,20,360)),cbind(matrix(0,20,140),true.adj1,matrix(0,20,340)),cbind(matrix(0,20,160),true.adj1,matrix(0,20,320)),cbind(matrix(0,20,180),true.adj1,matrix(0,20,300)),cbind(matrix(0,20,200),true.adj1,matrix(0,20,280)),cbind(matrix(0,20,220),true.adj1,matrix(0,20,260)),cbind(matrix(0,20,240),true.adj1,matrix(0,20,240)),cbind(matrix(0,20,260),true.adj1,matrix(0,20,220)),cbind(matrix(0,20,280),true.adj1,matrix(0,20,200)),cbind(matrix(0,20,300),true.adj1,matrix(0,20,180)),cbind(matrix(0,20,320),true.adj1,matrix(0,20,160)),cbind(matrix(0,20,340),true.adj1,matrix(0,20,140)),cbind(matrix(0,20,360),true.adj1,matrix(0,20,120)),cbind(matrix(0,20,380),true.adj1,matrix(0,20,100)),cbind(matrix(0,20,400),true.adj1,matrix(0,20,80)),cbind(matrix(0,20,420),true.adj1,matrix(0,20,60)),cbind(matrix(0,20,440),true.adj1,matrix(0,20,40)),cbind(matrix(0,20,460),true.adj1,matrix(0,20,20)),cbind(matrix(0,20,480),true.adj1))

Nk<- 25
p<- 20*Nk + 1

A<- rbind(cbind(A, rep(0,500)), rep(0,501))

A[7,p]<- A[8,p]<- A[p,5]<- 1

true.adj<- A

######### p=210 ###########
p<- 210

true.adj<- matrix(0,20,20)
true.adj1<- matrix(0,20,20)

true.adj[1,4]<- true.adj[2,3]<- true.adj[1,7]<- true.adj[2,8] <- true.adj[5,6]<- true.adj[11,14]<- true.adj[12,13]<- true.adj[11,17]<- true.adj[12,18]<- true.adj[17,15]<-true.adj[18,15]<- true.adj[15,16]<-true.adj[4,14]<-true.adj[3,13]<- 1

true.adj1[1,4]<- true.adj1[2,3]<- true.adj1[1,7]<- true.adj1[2,8] <- true.adj1[7,5]<-true.adj1[8,5]<- true.adj1[5,6]<- true.adj1[11,14]<- true.adj1[12,13]<- true.adj1[11,17]<- true.adj1[12,18]<- true.adj1[17,15]<-true.adj1[18,15]<- true.adj1[15,16]<-true.adj1[4,14]<-true.adj1[3,13]<- 1

A<- rbind(cbind(true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20), true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj, matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)), 
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20), true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),true.adj,matrix(0,20,20),matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),true.adj,matrix(0,20,20),matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),true.adj,matrix(0,20,20)),
          cbind(matrix(0,20,20), matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20), matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),matrix(0,20,20),true.adj))

Nk<- 10
p<- 20*Nk + Nk
Ay<- matrix(0,10,10)

A<- rbind(cbind(A, matrix(0,200,10)), cbind(matrix(0,10,200),Ay))

for(k in 1:10){
  
  foo<- 20*Nk+k
  temp1<- 20*(k-1)+5
  temp2<- 20*(k-1)+7
  temp3<- 20*(k-1)+8
  
  A[temp2,foo]<- A[temp3,foo]<- A[foo,temp1]<- 1
}

true.adj<- A

############ Calculate power and FDR for skeleton edges #################################

true.ske=(true.adj==1)|(t(true.adj)==1)
sum(true.adj)
sum(true.ske)/2

nrep=dim(adj.res.m)[3]
power.ske=numeric(nrep)
fdr.ske=numeric(nrep)

for (rep in 1:nrep){
  adj.cur=adj.res.m[,,rep]
  ske.cur=(adj.cur==T)|(t(adj.cur)==T)
  power.ske[rep]=sum(ske.cur==T&true.ske==T)/sum(true.ske==T) 
  fdr.ske[rep]=sum(ske.cur==T&true.ske==F)/sum(ske.cur==T) 
}

c(mean(power.ske), mean(fdr.ske)) ### power and FDR

############ Calculate power and FDR for directed edges #################################

true.ske=true.adj
sum(true.ske)

nrep=dim(adj.res.m)[3]
power.dir=numeric(nrep)
fdr.dir=numeric(nrep)

for (rep in 1:nrep){
  adj.cur=adj.res.m[,,rep]
  diag(adj.cur[1:p,1:p])<-0
  ske.cur=adj.cur
  power.dir[rep]=sum(ske.cur==T&true.ske==T)/sum(true.ske==T) 
  fdr.dir[rep]=sum(ske.cur==T&true.ske==F)/sum(ske.cur==T) #FDR of 
}

c(mean(power.dir, na.rm = T), mean(fdr.dir, na.rm = T)) 

########## power and fdr for directed edges between continuous and binary nodes for p=11,21,61,121,201,501 ###########

st.m<-NULL
sf.m<- NULL
for(k in 1:100)
{
  st.m<-c(st.m,((adj.res.m[7,p,k]+adj.res.m[8,p,k]+adj.res.m[p,5,k])/3))
  sf.m<- c(sf.m, (sum(adj.res.m[-c(7,8),p,k])+sum(adj.res.m[p,-c(5),k]))/(sum(adj.res.m[,p,k])+sum(adj.res.m[p,,k])))
}

##### change st.m and sf.m to st.c and sf.c (for DAGBagC), st.d and sf.d (for bnlearnd) and st.md and sf.md (for mDAG) and change adj.res.m accordingly.

m.comb.tp<- mean(st.m)
m.comb.fp<- mean(sf.m)

c.comb.tp<- mean(st.c)
c.comb.fp<- mean(sf.c)

md.comb.tp<- mean(st.md)
md.comb.fp<- mean(sf.md)

d.comb.tp<- mean(st.d)
d.comb.fp<- mean(sf.d)

df.comb<- cbind.data.frame("TP"=c(m.comb.tp,c.comb.tp,md.comb.tp,d.comb.tp),"FP"=c(m.comb.fp,c.comb.fp,md.comb.fp,d.comb.fp),"method"=c("DAGBagM","DAGBagC","mDAG","bnlearnD")) 
######## power and FDR for directed edges between continuous and binary nodes for p=210 ############
Nk<- 10
st.m<-NULL
sf.m<- NULL
for(k in 1:100)
{
  s1<-0
  s3<-0
  s3.f<-0
  for(i in 1:10){
    foo<- 20*Nk+i
    temp1<-20*(i-1)+5
    temp2<-20*(i-1)+7
    temp3<-20*(i-1)+8
    
    s1<- s1+adj.res.m[temp2,foo,k]+adj.res.m[temp3,foo,k]+adj.res.m[foo,temp1,k]
    s3<- s3+sum(adj.res.m[-c(temp2,temp3),foo,k])+sum(adj.res.m[foo,-c(temp1),k])
    
    s3.f<- s3.f + sum(adj.res.m[,foo,k])+sum(adj.res.m[foo,,k])
  }
  
  st.m<- c(st.m,s1/30)
  sf.m<- c(sf.m,s3/s3.f)
}

##### change st.m and sf.m to st.c and sf.c (for DAGBagC), st.d and sf.d (for bnlearnd) and st.md and sf.md (for mDAG) and change adj.res.m accordingly.

m.comb.tp<- mean(st.m)
m.comb.fp<- mean(sf.m)

c.comb.tp<- mean(st.c)
c.comb.fp<- mean(sf.c)

md.comb.tp<- mean(st.md)
md.comb.fp<- mean(sf.md)

d.comb.tp<- mean(st.d)
d.comb.fp<- mean(sf.d)

df.comb<- cbind.data.frame("TP"=c(m.comb.tp,c.comb.tp,md.comb.tp,d.comb.tp),"FP"=c(m.comb.fp,c.comb.fp,md.comb.fp,d.comb.fp),"method"=c("DAGBagM","DAGBagC","mDAG","bnlearnD")) 