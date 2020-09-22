data.p11<- function(n){
  for(rep in 1:nrep){
    print(paste("replicate ",  rep))
    
    k=500
    set.seed(k*1001+rep)
    x1<-rnorm(n,0,0.5)
    x2<-rnorm(n,0,0.5)
    x9<-rnorm(n,0,0.5)
    x10<-rnorm(n,0,0.5)
    
    x1.s<-(x1-mean(x1))/sd(x1) ### standardize x1
    x2.s<-(x2-mean(x2))/sd(x2) ### standardize x2
    x9.s<-(x9-mean(x9))/sd(x9) ### standardize x9
    x10.s<-(x10-mean(x10))/sd(x10) ### standardize x10
    
    yhat<- 0 + x1*0.5
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x7<- yhat+error
    x7.s<-(x7-mean(x7))/sd(x7)
    
    yhat<- 0 + x2*0.45
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x8<- yhat+error
    x8.s<-(x8-mean(x8))/sd(x8)
    
    yhat<- 0 + x2*0.4
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x3<- yhat+error
    x3.s<-(x3-mean(x3))/sd(x3)
    
    yhat<- 0 + x1*0.45
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x4<- yhat+error
    x4.s<-(x4-mean(x4))/sd(x4)
    
    beta0.y<-2
    beta1.y<-8
    beta2.y<-2
    
    p.y<-exp(beta0.y + beta1.y*x7.s + beta2.y*x8.s)/(1+exp(beta0.y + beta1.y*x7.s + beta2.y*x8.s))
    Y.disc<-vapply(p.y, rbinom, FUN.VALUE = 0,n=1, size=1)
    
    yhat<- 0 + Y.disc*0.55 
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x5<- yhat+error
    x5.s<-(x5-mean(x5))/sd(x5)
    
    yhat<- 0 + x5*1.2 
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x6<- yhat+error
    x6.s<-(x6-mean(x6))/sd(x6)
    
    Y.n<-cbind("x1"=x1.s,"x2"=x2.s, "x3"=x3.s, "x4"=x4.s, "x5"=x5.s, "x6"=x6.s, "x7"=x7.s, "x8"=x8.s, "x9"=x9.s, "x10"=x10.s, "y"=Y.disc)
    
    file.name<- paste0("data_n",n,"_p",p,"_rep",rep,".RData")
    save(Y.n, file=file.name)
}
}

data.p21<- function(n){
  for(rep in 1:nrep){
    print(paste("replicate ",  rep))
    
    k=500
    set.seed(k*1001+rep)
    
    ######################### 11 nodes ############################
    x1<-rnorm(n,0,0.5)
    x2<-rnorm(n,0,0.5)
    x9<-rnorm(n,0,0.5)
    x10<-rnorm(n,0,0.5)
    
    x1.s<-(x1-mean(x1))/sd(x1) ### standardize x1
    x2.s<-(x2-mean(x2))/sd(x2) ### standardize x2
    x9.s<-(x9-mean(x9))/sd(x9) ### standardize x9
    x10.s<-(x10-mean(x10))/sd(x10) ### standardize x10
    
    yhat<- 0 + x1*0.5
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x7<- yhat+error
    x7.s<-(x7-mean(x7))/sd(x7)
    
    yhat<- 0 + x2*0.45
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x8<- yhat+error
    x8.s<-(x8-mean(x8))/sd(x8)
    
    yhat<- 0 + x2*0.4
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x3<- yhat+error
    x3.s<-(x3-mean(x3))/sd(x3)
    
    yhat<- 0 + x1*0.45
    ysd=sd(as.vector(yhat))
    esd.v=ysd/.75
    error=rnorm(n, mean=0, sd=esd.v)
    x4<- yhat+error
    x4.s<-(x4-mean(x4))/sd(x4)
    
    beta0.y<-2
    beta1.y<-8
    beta2.y<-2
    
    p.y<-exp(beta0.y + beta1.y*x7.s + beta2.y*x8.s)/(1+exp(beta0.y + beta1.y*x7.s + beta2.y*x8.s))
    Y.disc<-vapply(p.y, rbinom, FUN.VALUE = 0,n=1, size=1)
    
    yhat<- 0 + Y.disc*0.55 
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.85
    error=rnorm(n, mean=0, sd=esd.v)
    x5<- yhat+error
    x5.s<-(x5-mean(x5))/sd(x5)
    
    yhat<- 0 + x5*1.2 
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x6<- yhat+error
    x6.s<-(x6-mean(x6))/sd(x6)
    
    ############ additional nodes #############
    
    x11<-rnorm(n,0,0.5)
    x12<-rnorm(n,0,0.5)
    x19<-rnorm(n,0,0.5)
    x20<-rnorm(n,0,0.5)
    
    x11.s<-(x11-mean(x11))/sd(x11) ### standardize x1
    x12.s<-(x12-mean(x12))/sd(x12) ### standardize x2
    x19.s<-(x19-mean(x19))/sd(x19) ### standardize x9
    x20.s<-(x20-mean(x20))/sd(x20) ### standardize x10
    
    yhat<- 0 + x11*0.5
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x17<- yhat+error
    x17.s<-(x17-mean(x17))/sd(x17)
    
    yhat<- 0 + x12*0.45
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x18<- yhat+error
    x18.s<-(x18-mean(x18))/sd(x18)
    
    yhat<- 0 + x12*0.4 + x3*.5
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x13<- yhat+error
    x13.s<-(x13-mean(x13))/sd(x13)
    
    yhat<- 0 + x11*0.45 + x4*.55
    ysd=sd(as.vector(yhat))
    esd.v=ysd/.75
    error=rnorm(n, mean=0, sd=esd.v)
    x14<- yhat+error
    x14.s<-(x14-mean(x14))/sd(x14)
    
    yhat<- 0 + x17*0.55 + x18*.5
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.75
    error=rnorm(n, mean=0, sd=esd.v)
    x15<- yhat+error
    x15.s<-(x15-mean(x15))/sd(x15)
    
    yhat<- 0 + x15*1.2 
    ysd=sd(as.vector(yhat))
    esd.v=ysd/0.8
    error=rnorm(n, mean=0, sd=esd.v)
    x16<- yhat+error
    x16.s<-(x16-mean(x16))/sd(x16)
    
    ###########################################  
    Y.n<-cbind("x1"=x1.s,"x2"=x2.s, "x3"=x3.s, "x4"=x4.s, "x5"=x5.s, "x6"=x6.s, "x7"=x7.s, "x8"=x8.s, "x9"=x9.s, "x10"=x10.s,"x11"=x11.s,"x12"=x12.s, "x13"=x13.s, "x14"=x14.s, "x15"=x15.s, "x16"=x16.s, "x17"=x17.s, "x18"=x18.s, "x19"=x19.s, "x20"=x20.s, "y"=Y.disc)
    
    file.name<- paste0("data_n",n,"_p",p,"_rep",rep,".RData")
    save(Y.n, file=file.name)
  }
}

data.p<- function(n,p){
  for(rep in 1:nrep){
    print(paste("replicate ",  rep))
    
    s=500
    set.seed(s*1001+rep)
    
    ######################### 11 nodes ############################
    beta0.y<-2
    beta1.y<-8
    beta2.y<-2
    
    Nk<- (p-1)/20
    
    Y<- list() 
    for(k in 1:Nk){
      
      temp=c(20*(k-1)+1, 20*(k-1)+2,20*(k-1)+9, 20*(k-1)+10,20*(k-1)+11, 20*(k-1)+12,20*(k-1)+19, 20*(k-1)+20)
      
      for(i in temp){
        t<- rnorm(n,0,0.5)
        assign(paste0("x",i),t)
        
        foo<- get(paste0("x",i))
        assign(paste0("x",i,".s"), (foo-mean(foo))/sd(foo))
      }    
      
      ################################ 7 and 4 ###################        
      
      foo1<- get(paste0("x",temp[1]))
      yhat<- 0 + foo1*0.5
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[1]+6), yhat+error)
      foo11<- get(paste0("x",temp[1]+6))
      assign(paste0("x",temp[1]+6,".s"), (foo11-mean(foo11))/sd(foo11))
      
      yhat<- 0 + foo1*0.45
      ysd=sd(as.vector(yhat))
      esd.v=ysd/.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[1]+3), yhat+error)
      foo11<- get(paste0("x",temp[1]+3))
      assign(paste0("x",temp[1]+3,".s"), (foo11-mean(foo11))/sd(foo11))
      
      ################################ 8 and 3 ################### 
      foo2<- get(paste0("x",temp[2]))
      
      yhat<- 0 + foo2*0.45
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.8
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[2]+6), yhat+error)
      foo11<- get(paste0("x",temp[2]+6))
      assign(paste0("x",temp[2]+6,".s"), (foo11-mean(foo11))/sd(foo11))
      
      yhat<- 0 + foo2*0.4
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[2]+1), yhat+error)
      foo11<- get(paste0("x",temp[2]+1))
      assign(paste0("x",temp[2]+1,".s"), (foo11-mean(foo11))/sd(foo11))
      
      if(k==1){
        
        foo7<- get(paste0("x",7,".s"))
        foo8<- get(paste0("x",8,".s"))
        ################################# y ############################ 
        
        p.y<-exp(beta0.y + beta1.y*foo7 + beta2.y*foo8)/(1+exp(beta0.y + beta1.y*foo7 + beta2.y*foo8))
        Y.disc<-vapply(p.y, rbinom, FUN.VALUE = 0,n=1, size=1)
        
        ################################# 5 and 6 #######################
        
        yhat<- 0 + Y.disc*0.55 
        ysd=sd(as.vector(yhat))
        esd.v=ysd/0.85
        error=rnorm(n, mean=0, sd=esd.v)
        assign(paste0("x",temp[1]+4), yhat+error)
        foo11<- get(paste0("x",temp[1]+4))
        assign(paste0("x",temp[1]+4,".s"), (foo11-mean(foo11))/sd(foo11))
      }else{
        ################################# 5 and 6 #######################
        foo7<- get(paste0("x",temp[1]+6))
        foo8<- get(paste0("x",temp[2]+6))
        
        yhat<- 0 + foo7*0.55 + foo8*0.5 
        ysd=sd(as.vector(yhat))
        esd.v=ysd/0.85
        error=rnorm(n, mean=0, sd=esd.v)
        assign(paste0("x",temp[1]+4), yhat+error)
        foo11<- get(paste0("x",temp[1]+4))
        assign(paste0("x",temp[1]+4,".s"), (foo11-mean(foo11))/sd(foo11))
      }  
      
      foo5<- get(paste0("x",temp[1]+4))
      yhat<- 0 + x5*1.2 
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.8
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[1]+5), yhat+error)
      foo11<- get(paste0("x",temp[1]+5))
      assign(paste0("x",temp[1]+5,".s"), (foo11-mean(foo11))/sd(foo11))    
      
      ########################### 17 and 18 ############################ 
      
      f11<- get(paste0("x",temp[5]))
      yhat<- 0 + foo11*0.5
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[5]+6), yhat+error)
      foo11<- get(paste0("x",temp[5]+6))
      assign(paste0("x",temp[5]+6,".s"), (foo11-mean(foo11))/sd(foo11))
      
      
      foo12<- get(paste0("x",temp[6]))
      yhat<- 0 + foo12*0.45
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.8
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[6]+6), yhat+error)
      foo11<- get(paste0("x",temp[6]+6))
      assign(paste0("x",temp[6]+6,".s"), (foo11-mean(foo11))/sd(foo11))
      
      ######################## 13 and 14 #########################
      
      foo3<- get(paste0("x",temp[1]+2))
      yhat<- 0 + foo12*0.4 + foo3*.5
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[5]+2), yhat+error)
      foo11<- get(paste0("x",temp[5]+2))
      assign(paste0("x",temp[5]+2,".s"), (foo11-mean(foo11))/sd(foo11))
      
      foo4<- get(paste0("x",temp[2]+2))
      yhat<- 0 + f11*0.45 + foo4*.55
      ysd=sd(as.vector(yhat))
      esd.v=ysd/.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[6]+2), yhat+error)
      foo11<- get(paste0("x",temp[6]+2))
      assign(paste0("x",temp[6]+2,".s"), (foo11-mean(foo11))/sd(foo11))
      
      ############################ 15 and 16 #############################    
      
      foo17<- get(paste0("x",temp[5]+6))
      foo18<- get(paste0("x",temp[6]+6))
      yhat<- 0 + foo17*0.55 + foo18*.5
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.75
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[5]+4), yhat+error)
      foo11<- get(paste0("x",temp[5]+4))
      assign(paste0("x",temp[5]+4,".s"), (foo11-mean(foo11))/sd(foo11))
      
      foo15<- get(paste0("x",temp[5]+4))
      yhat<- 0 + foo15*1.2 
      ysd=sd(as.vector(yhat))
      esd.v=ysd/0.8
      error=rnorm(n, mean=0, sd=esd.v)
      assign(paste0("x",temp[5]+5), yhat+error)
      foo11<- get(paste0("x",temp[5]+5))
      assign(paste0("x",temp[5]+5,".s"), (foo11-mean(foo11))/sd(foo11))
      
      ###########################################  
      mat<- NULL
      
      start<- (k-1)*20 +1
      end<- 20*k
      
      for(i in start:end){
        var.name<- get(paste0("x",i,".s"))
        mat<- cbind(mat,var.name)               
      }
      colnames(mat)<- paste0("x",start:end)
      Y[[k]]<- mat
    }  
    
    Y.n<- NULL  
    for(k in 1:Nk){
      Y.n<- cbind(Y.n, Y[[k]])
    }
    Y.n<- cbind(Y.n, "y"=Y.disc)
    file.name<- paste0("data_n",n,"_p",p,"_rep",rep,".RData")
    save(Y.n, file=file.name)  
  }
}

data.gen<- function(n,p){
  
  if(p==11){
    data.p11(n=n)
  }else if(p==21){
    data.p21(n=n)
  }else if(p>21){
    data.p(n=n,p=p)
  }
}