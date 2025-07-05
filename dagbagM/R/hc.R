hc <- function(Y, nodeType=NULL, whiteList=NULL, blackList=NULL, standardize= TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L, verbose = FALSE) {
##Y: n by p data matrix;  nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##standardize: whether standardize continuous nodes to have mean zero and sd. 1
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for breaking score tie;
##return: HC learned adjacency matrix on Y

	p=ncol(Y)
	if(is.null(whiteList)){
		whiteList=matrix(FALSE, p,p)
	}else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
      	stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

	}

	if(is.null(blackList)){
		blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
	}else{##check black list format 
		if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
      	stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
	}	

	if(is.null(nodeType)){
		nodeType=rep("c",p)
	}else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
       		stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

	}

  if(standardize){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes
      Y[,i]=(Y[,i]-mean(Y[,i]))/sd(Y[,i])
     }
    }
  }

    hc_(Y, nodeType, whiteList, blackList, tol, maxStep, restart, seed, verbose)
}


hc_boot<-function(Y, n.boot=1, nodeType=NULL,  whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L,  nodeShuffle=TRUE, verbose = FALSE){
##Y: n by p data matrix;  n.boot: integer, number of boostrap resamples; nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##standardize:  whether standardize continuous nodes to have mean zero and sd. 1
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for boostrap resampling, nodeShuffling
##nodeShuffle: whether (T) or not (F) to perform node shuffle (to avoid bias due to search order and score ties)  
##return: HC learned adjacency matrices on bootstrap resamples 
   p=ncol(Y)
   n=nrow(Y)

   ##check argument type 
   if(is.null(whiteList)){
		whiteList=matrix(FALSE, p,p)
	}else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
      	stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

	}

	if(is.null(blackList)){
		blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
	}else{##check black list format 
		if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
      	stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
	}	

	if(is.null(nodeType)){
		nodeType=rep("c",p)
	}else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
       		stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

	}

  if(standardize){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes
      Y[,i]=(Y[,i]-mean(Y[,i]))/sd(Y[,i])
     }
    }
  }

   ## hc on bootstrap resamples  
   adj=array(NA,c(p,p,n.boot))

   for(i in 1:n.boot){
   print(paste("fit bootstrap sample ", i, "..."))
   set.seed(i*1001+seed)
   s.pick=sample(1:n, n, replace=TRUE) ##resample data 
  
   if(nodeShuffle){##shuffle node order 
   	node.rand=sample(1:p, p, replace=FALSE)  
   	node.index=sort(node.rand,decreasing=F,ind=T)$ix 
   }else{##not shuffle node order
   	node.rand=1:p
   	node.index=1:p
   }

   Y.B=Y[s.pick, node.rand]
   node.B=nodeType[node.rand]
   whiteList.B=whiteList[node.rand, node.rand]
   blackList.B=blackList[node.rand, node.rand]

   curRes=hc_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep, restart, i*11+seed, verbose)
   adjRes=curRes$adjacency
   adj[,,i]=adjRes[node.index, node.index]
   }

   return(adj)
}

############################
##parallel version of hc_boot: need "doParallel" and "foreach" packages
hc_boot_parallel<-function(Y, n.boot=1, nodeType=NULL,  whiteList=NULL, blackList=NULL, standardize=TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L,  nodeShuffle=TRUE, numThread=2, verbose = FALSE){
##Y: n by p data matrix;  n.boot: integer, number of boostrap resamples; nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##standardize:  whether standardize continuous nodes to have mean zero and sd. 1
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for boostrap resampling, nodeShuffling
##nodeShuffle: whether (T) or not (F) to perform node shuffle (to avoid bias due to search order and score ties)  
##numThread: number of threads to use in parallel computing
##return: HC learned adjacency matrices on bootstrap resamples 
  #
  library(parallel)
  library(foreach)
  library(future)
  library(doFuture)
  
  #
   p=ncol(Y)
   n=nrow(Y)

   ##check argument type 
   if(is.null(whiteList)){
    whiteList=matrix(FALSE, p,p)
    }else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
        stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

  }

  if(is.null(blackList)){
    blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
    }else{##check black list format 
      if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
        stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
  } 

  if(is.null(nodeType)){
    nodeType=rep("c",p)
    }else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
          stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

  }

  if(standardize){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes
      Y[,i]=(Y[,i]-mean(Y[,i]))/sd(Y[,i])
     }
    }
  }

   ## hc on bootstrap resamples  
   plan(multisession, workers = min(numThread, detectCores()-1)) ## use plan(multisession): create independent R sessions, compatible with Rcpp calls
   
   
   result<-foreach(i = 1:n.boot,  .errorhandling = "stop",
                   .options.future = list(seed = TRUE, packages=c("dagbagM"))) %dofuture% {
   
       set.seed(i*1001+seed)
       s.pick=sample(1:n, n, replace=TRUE) ##resample data 
      
       if(nodeShuffle){##shuffle node order 
        node.rand=sample(1:p, p, replace=FALSE)  
        node.index=sort(node.rand,decreasing=F,ind=T)$ix 
       }else{##not shuffle node order
        node.rand=1:p
        node.index=1:p
       }
  
       Y.B=Y[s.pick, node.rand]
       node.B=nodeType[node.rand]
       whiteList.B=whiteList[node.rand, node.rand]
       blackList.B=blackList[node.rand, node.rand]
  
       curRes=hc_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep, restart, i*11+seed, verbose)
       adjRes=curRes$adjacency
       adjRes=adjRes[node.index, node.index]
       adjRes
   }
   plan(sequential)  # reset future plan to sequential 
   gc()              # force garbage collection
   
   ##return results as a p by p by n.boot array  
   return(array(as.numeric(unlist(result)), dim=c(p, p, n.boot)))
}


