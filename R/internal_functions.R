#=====Estimation of Proficiency Class via Hamming distance=====
hamming = function(Ideal,Y){
  M=nrow(Ideal)
  J=ncol(Ideal)
  N=nrow(Y)

  H.class=NULL
  h.ntie=0
  for (i in 1:N)
  {
    ham=apply(abs(matrix(rep(Y[i,], M), M, J, byrow=TRUE)-Ideal), 1, sum)
    min.ham=which(ham==min(ham))
    if (length(min.ham)!=1)
    {
      h.ntie=h.ntie+1
      min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties
    }else h.ntie=h.ntie
    H.class=c(H.class, min.ham)
  }
  return(H.class)
}

#=====Estimation of Proficiency Class via weighted Hamming distance=====
whamming=function(Ideal,Y){
  M=nrow(Ideal)
  J=ncol(Ideal)
  N=nrow(Y)

  p.bar=apply(Y,2,mean)
  weight=1/(p.bar*(1-p.bar))
  WH.class=NULL
  ntie=0
  H=NULL
  for (i in 1:N)
  {
    ham=apply(matrix(rep(weight, M), M, J, byrow=TRUE)*abs(matrix(rep(Y[i,], M), M, J, byrow=TRUE)-Ideal), 1, sum)
    H=rbind(H,ham)
    min.ham=which(ham==min(ham))
    if (length(min.ham)!=1)
    {
      ntie=ntie+1
      min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties
    }else ntie=ntie
    WH.class=c(WH.class, min.ham)
  }
  return(WH.class)
}

class.generate = function(K){
  M <- diag(K)
  for (i in 2:K){
    M <- rbind(M,t(apply(utils::combn(K,i),2,function(x){apply(M[x,],2,sum)})))
  }
  M <- rbind(0,M)
  return(M)
}

epc.generate = function(mcq,O,J,K,key,Class){
  Q <- mcq[, -c(1:2)]
  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  no.options <- rep(O, J)

  # "scored" option
  eta.class <- matrix(0,J,2^K)
  for(j in 1:J){

    Qj <- Q[which(item.no==j),,drop=FALSE]  # won't change data type
    Kj <- rowSums(Qj)
    kj.order <- order(Kj,decreasing = TRUE)  #get location;
    coded.op.j <- coded.op[which(item.no==j)]

    if (num.coded[j]>1){

      key.loc <- which(coded.op.j==key[j])
      if(key.loc!=kj.order[1]){
        kj.order <- c(key.loc,setdiff(kj.order,key.loc))  # make the answer goes first
      }


      Qj <- Qj[kj.order,]
      #et.label[[j]] <- c(apply(Qj,1,paste,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj))))))
      #matrix(unlist(lapply(apply(B,2,function(x){which(x==max(x))}),function(y){if (length(y)<O){c(y,rep(0,O-length(y)))}else{y=y}})),nrow=2^K,byrow=TRUE)

      #tmp <- eta(Qj)
      eta.j <- apply(B,2,which.max)
      max.j <- max(eta.j)
      eta.j <- eta.j - 1
      eta.j[eta.j==0] <- max.j
      eta.j <- num.coded[j]+1-eta.j


    }else{

      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj))))))

      eta.j <- apply(B,2,which.max)
      eta.j <- eta.j-1
    }

    eta.class[j,] <- eta.j

  }
  return(eta.class)
}

score.option = function(mcq,O){


  J = max(mcq[,1])
  Q <- mcq[, -c(1:2)]

  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m

  op <- vector(mode="list",length=J)

  mc.q = as.data.frame(mcq)
  save.m <- c()
  key <- c()
  for (i in 1:J){
    key[i] <- mc.q$Option[min(which(mc.q$Item==i))]
    save.m <- c(save.m,length(which(mc.q$Item==i)))
  }


  for (i in 1:J){

    Qj <- Q[which(item.no==i),,drop=FALSE]  # won't change data type
    Kj <- rowSums(Qj)
    kj.order <- order(Kj,decreasing = TRUE)  #get location;
    coded.op.j <- coded.op[which(item.no==i)]

    sub.info <- Item.info[which(Item.info[,1]==i),,drop=FALSE]

    if (num.coded[i]>1){

      key.loc <- which(coded.op.j==key[i])
      if(key.loc!=kj.order[1]){
        kj.order <- c(key.loc,setdiff(kj.order,key.loc))  # make the answer goes first
      }

      if (num.coded[i]==O){
        g = c(O:1)
      } else {g = c(c(nrow(sub.info):1),rep(0,O-nrow(sub.info)))}


      ans = c(sub.info[,2][kj.order],setdiff(c(1:O),sub.info[,2][kj.order]))
      score <- cbind(ans,g)


    } else {
      ans =c(unname(sub.info[,2]),setdiff(c(1:O),sub.info[,2]))
      g = c(1,rep(0,(O-1)))
      score <- cbind(ans,g)
    }

    op[[i]] <- score

  }

  return(op)
}

