#' @title Estimation of examinees' attribute profiles using the GNPC method
#'
#' @description
#' Function \code{GNPC} is used to estimate examinees' attribute profiles using
#'     the general nonparametric classification (GNPC) method
#'     (Chiu, Sun, & Bian, 2018; Chiu & Koehn, 2019). It can be
#'     used with data conforming to any CDMs.
#'
#' @section GNPC algorithm:
#' A weighted ideal response \eqn{\eta^{(w)}}, defined as the convex combination
#' of \eqn{\eta^(c)} and \eqn{\eta^(d)}, is proposed.
#' Suppose item j requires \eqn{K_{j}^* \leq {K}} attributes that, without loss of
#' generality, have been permuted to the first \eqn{K_{j}^*} positions of the item
#' attribute vector \eqn{\boldsymbol{q_j}}. For each item j and \eqn{\mathcal{C}_{l}},
#' the weighted ideal response \eqn{\eta_{ij}^{(w)}} is defined as the convex combination
#' \eqn{\eta_{ij}^{(w)} = w _{lj} \eta_{lj}^{(c)}+(1-w_{lj})\eta_{lj}^{(d)}}
#' where \eqn{0\leq w_{lj}\leq 1}. The distance between the observed responses
#' to item j and the weighted ideal responses \eqn{w_{lj}^{(w)}} of examinees
#' in \eqn{\mathcal{C}_{l}} is defined as the sum of squared deviations:
#' \eqn{d_{lj} = \sum_{i \in \mathcal {C}_{l}} (y_{ij} - \eta_{lj}^{(w)})^2=\sum_{i \in \mathcal {C}_{l}}(y_{ij}-w_{lj}\eta_{lj}^{(c)}-(1-w_{lj})\eta_{lj}^{(d)})}
#' Thus, \eqn{\widehat{w_{lj}}} can be minimizing \eqn{d_{lj}}:
#' \eqn{\widehat{w_{lj}}=\frac{\sum_{i \in \mathcal {C}_{l}}(y_{ij}-\eta_{lj}^{(d)})}{\left \| \mathcal{C}_{l} \right \|(\eta_{lj}^{(c)}-\eta_{lj}^{(d)})}}
#'
#' As a viable alternative to \eqn{\boldsymbol{\eta^{(c)}}} for obtaining initial
#' estimates of the proficiency classes, Chiu et al. (2018) suggested to
#' use an ideal response with fixed weights defined as
#' \eqn{\eta_{lj}^{(fw)}=\frac{\sum_{k=1}^{K}\alpha_{k}q_{jk}}{K}\eta_{lj}^{(c)}+(1-\frac{\sum_{k=1}^{K}\alpha_{k}q_{jk}}{K})\eta_{lj}^{(d)}}
#'
#' @param Y A \eqn{N \times J} binary data matrix consisting of the responses
#'     from \eqn{N} examinees to \eqn{J} items.
#' @param Q A \eqn{J \times K} binary Q-matrix where the entry \eqn{q_{jk}}
#'     describing whether the \eqn{k}th attribute is required by the \eqn{j}th item.
#' @param initial.dis The type of distance used in the \code{AlphaNP} to carry
#'     out the initial attribute profiles for the GNPC method.
#'     Allowable options are \code{"hamming"} and \code{"whamming"} representing
#'     the Hamming and the weighted Hamming distances, respectively.
#' @param initial.gate The type of relation between examinees' attribute profiles
#'     and the items.
#'     Allowable relations are \code{"AND"}, \code{"OR"}, and \code{"Mix"},
#'     representing the conjunctive, disjunctive, and mixed relations, respectively
#'
#' @return The function returns a series of outputs, including
#' \describe{
#' \item{att.est}{The estimates of examinees' attribute profiles}
#' \item{class}{The estimates of examinees' class memberships}
#' \item{weighted.ideal}{The weighted ideal responses}
#' \item{weight}{The weights used to compute the weighted ideal responses}
#' }
#'
#' @export



GNPC=function(Y, Q, initial.dis= c("hamming", "whamming"), initial.gate = c("AND", "OR", "Mix")){

  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)
  if (!is.null(check)) return(warning(check))

  N=dim(Y)[1]
  K=dim(Q)[2]
  J=dim(Q)[1]
  M=2^K

  pattern <- diag(K)
  for (l in 2:K){
    pattern <- rbind(pattern,t(apply(utils::combn(K,l),2,function(x){apply(pattern[x,],2,sum)})))
  }
  pattern <- rbind(0,pattern)

  #============================
  # Conjunctive Ideal Response
  #============================
  Ideal=pattern%*%t(Q) #M*K %*% K*J = M * J
  Ideal.conj=1*(Ideal==(matrix(1,M)%*%t(rowSums(Q))))

  #============================
  # Disjunctive Ideal Response
  #============================
  Ideal.dis=1*(Ideal>=1)

  #=================================
  # Assigning Initial Weight
  # Initial Weighted Ideal Response
  #=================================
  weight=Ideal/matrix(rep(colSums(t(Q)),M),M,J,T)
  Ideal.mix=Ideal.conj+(Ideal.dis-Ideal.conj)*weight

  #===========================================
  # Classify Initial Latent Class
  #===========================================
  if(initial.gate=="AND") {
    Ideal=Ideal.conj}
  else if(initial.gate=="OR") {
      Ideal=Ideal.dis}
  else if(initial.gate=="Mix") {
        Ideal=Ideal.mix}
  else{
    return(warning("Gate specification not valid."))
  }

  if (initial.dis=="hamming") {
    initial.class = hamming(Ideal,Y)
  }
  else if (initial.dis=="whamming") {
    initial.class = whamming(Ideal,Y)
  }
  else{
    return(warning("Initial distance method specification not valid."))
  }

  #============================================
  #  Iteration Starts
  #============================================
  d=1
  time=0
  while(d>0.001)
  { #print(time)
    time=time+1
    #===================================================
    # Compute the general weights using the closed form
    #===================================================
    w=NULL
    temp=matrix(NA,M,1)
    Ideal.comb=Ideal.dis+Ideal.conj
    for (j in 1:J) {
      pQ=pattern*matrix(Q[j,]==1,M,K,T)
      upQ=unique(pQ)
      ng=nrow(upQ)
      for (g in 1:ng){
        match <- apply(pQ, 1, identical, upQ[g,])
        m=which(match)
        c= which(initial.class %in% m)
        c=c[!is.na(c)]
        if (length(c)==0){temp[m,]=0.01}
        else if (length(c)!=0){
          temp[m,]=sum(Y[c,j])/sum((1-Ideal.conj[initial.class[c],j])^2)}
      }
      temp[which(Ideal.comb[,j]==0)]=.01
      temp[which(Ideal.comb[,j]==2)]=.99
      w=cbind(w,temp)
    }
    Ideal.w=(1-w)*Ideal.conj+w*Ideal.dis
    w.class=NULL
    for (i in 1:N){
      ham=rowSums(((matrix(Y[i,], M, J, byrow=TRUE)-Ideal.w))^2)
      min.ham=which(ham==min(ham))
      if (length(min.ham)!=1) min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties

      w.class=c(w.class, min.ham)
    }
    d=length(which(w.class-initial.class!=0))/N
    initial.class=w.class
  }
  att.est=pattern[w.class,]
  output = list(att.est = att.est, class = w.class, weighted.ideal = Ideal.w, weight = w)
  #class(output) = "GNPC"
  return(output)
}
