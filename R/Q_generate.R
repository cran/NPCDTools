#' @title Generation of Dichotomous Q-Matrix
#'
#' @description The function generates a complete Q-matrix based on a
#'     pre-specified probability of getting a one.
#'
#' @param K The number of attributes
#' @param J The number of items
#' @param p The probability of getting a one in the Q-matrix
#' @param single.att Whether all the single attribute patterns are included.
#' If \code{T}, the completeness of the Q-matrix is guaranteed.
#'
#' @return The function returns a complete dichotomous Q-matrix
#'
#' @export
#' @examples
#' q = Q.generate(3,20,0.5,single.att = TRUE)
#' q1 = Q.generate(5,30,0.6,single.att = FALSE)
#'
Q.generate = function(K, J, p, single.att = TRUE){
  if (K>J){
    stop("The number of items should be larger than the number of attributes.")
  }

  if (single.att == TRUE){
    J_len = ceiling((K-1)/(1-p))
    if (K==J){
      if (J < J_len){
        warning(paste("The number of items should be at least",J_len,"to get the input p."))
      }
      Q = diag(K)
    }

    if (K < J){
      Q0 = diag(K)
      n=1

      if(J<J_len){
        m=J-K
        warning(paste("The number of items should be at least",J_len,"to get the input p."))
        while (n > 0){
          Q1 = matrix(sample(c(0,1), m*K, replace=T, prob = c(1-p, p)), m, K)
          n = length(which(rowSums(Q1)==0))
        }
        Q = rbind(Q0, Q1)
      }

      else if(J == J_len){
        m=J-K
        warning(paste("There will be",m,"items having all attributes."))
        Q1 = matrix(sample(1, m*K, replace=T), m, K)
        Q = rbind(Q0,Q1)
      }

      else if(J_len < J && J<=(2*J_len)){
        m=J-K
        p.new = (p*J-1)/(J-K)
        while (n > 0){
          Q1 = matrix(sample(c(0,1), m*K, replace=T, prob = c(1-p.new, p.new)), m, K)
          n = length(which(rowSums(Q1)==0))
        }
        Q = rbind(Q0,Q1)
      }

      else{
        num_part = ceiling(J/(2*J_len))-2
        p.new = (p*2*J_len-1)/(2*J_len-K)
        m=2*J_len-K

        #get the first small matrix with p.new
        while (n > 0){
          Q1 = matrix(sample(c(0,1), m*K, replace=T, prob = c(1-p.new, p.new)), m, K)
          n = length(which(rowSums(Q1)==0))
        }

        #get the last small matrix with p of which length is not m.
        n.tail = 1
        m.tail = J%%(2*J_len)

        if(m.tail==0){
          m.tail = 2*J_len
        }

        while (n.tail > 0){
          Q2 = matrix(sample(c(0,1), m.tail*K, replace=T, prob = c(1-p, p)), m.tail, K)
          n.tail = length(which(rowSums(Q2)==0))
        }

        Q = rbind(Q0,Q1,Q2)

        while (num_part>0) {
          n.rest = 1
          while (n.rest >0){
            Q3 = matrix(sample(c(0,1), (m+K)*K, replace=T, prob = c(1-p, p)), m+K, K)
            n.rest = length(which(rowSums(Q3)==0))
          }
          Q = rbind(Q,Q3)
          num_part = num_part-1
        }
      }
    }
  }
  else{
    n=1
    while (n > 0){
      Q = matrix(sample(c(0,1), J*K, replace=T, prob = c(1-p, p)), J, K)
      n = length(which(rowSums(Q)==0))
    }
  }
  return(Q)
}
