#' @title Two-step Q-matrix estimation method
#'
#' @description The function estimates the Q-matrix based on the
#'     response data using the two-step Q-matrix estimation method.
#'
#' @section Details:
#'
#' The TSQE method estimates a Q-matrix by integrating the provisional attribute extraction (PAE) algorithm
#' with a Q-matrix refinement-and-validation method,
#' such as the Q-Matrix Refinement (QR) method and the G-DINA Model
#' Discrimination Index (GDI). Specifically, the PAE algorithm relies on
#' classic exploratory factor analysis (EFA) combined with a unique stopping
#' rule for identifying a provisional Q-matrix, and the resulting provisional
#' Q-Matrix is "polished" by a refinement method to derive the finalized
#' estimation of Q-matrix.
#'
#'
#' The PAE Algorithm starts with computing the 
#' inter-item tetrachoric correlation matrix. The reason for using
#' tetrachoric correlation is that the examinee responses are binary, so it
#' is more appropriate than the Pearson product moment correlation coefficient.
#' See \enc{Köhn}{Koehn} et al. (2025) for details. The next step is to use factor analysis
#' on the item-correlation matrix, and treat the extracted factors as proxies
#' for the latent attributes. The third step concerns the identification of specific
#' attributes required for each item. The detailed algorithm is described below:
#'
#' \describe{
#'   \item{(1)}{Initialize the item index as \eqn{j = 1}.}
#'   \item{(2)}{Let \eqn{l_{jk}} denote the loading of item \eqn{j} on factor \eqn{k}, where \eqn{k = 1,2,...,K}.}
#'   \item{(3)}{Arrange the loadings in descending order. Define a mapping
#'       function \eqn{f(k) = t}, where \eqn{t} is the order index.
#'       Hence, \eqn{l_{j(1)}} will indicate the maximum loading,
#'       while \eqn{l_{j(K)}} will indicate the minimum loading.}
#'   \item{(4)}{Define \deqn{p_j(t) = \frac{\sum_{h=1}^t l_{j(h)}^2}{\sum_{k=1}^K l_{jk}^2}}
#'       as the proportion of the communality of item \eqn{j} accounted for
#'       by the first \eqn{t} factors.}
#'   \item{(5)}{Define \deqn{K_j = \min \{ t \mid p_j(t) \geq \lambda \}},
#'       where \eqn{\lambda} is the cut-off value for the desired proportion
#'       of item variance-accounted-for. Then, the ordered entries of the
#'       provisional q-vector of item \eqn{j} are obtained as
#'       \deqn{q_{j(t)}^* = \begin{cases}
#'       1 & \text{if } t \leq K_j \\
#'       0 & \text{if } t > K_j
#'       \end{cases}}.}
#'   \item{(6)}{Identify \eqn{q_j^* = (q_{j1}^*,q_{j2}^*,...,q_{jK}^*)}
#'       by rearranging the ordered entries of the q-vector using the inverse function \eqn{k = f^{-1}(t)}.}
#'   \item{(7)}{Set \eqn{j = j + 1} and repeat (2) to (6) until \eqn{j = J}.
#'       Then denote the provisional Q-matrix as \eqn{\mathbf{Q}^*}.}
#' }
#'
#' The provisional Q-matrix \eqn{\mathbf{Q}^*}{} is then refined by 
#' using either the \code{QR} or \code{GDI} method.
#'
#'
#' @param Y A \eqn{N \times J} binary data matrix consisting of responses
#'     from \eqn{N} examinees to \eqn{J} items
#' @param K The number of attributes in the Q-matrix
#' @param input.cor The type of correlation used as input for the
#'     provisional attribute extraction (PAE) algorithm. It could be the 
#'     \code{tetrachoric} or \code{pearson} correlation.
#' @param ref.method The refinement method used to polish the provisional
#'     Q-matrix obtained from the PAE. Currently available methods include
#'     the Q-matrix refinement (\code{QR}) method and the G-DINA discrimination index (\code{GDI}).
#' @param GDI.model The CDM used in the GDI algorithm to fit the data. Currently
#'     available models include the DINA model, the ACDM, the RRUM, and the
#'     G-DINA model.
#' @param cutoff The cutoff used to dichotomize the entries in the provisional
#'    Q-matrix. The default is 0.8.
#'
#' @return The function returns the estimated Q-matrix.
#' @references
#' Chiu, C. Y. (2013). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. 
#' \emph{Applied Psychological Measurement, 37(8)}, 598-618.
#' \doi{10.1177/0146621613488436}
#'
#' de la Torre, J., & Chiu, C.-Y. (2016). A general method of empirical Q-matrix validation. 
#' \emph{Psychometrika, 81}, 253-73.
#' \doi{10.1007/s11336-015-9467-8}
#' 
#' Köhn, H. F., Chiu, C.-Y., Oluwalana, O., Kim, H. & Wang, J. (2025). A two-step Q-matrix estimation 
#' method, \emph{Applied Psychological Measurement, 49}(1-2), 3-28.
#' \doi{10.1177/01466216241284418}
#'
#' @export
#' 
#' @seealso \code{\link{QR}}
#' 
#' @examples
#' \dontrun{
#' library(GDINA)
#' N = 1000
#' Q = sim30GDINA$simQ
#' J = nrow(Q)
#' K= ncol(Q)
#' gs = data.frame(guess=rep(0.2,J),slip=rep(0.2,J))
#' sim = simGDINA(N,Q,gs.parm = gs,model = "DINA")
#' Y = extract(sim,what = "dat")
#' 
#' ## Run TSQE method with QR
#' est.Q = TSQE(Y, K, input.cor = "tetrachoric", ref.method = "QR", cutoff = 0.8)
#' 
#' ## If the recovery rate is to be computed, the columns of the estimated Q-matrix 
#' ## should be permuted so that they align with those of the true Q-matrix. 
#' best.est.Q = bestQperm(est.Q, Q)
#' 
#' ## Compute the recovery rate
#' RR(best.est.Q, Q)
#' }

TSQE = function(Y, K, input.cor = c("tetrachoric", "pearson"), 
                      ref.method = c("QR", "GDI"), 
                      GDI.model = c("GDINA","DINA", "ACDM", "RRUM"), 
                      cutoff = 0.8){
  if (missing(ref.method)) {
    message("Note: No 'ref.method' specified. Defaulting to 'QR'.")
    message("      To use the GDI method, please specify ref.method = 'GDI'.\n")
  }
  ref.method <- match.arg(toupper(ref.method), c("QR", "GDI"))
  if (ref.method == "GDI") {
    if (missing(GDI.model)) {
      default_model <- GDI.model[1] 
      message("Note: 'GDI.model' not specified. Using default model: ", default_model, ".")
    }
    GDI.model <- match.arg(toupper(GDI.model), c("GDINA", "DINA", "ACDM", "RRUM"))
  }
  #  input.cor: Convert to lower case, match against valid options
  input.cor <- match.arg(tolower(input.cor), c("tetrachoric", "pearson"))
  
  N = nrow(Y)
  J = ncol(Y)

  if (input.cor == "pearson"){
    cor.data = stats::cor(Y)
    fa.out = stats::factanal(factors = K, covmat = cor.data)
  }else if (input.cor == "tetrachoric"){
    cor.data = suppressWarnings(psych::tetrachoric(Y,smooth = T)$rho)
    cor.data0 = cor.data
    cor.data = suppressWarnings(psych::cor.smooth(round(cor.data0, 5)))

    #=== If tetrachoric correlation doesn't work, switch to pearson ======
    fa.out = try(stats::factanal(factors = K, covmat = cor.data), silent = T)
    if (inherits(fa.out, "try-error")){
      warning("The tetrachoric correlation doesn't work, the correlation has switched to pearson.", call. = FALSE)
      cor.data = stats::cor(Y)
      fa.out = stats::factanal(factors = K, covmat = cor.data)
    }
  }

  #== Extract loadings and Factors ========
  loads = fa.out$loadings[,1:K]

  #==== Determine the provisional Q based on the given cutoff =======
  fQ = NULL
  for (j in 1:J){
    l = abs(loads[j,])
    ordered.l = l[order(l, decreasing = T)]
    ordered.r = ordered.l^2/sum(ordered.l^2)
    c.sum = cumsum(ordered.r)
    cut = (c.sum <= cutoff)*1
    cut[sum(cut)+1]=1
    if (sum(cut)==0){cut[1]=1}
    fq = cut[K-rank(l)+1]
    fQ = rbind(fQ, fq)
  }

  #====== Refine the provisional Q ============
  if (ref.method == "QR"){
    ref.Q = suppressMessages(QR(Y, fQ, gate = "AND", max.ite = 50)$modified.Q)
  } else if (ref.method == "GDI"){
    mod = GDINA::GDINA(dat = Y, Q = fQ, model = GDI.model, verbose = 0)
    ref.Q = GDINA::Qval(mod, method = "PVAF",eps = 0.95)$sug.Q
  }
  return(ref.Q)
}
