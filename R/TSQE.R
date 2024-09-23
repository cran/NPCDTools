#' @title Two-step Q-matrix Estimation Method
#'
#' @description The function is used to estimate the Q-matrix based on the
#'     data (responses) using the two-step Q-matrix estimation method.
#'
#' @section Estimation Method:
#'
#' The TSQE method merges the Provisional Attribute Extraction (PAE) algorithm
#' with a Q-matrix refinement-and-validation method
#' including the Q-Matrix Refinement (QR) Method and the G-DINA Model
#' Discrimination Index (GDI). Specifically, the PAE algorithm relies on
#' classic exploratory factor analysis (EFA) combined with a unique stopping
#' rule for identifying a provisional Q-matrix, and the resulting provisional
#' Q-Matrix will be "polished" with a refinement method to derive the final
#' estimation of Q-matrix.
#'
#' @section The Provisional Attribute Extraction (PAE) Algorithm:
#'
#' The initial step of the algorithm is to aggregating the collected Q-Matrix
#' into an inter-item tetrachoric correlation matrix. The reason for using
#' tetrachoric correlation is that the examinee responses are binary, so it
#' is more appropriate than the Pearson product moment correlation coefficient.
#' See Chiu et al. (2022) for details. The next step is to use factor analysis
#' on the item-correlation matrix, and treat the extracted factors as proxies
#' for the latent attributes. The third step concerns identifying which specific
#' attributes are required for which item:
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
#'@section The Q-Matrix Refinment (QR) Method:
#'
#' This function implements the Q-matrix refinement method developed by
#' Chiu (2013), which is also based on the aforementioned nonparametric
#' classification methods (Chiu & Douglas, 2013). This Q-matrix refinement
#' method corrects potential misspecified entries of the Q-matrix through
#' comparisons of the residual sum of squares computed from the observed
#' and the ideal item responses.
#'
#' The algorithm operates by minimizing the RSS. Recall that \eqn{Y_{ij}}
#' is the observed response and \eqn{\eta_{ij}} is the ideal response.
#' Then the RSS of item \eqn{j} for examinee \eqn{i} is defined as
#' \deqn{RSS_{ij} = (Y_{ij} - \eta_{ij})^2}.
#' The RSS of item \eqn{j} across all examinees is therefor
#' \deqn{RSS_{j} = \sum_{i=1}^{N} (Y_{ij} - \eta_{ij})^2 = \sum_{m=1}^{2^k} \sum_{i \in C_{m}} (Y_{ij} - \eta_{jm})^2}
#' where \eqn{C_m} is the latent proficiency-class \eqn{m},
#' and \eqn{N} is the number of examinees.
#' Chiu(2013) proved that the expectation of \eqn{RSS_j} is minimized for
#' the correct q-vector among the \eqn{2^K - 1} candidates. Please see  the
#' paper for the justification.
#'
#'@section The G-DINA Model Discrimination Index (GDI):
#'
#' The GDI is an extension of de la Torre's (2008) \eqn{\delta}-method,
#' which has a limitation that it cannot be used with CDMs that
#' devide examinees into more than two groups. In response to the limitation,
#' de la Torre and Chiu (2016) porposed to select that item attribute vector
#' which maximizes the weighted variance of the probabilities of a correct
#' response for the different groups defined as
#' \deqn{\zeta_j^2 = \sum_{l=1}^{2^{K_j}} P(\alpha_{lj}) \left[ P(Y_{ij} = 1 \mid \alpha_{lj}) - \bar{P}_{j} \right]^2}
#' where \eqn{P(\alpha_{lj})} is the posterior probability for the proficiency class
#' \eqn{\alpha_{lj}}, and \eqn{\bar{P}_{j} =  \sum_{l=1}^{2^{K_j}} P(\alpha_{lj})P(Y_{ij} = 1 \mid \alpha_{lj})},
#' where \eqn{l = 1,2,...,2^{K_j}}. De la Torre and Chiu (2016) called \eqn{\zeta^2}
#' the GDI, which can be applied to any CDM that can be reparameterized in
#' terms of the G-DINA model.
#'
#' @param Y A \eqn{N \times J} binary data matrix consisting of responses
#'     from \eqn{N} examinees to \eqn{J} items
#' @param K The number of attributes in the Q-matrix
#' @param input.cor The type of correlation used to compute the input for the
#'     exploratory factor analysis. It could be the tetrachoric or Pearson correlation.
#' @param ref.method The refinement method use to polish the provisional
#'     Q-matrix obtained from the EFA. Currently available methods include
#'     the Q-matrix refinement (QR) method and the G-DINA discrimination index (GDI).
#' @param GDI.model The CDM used in the GDI algorithm to fit the data. Currently
#'     available models include the DINA model, the ACDM, the RRUM, and the
#'     G-DINA model
#' @param cutoff The cutoff used to dichotomize the entries in the provisional
#'    Q-matrix
#'
#' @return The function returns the estimated Q-matrix
#' @references
#' Chiu, C. Y. (2013). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Applied Psychological Measurement, 37(8)}, 598-618.
#'
#' Chiu, C. Y., & Douglas, J. A. (2013). A nonparametric approach to cognitive diagnosis by proximity to ideal response patterns. \emph{Journal of Classification 30(2)}, 225-250.
#'
#' de la Torre, J., & Chiu, C.-Y. (2016) A general method of empirical Q-matrix validation. \emph{Psychometrika, 81}, 253-73.
#'
#' de la Torre, J. (2008). An empirically based method of Q-matrix validation for the DINA model: Development and applications. \emph{Journal of Educational Measurement, 45}, 343-362.
#'
#' @export
#'
#'

TSQE = function(Y, K, input.cor = c("tetrachoric", "Pearson"), ref.method = c("QR", "GDI"), GDI.model = c("DINA", "ACDM", "RRUM", "GDINA"), cutoff = 0.8){
  N = nrow(Y)
  J = ncol(Y)



  if (input.cor == "Pearson"){
    cor.data = stats::cor(Y)
    fa.out = stats::factanal(factors = K, covmat = cor.data)
  }else if (input.cor == "tetrachoric"){
    cor.data = suppressWarnings(psych::tetrachoric(Y,smooth = T)$rho)
    cor.data0 = cor.data
    cor.data = suppressWarnings(psych::cor.smooth(round(cor.data0, 5)))

    #=== If tetrachoric correlation doesn't work, switch to Pearson ======
    fa.out = try(stats::factanal(factors = K, covmat = cor.data), silent = T)
    if (inherits(fa.out, "try-error")){
      warning("The tetrachoric correlation doesn't work, the correlation has switched to Pearson.", call. = FALSE)
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
    ref.Q = SimDesign::quiet(NPCD::Qrefine(Y, fQ, gate = "AND", max.ite = 50)$modified.Q)
  } else if (ref.method == "GDI"){
    mod = GDINA::GDINA(dat = Y, Q = fQ, model = GDI.model, verbose = 0)
    ref.Q = GDINA::Qval(mod, method = "PVAF",eps = 0.95)$sug.Q
  }
  return(ref.Q)
}
