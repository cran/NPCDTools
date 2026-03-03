#' @title Plot Diagnostics for GNPC
#'
#' @description
#' This function gives two types of diagnostic plots for the outcomes of the GNPC algorithm. 
#' 
#' The type = "\code{convergence}" option gives two graphs: The upper panel displays the trajectory of the proportion of 
#' membership switches and the lower panel shows the total squared distance along with the iterations. 
#' They illustrate the detailed information about how and whether the 
#' algorithm has converged. 
#' 
#' The type = "\code{individual}" option returns a sequence of squared distances for a 
#' single examinee and the estimates of the examinee's attribute profile 
#' along with the iterations. The plots allow users to investigate how the algorithm
#' arrives at its final classifications. In simulation studies, the true
#' attribute profile can be provided as a reference.
#' 
#' @param x An object of class \code{"GNPC"} when \code{track.convergence = TRUE}.
#' @param type \code{"convergence"} (default) for overall convergence diagnostics,
#'   or \code{"individual"} for a single examinee's estimation trajectory.
#' @param examinee.id An integer indicating which examinee to be plotted. This argument is required 
#'   if \code{"individual"} is specified.
#' @param true.alpha  A numeric vector of length \eqn{K}{} when the true attribute profile
#'   of the examinee is available (optional; usually used for simulation studies).
#' @param top.n.pattern An integer specifying the maximum number of patterns to be displayed.
#'   The default is \code{min(2^K, 6)}. The GNPC estimate and true pattern
#'   (if provided) are always included.
#' @param \dots Additional arguments passed to \code{\link{plot}}.
#'
#'
#' @details
#' For \code{"individual"} plots, the visual elements are:
#' \itemize{
#'   \item \strong{Red line}: the attribute pattern ultimately selected by GNPC.
#'   \item \strong{Black line}: the true attribute pattern (only when
#'     \code{true.alpha} is provided). A small vertical jitter is applied
#'     when it overlaps with the red line.
#'   \item \strong{Other colored lines}: the most competitive candidate
#'     patterns, selected by proximity at the final iteration.
#'   \item \strong{Filled circle at each iteration}: indicates which attribute profile
#'     GNPC assigned to the examinee at each iteration, drawn in the
#'     corresponding line's color. The line for the true attribute profile always has black
#'     circles at every iteration as a fixed reference.
#' }
#'
#' @importFrom graphics layout plot.new
#' @seealso \code{\link{GNPC}}
#'
#' @method plot GNPC
#' @rdname plot.GNPC
#'
#' @examples
#' \dontrun{
#' library(GDINA)
#' set.seed(123)
#' N <- 500
#' Q <- sim30GDINA$simQ
#' gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA")
#' Y <- extract(sim, what = "dat")
#' alpha <- extract(sim, what = "attribute")
#' 
#' # Analyze data using GNPC
#' result <- GNPC(Y, Q, initial.dis = "hamming", initial.gate = "AND")
#'
#' # Convergence
#' plot(result)
#'
#' # Individual with true attribute profile (simulation)
#' plot(result, type = "individual", examinee.id = 1, true.alpha = alpha[1, ])
#'
#' # Individual without true attribute profile (real data)
#' plot(result, type = "individual", examinee.id = 1)
#' }
#' 
#' @importFrom graphics par plot abline grid lines points legend
#' @export
plot.GNPC <- function(x,
                      type = c("convergence", "individual"),
                      examinee.id = NULL,
                      true.alpha = NULL,
                      top.n.pattern = NULL,
                      ...) {
  
  if (!inherits(x, "GNPC")) {
    stop("'x' must be an object of class 'GNPC'.")
  }
  if (is.null(x$convergence)) {
    stop("No convergence information. Run GNPC with track.convergence = TRUE.")
  }
  
  type <- match.arg(type)
  
  # ================================================================
  # Convergence
  # ================================================================
  if (type == "convergence") {
    conv <- x$convergence
    iter <- conv$iteration
    prop.change <- conv$prop.change
    total.dist  <- conv$total.distance
    
    if (is.null(iter) || is.null(prop.change) || is.null(total.dist)) {
      stop("Convergence data is incomplete.")
    }
    
    old.par <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
    on.exit(par(old.par))
    
    plot(iter, prop.change, type = "l",
         xlab = "Iteration", ylab = "Proportion of Changes",
         main = "Proportion of Changed Classifications",
         col = "steelblue", lwd = 2, ...)
    abline(h = 0, lty = 2, col = "gray50")
    grid()
    
    plot(iter, total.dist, type = "l",
         xlab = "Iteration", ylab = "Total Squared Distance",
         main = "Total Distance",
         col = "firebrick", lwd = 2, ...)
    grid()
    
    return(invisible(NULL))
  }
  
  # ================================================================
  # Individual
  # ================================================================
  conv <- x$convergence
  
  if (is.null(conv$ideal.response.history))
    stop("Individual plot requires 'ideal.response.history' in convergence output.")
  if (is.null(conv$classifications))
    stop("Individual plot requires 'classifications' in convergence output.")
  if (is.null(conv$attribute.patterns))
    stop("Individual plot requires 'attribute.patterns' in convergence output.")
  if (is.null(x$Y))
    stop("Individual plot requires response data 'Y' in the GNPC output.")
  
  # Validate examinee.id
  if (is.null(examinee.id))
    stop("'examinee.id' is required when type = 'individual'.")
  if (length(examinee.id) != 1)
    stop("'examinee.id' must be a single integer.")
  
  N <- nrow(x$Y)
  if (examinee.id < 1 || examinee.id > N)
    stop("'examinee.id' must be between 1 and ", N, ".")
  examinee.id <- as.integer(examinee.id)
  
  patterns <- conv$attribute.patterns
  M <- nrow(patterns)
  K <- ncol(patterns)
  n.iter <- conv$n.iter
  
  # Validate true.alpha
  true.idx <- NULL
  if (!is.null(true.alpha)) {
    true.alpha <- as.numeric(true.alpha)
    if (length(true.alpha) != K)
      stop("'true.alpha' must have length ", K, ".")
    for (m in 1:M) {
      if (all(patterns[m, ] == true.alpha)) {
        true.idx <- m
        break
      }
    }
    if (is.null(true.idx))
      warning("'true.alpha' does not match any attribute pattern. Ignoring.")
  }
  
  # Default top.n.pattern
  if (is.null(top.n.pattern)) {
    top.n.pattern <- min(M, 6)
  } else {
    if (!is.numeric(top.n.pattern) || length(top.n.pattern) != 1 || top.n.pattern < 2)
      stop("'top.n.pattern' must be an integer >= 2.")
    top.n.pattern <- min(as.integer(top.n.pattern), M)
  }
  
  # Compute distances at each iteration
  Y.i <- x$Y[examinee.id, ]
  distances <- matrix(NA, nrow = n.iter, ncol = M)
  for (t in 1:n.iter) {
    Ideal.w.t <- conv$ideal.response.history[[t]]
    distances[t, ] <- colSums((t(Ideal.w.t) - Y.i)^2)
  }
  
  # Select patterns to display
  gnpc.idx <- conv$classifications[n.iter, examinee.id]
  assigned.classes <- conv$classifications[, examinee.id]
  
  final.dists <- distances[n.iter, ]
  ranked <- order(final.dists)
  show.patterns <- ranked[1:min(top.n.pattern, M)]
  
  # Guarantee GNPC final estimate
  if (!(gnpc.idx %in% show.patterns))
    show.patterns <- c(show.patterns, gnpc.idx)
  # Guarantee true pattern
  if (!is.null(true.idx) && !(true.idx %in% show.patterns))
    show.patterns <- c(show.patterns, true.idx)
  # Guarantee all ever-assigned patterns
  ever.assigned <- unique(assigned.classes)
  for (ea in ever.assigned) {
    if (!(ea %in% show.patterns))
      show.patterns <- c(show.patterns, ea)
  }
  
  show.patterns <- sort(unique(show.patterns))
  n.show <- length(show.patterns)
  
  # Labels
  class.labels <- apply(patterns, 1, function(p)
    paste0("(", paste(p, collapse = ","), ")"))
  
  # ---- Color assignment ----
  other.palette <- c("#377EB8", "#4DAF4A", "#984EA3",
                     "#FF7F00", "#A65628", "#F781BF", "#66C2A5",
                     "#8DD3C7", "#BEBADA", "#FB8072")
  colors <- character(n.show)
  other.idx <- 1
  for (i in seq_along(show.patterns)) {
    m <- show.patterns[i]
    if (m == gnpc.idx) {
      colors[i] <- "red"
    } else if (!is.null(true.idx) && m == true.idx) {
      colors[i] <- "black"
    } else {
      colors[i] <- other.palette[((other.idx - 1) %% length(other.palette)) + 1]
      other.idx <- other.idx + 1
    }
  }
  
  gnpc.true.overlap <- !is.null(true.idx) && gnpc.idx == true.idx
  
  # ---- Plot setup ----
  y.range <- range(distances[, show.patterns], na.rm = TRUE)
  y.pad <- diff(y.range) * 0.1
  ylim <- c(y.range[1] - y.pad, y.range[2] + y.pad)
  
  # Jitter only when GNPC and TRUE overlap
  true.jitter <- 0
  if (gnpc.true.overlap) {
    true.jitter <- diff(ylim) * 0.008
  }
  
  main.title <- paste0("Examinee ", examinee.id)
  if (!is.null(true.idx)) {
    main.title <- paste0(main.title,
                         if (gnpc.idx == true.idx) " (correct)" else " (incorrect)")
  }
  
  # ---- Layout: plot area + legend area ----
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1.2))
  
  # Panel 1: main plot
  par(mar = c(4, 4, 3, 1))
  
  plot(1:n.iter, distances[, show.patterns[1]], type = "n",
       xlab = "Iteration", ylab = "Squared Distance",
       main = main.title,
       xlim = c(1, n.iter), ylim = ylim, ...)
  grid()
  
  # Draw lines
  for (i in seq_along(show.patterns)) {
    m <- show.patterns[i]
    lines(1:n.iter, distances[, m], col = colors[i], lwd = 1.5)
  }
  
  # Draw TRUE line on top in black with jitter (when overlapping)
  if (gnpc.true.overlap) {
    lines(1:n.iter, distances[, true.idx] + true.jitter,
          col = "black", lwd = 1.5)
  }
  
  
  # Filled circles: GNPC choice at each iteration
  for (t in 1:n.iter) {
    assigned.m <- assigned.classes[t]
    if (assigned.m %in% show.patterns) {
      color.idx <- which(show.patterns == assigned.m)
      points(t, distances[t, assigned.m],
             pch = 19, col = colors[color.idx], cex = 1.2)
    }
  }
  
  # Filled circles: TRUE pattern (every iteration, with jitter)
  if (!is.null(true.idx) && true.idx %in% show.patterns) {
    for (t in 1:n.iter) {
      points(t, distances[t, true.idx] + true.jitter,
             pch = 15, col = "black", cex = 1.2)
    }
  }
  
  # ---- Panel 2: legend ----
  par(mar = c(4, 0, 3, 1))
  plot.new()
  
  # Build legend entries
  legend.labels <- character(n.show)
  legend.lwd   <- numeric(n.show)
  legend.pch   <- rep(NA_integer_, n.show)
  legend.col   <- colors
  
  for (i in seq_along(show.patterns)) {
    m <- show.patterns[i]
    label <- class.labels[m]
    
    if (m == gnpc.idx) {
      label <- paste0(label, "\n[GNPC]")
      legend.pch[i] <- 19
    } else if (!is.null(true.idx) && m == true.idx) {
      label <- paste0(label, "\n[TRUE]")
      legend.pch[i] <- 15
    }
    
    legend.labels[i] <- label
    legend.lwd[i] <- 1.5
  }
  
  # When GNPC and TRUE overlap, add separate black TRUE entry
  if (gnpc.true.overlap) {
    legend.labels <- c(legend.labels, paste0(class.labels[true.idx], "\n[TRUE]"))
    legend.col   <- c(legend.col, "black")
    legend.lwd   <- c(legend.lwd, 1.5)
    legend.pch   <- c(legend.pch, 15)
  }
  
  legend("center",
         title = "Attribute Pattern",
         legend = legend.labels,
         col = legend.col, lwd = legend.lwd, pch = legend.pch,
         bty = "o", bg = "white",
         cex = 0.9, pt.cex = 0.8, y.intersp = 1.5)
  
  invisible(NULL)
}