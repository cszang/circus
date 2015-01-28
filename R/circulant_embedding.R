##' Exact simulation by circulant embedding
##'
##' Circulant embedding is used to simulate the time series data n
##' times as time series with the same frequency characteristics like
##' the original time series (Percival & Constantine, 2006). The
##' resulting simulated time series are scaled to the mean and
##' standard deviation of the original time series.
##' @param x original time series
##' @param n number of simulations (must be a multiple of 2)
##' @return a data.frame with n simulated times series in columns.
##' @examples
##' x <- rnorm(100)
##' xsim <- circ_embed(x)
##' @export
circ_embed <- function(x, n = 1000) {

  if (!is.null(dim(x))) {
    stop("x must be a vector.")
  }
  
  if (length(x) < 31) {
    stop("x must have a length of at least 32.")
  }

  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }

  if (n %% 1 != 0) {
    stop("n must be an integer value.")
  }

  if (n %% 2 != 0) {
    stop("n must be a multiple of 2.")
  }
  
  nh <- n/2
  
  ## subtract series mean
  trm <- x - mean(x)
  
  ## store original variance and sd
  varo <- var(x)
  xsd <- sd(x)
  
  ## taper with cosine bell
  trt <- spec.taper(trm, 0.05)
  
  ## rescale so that mean is exactly zero
  trt <- trt - mean(trt)
  
  ## find power of 2 greater than double length of data
  ll <- length(trt)
  pow2 <- 2^c(1:12)
  pow <- pow2[which(pow2 > 2 * ll)[1]]
  
  ## pad data with 0 to length of next power of 2
  padlen <- pow - ll
  trp <- c(trt, rep(0, padlen))
  
  ## scale variance of trp to original variance of u
  sdx <- sd(x)
  sdp <- sd(trp)
  meanp <- mean(trp)
  trtemp <- trp - meanp
  trtemp <- trtemp * (sdx / sdp)
  trp <- trtemp + meanp
  
  ## compute discrete fourier transform on tapered and padded series
  z <- fft(trp)
  
  ## compute periodogram; this differs from MATLAB version according
  ## to a scaling factor of ~ 2
  Pyy  <- Re(z * Conj(z))/padlen
  
  ## sample Gaussian noise and compute mu for n simulation runs
  M <- length(trp)/2
  Z <- matrix(rnorm(nh * length(trp) * 2), ncol = nh)
  k <- t(1:length(trp))
  i1 <- 2 * k - 1
  i2 <- 2 * k
  
  term1 <- matrix(complex(pow, Z[i1,], Z[i2,]), ncol = nh)
  term2 <- sqrt(Pyy/length(trp))
  term2 <- matrix(rep(term2, nh), ncol = nh)
  mu <- term1 * term2
  
  V <- fft(mu)
  
  Vr <- Re(V)
  Vi <- Im(V)
  D1 <- Vr[1:M,]
  D2 <- Vi[1:M,]
  D <- cbind(D1, D2)
  D <- D[1:ll, 1:n]
  
  ## scale to original mean and variance
  Dmean <- matrix(rep(colMeans(D), each = ll), ncol = n)
  Dsd <- matrix(rep(apply(D, 2, sd), each = ll), ncol = n)
  Dz <- (D - Dmean) / Dsd
  out_x <- Dz * matrix(xsd, nrow = ll, ncol = n) + mean(x)
  out_x
}

