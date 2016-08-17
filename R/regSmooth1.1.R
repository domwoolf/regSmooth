#' Create divided difference matrix
#' @param x vector to process
#' @param d Dimension
#' @export
#' @importClassesFrom Matrix dgCMatrix
#' @importMethodsFrom Matrix diff
#' @return Output is divided difference matrix for use internally in regSmooth
#' @keywords internal
divided.diff <- function (x, d){
  m  <- length(x)
  if (d == 0) {
    D <- Matrix::sparseMatrix(1:m,  1:m,  x  =  1,  dims  =  c(m,m))}
  else {
    dx  <-  x[(d + 1):m] - x[1:(m - d)]
    dl <- length(dx)
    V <- Matrix::sparseMatrix(1:dl,  1:dl,  x  =  1/dx,  dims  =  c(dl,dl))
    D  <-  d * V %*% diff(divided.diff(x, d - 1))}
  return (D)
}


#' Data smoothing by regularisation.
#'
#' Provides data smoothing by Tikhonov regularisation.
#'
#' @details
#' This function is based on the method of Eilers (2003), and its implementation by Stickel (2010).
#' @param x x values of data to smooth. Numeric vector.If xhat is not supplied, then x must be strictly monotonic increasing.
#' @param y y values of data to smooth. Numeric vector. Must be same length as x.
#' @param lambda smoothing parameter (positive real numeric). Higher values give greater smoothing. Smaller values follow data more closely.
#' @param d Optional order of smoothing derivative. Default=2 (2nd order).
#' @param xhat An optional numeric vector of x-values to use for the smooth curve. If provided, xhat must be strictly monotonically increasing, and must at least span the data. If xhat is not provided, then the smoothed y values will be have x coordinates at x.
#' @param w Optional weighting values for fitting each point in the data. Numeric vector. If not provided, then all data points are given equal weighting. Must have same length as x and y.
#' @param relative Logical (default = False) flag to use relative differences for the goodness of fit term.  Conflicts with the weights (w) optional argument. If relative=True and weights are provided as vector w, then a warning will be issued, and relative set to False.
#' @param midpr Logical (default = False) flag to use the midpoint rule for the integration terms rather than a direct sum. This option conflicts with the option "xhat". If midpr=True and xhat is provided, then a warning will be issued, and midpr set to False.
#' @importClassesFrom Matrix dgCMatrix
#' @importMethodsFrom Matrix t diag solve
#' @export
#' @return Provides a list of two elements:
#'  1) yhat = A numeric vector (yhat) of smoothed y values. Vector is of same length as xhat (if provided), or x if not. If xhat is provided, then the smoothed function is given by the data points (xhat,yhat), otherwise, the smoothed curve is (x,yhat).
#'  2) variance = the variance of the fitted (smoothed model).
#' @seealso \code{\link{regSmoothAuto}}
#' @references
#'   Eilers, P., (2003). A Perfect Smoother. anal. chem. 75, 3631-3636.
#'
#'   Stickel, J. (2010). Data smoothing and numerical differentiation by a regularization method. Comp. chem. eng. 34.4 467-475.
#' @author
#'   Dominic Woolf.
#'   d.woolf@cornell.edu
#' @keywords Smoothing
#' @examples
#'   x <- 1:1000
#'   y <- sin(x/100) + 0.3*cos(x/10)                       # Generate a curve
#'   y.noisy <- y + rnorm(1000, sd=0.1)                    # Add some noise
#'   y.smoothed <- regSmooth(x=x, y=y.noisy, lambda=1e-6)  # Smooth the noisy data
#'   plot  (y.noisy ~ x)                                   # Plot noisy data
#'   lines (y.smoothed$yhat ~ x, col="red", lw=2)     # Draw smoothed curve in red
#'   lines (y ~ x, col="darkgreen", lw=2, lt=1)       # Original function (without noise) in green

regSmooth <- function (x, y, lambda, d=2, xhat, w, relative, midpr) {
  # check for required parameters (x, y, and smoothing parameter "lambda")
  if (missing(x)) stop("Need to specify x values")
  if (missing(y)) stop("Need to specify y values")
  if (missing(lambda)) stop("Need to specify lambda")

  ## Defaults if not provided
  #if (missing(d)) d <- 3
  if (missing(midpr)) midpr <- F
  xhat.provided  <-  T
  if (missing(xhat)) {
    xhat <- x
    xhat.provided  <-  F}
  if (xhat.provided & midpr) {
      warning("midpointrule is not implemented if xhat is provided. midpr set to false.")
      midpr <- F  }
  weights.provided <- T
  if (missing(w)) weights.provided <- F
  if (missing(relative)) relative <- F
  if (relative & weights.provided) {
    warning("relative differences are not used if a weighting vector is provided")
    relative <- F}
  N <- length(x)
  Nhat <- length(xhat)

  if (!(all(diff(xhat)>0))) { #test that x is strict monotonic increasing
    if (xhat.provided) {
      stop ("xhat must be monotonically increasing")
      } else {
        stop ("x must be monotonically increasing if xhat is not provided") }}

  # test that xhat spans x
  if ( (min(x) < min(xhat)) | (max(xhat) < max(x)) )  stop("xhat must at least span the data")

  ########################################
  # start function proper
  ########################################
  M <- Matrix::sparseMatrix(1:N,  1:N,  x  =  1,  dims  =  c(N,N))
  idx <- unlist(lapply(xhat, function(xi) which.min(abs(xi - x)))) #  for unequally spaced xhat
  M <- M[,idx]
  D <- divided.diff(xhat,d)

  ## construct weighting matrices W and U
  if (weights.provided) { # use weighting provided
    W <- diag(w)
  } else {
    if (relative) {
      yl <- length(y)
      Yinv <- Matrix::sparseMatrix(1:yl,  1:yl,  x  =  1/y,  dims  =  c(yl,yl))
      W <- Yinv^2
    } else {
      W <- Matrix::sparseMatrix(1:N,  1:N,  x  =  1,  dims  =  c(N,N))
      }
    }

  if (midpr) { ## use midpoint rule integration (rather than simple sums)
    Bhat <- Matrix::Matrix(0,nrow=N,ncol=N)
    Bhat[1:N-1,2:N] <- Matrix::Diagonal(N-1)
    Bhat[2:N,1:N-1] <- Bhat[2:N,1:N-1]-Matrix::Diagonal(N-1)
    Bhat[1,1] <- -1
    Bhat[N,N] <- 1
    B <- Matrix::sparseMatrix(1:N, 1:N, x = (Bhat %*% x)/2, dims=c(N,N))
    if ( floor(d/2) == d/2 ) { # test if d is even
      dh <- d/2
      Btilda <- B[(dh+1):(N-dh) , (dh+1):(N-dh)]
      } else { # d is odd
        dh <- ceiling(d/2)
        Btilda <- B[dh:(N-dh) , dh:(N-dh)]  }
    W <- W%*%B
    U <- Btilda
  } else {
    # W = W*speye(Nhat)
    U <- Matrix::sparseMatrix(1:(Nhat-d),  1:(Nhat-d),  x  =  1,  dims  =  c((Nhat-d),(Nhat-d)))
    }

  ## Do the smoothing
  delta <- sum(Matrix::diag(Matrix::t(D)%*%D)) / Nhat^(2+d)
  yhatA <- (Matrix::t(M)%*%W%*%M + lambda*delta^(-1)*Matrix::t(D) %*% U %*% D)
  yhatB <- Matrix::t(M)%*%W%*%y
  yhat <- solve(yhatA) %*% yhatB
  H1 = Matrix::t(M)%*%W%*%M + (lambda/delta)*Matrix::t(D)%*%U%*%D
  H = M %*% ((solve(H1) %*% Matrix::t(M)) %*% W)
  v = Matrix::t(M%*%yhat - y) %*% (M%*%yhat - y) / N / (1 - sum(Matrix::diag(H))/N)^2
  return(list(yhat=yhat[,1], variance=v[1,1]))
}


#' Data smoothing by regularisation.
#'
#' Function to smooth data using Tikhonov regularisation,
#' with optional optimisation of the regularisation parameter by generalised cross validation.
#' @details
#' This function is based on the method of Eilers (2003), and its implementation by Stickel (2010).
#' @param x x values of data to smooth. Numeric vector.If xhat is not supplied, then x must be strictly monotonic increasing.
#' @param y y values of data to smooth. Numeric vector. Must be same length as x.
#' @param lambda smoothing parameter (positive real numeric). Higher values give greater smoothing. Smaller values follow data more closely.
#' @param d Optional order of smoothing derivative. Default=2 (2nd order derivatives).
#' @param relative Logical (default = False) flag to use relative differences for the goodness of fit term.  Conflicts with the weights (w) optional argument. If relative=True and weights are provided as vector w, then a warning will be issued, and relative set to False.
#' @param midpr Logical (default = False) flag to use the midpoint rule for the integration terms rather than a direct sum. This option conflicts with the option "xhat". If midpr=True and xhat is provided, then a warning will be issued, and midpr set to False.
#' @param maxit Maximum number of iterations to use when solving for lambda.
#' @param guess.lambda Starting value for lambda optimisation.
#' @param ... Additional optional parameters to pass to regSmooth:
#' w =  Optional weighting values for fitting each point in the data (numeric vector). If not provided, then all data points are given equal weighting. Must have same length as x and y.
#' xhat = An optional numeric vector of x-values to use for the smooth curve. If provided, xhat must be strictly monotonically increasing, and must at least span the data. If xhat is not provided, then the smoothed y values will be have x coordinates at x.
#' @export
#' @return Provides a list of two elements:
#'  1) yhat = A numeric vector (yhat) of smoothed y values. Vector is of same length as xhat (if provided), or x if not. If xhat is provided, then the smoothed function is given by the data points (xhat,yhat), otherwise, the smoothed curve is (x,yhat).
#'  2) variance = the variance of the fitted (smoothed model).
#' @seealso \code{\link{regSmooth}}
#' @references
#'   Eilers, P. (2003). A Perfect Smoother. Anal. chem. 75, 3631-3636.
#'
#'   Stickel, J. (2010). Data smoothing and numerical differentiation by a regularization method. Comp. chem. eng. 34.4 467-475.
#' @author
#'   Dominic Woolf.
#'   d.woolf@cornell.edu
#' @keywords Smoothing
#' @examples
#'   x <- 1:1000
#'   y <- sin(x/100) + 0.3*cos(x/10)              # Generate a curve
#'   y.noisy <- y + rnorm(1000, sd=0.1)           # Add some noise
#'   y.smoothed <- regSmoothAuto(x=x, y=y.noisy)  # Smooth the noisy data
#'   plot  (y.noisy ~ x)                          # Plot the noisy data
#'   lines (y.smoothed$yhat ~ x, col="red", lw=2) # Draw the smoothed curve in red
#'   lines (y ~ x, col="darkgreen", lw=2, lt=1)   # Original function (without noise) in green

regSmoothAuto <- function (x, y, lambda, d=2, relative=F, midpr=F, maxit=1000, guess.lambda=1e-3, ...){
  payoff <- function(try.lambda){
    return(regSmooth(x=x, y=y, lambda=try.lambda, relative=relative, midpr=midpr)$variance)
  }
  if (missing(lambda)) {# solve for optimal lambda if user-defined value not provided
    lambda <- nlm(f=payoff, p=guess.lambda, iterlim=maxit)$estimate
  }
  smooth.result  <-  regSmooth(x=x, y=y, d=d, lambda=lambda, relative=relative, midpr=midpr,...)
  return (list(yhat=smooth.result$yhat, lambda=lambda, variance=smooth.result$variance))
}
