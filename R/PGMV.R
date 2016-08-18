##
## Global Minimum Variance Portfolio
##
PGMV <- function(Returns, percentage = TRUE, optctrl = ctrl(), ...){
  if(is.null(dim(Returns))){
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()  
  V <- cov(Returns, ...)
  N <- ncol(Returns)
  ## Budget
  A <- matrix(rep(1, N), nrow = 1)
  b <- matrix(1)
  ## Nonnegativity constraint (inequality)
  nno1 <- nnoc(G = -diag(N), h = matrix(rep(0, N)))
  ## Call to cccp
  opt <- cccp(P = V, q = rep(0, N), A = A, b = b, cList = list(nno1), optctrl = optctrl)
  ## Recovering weights for assets
  w <- drop(getx(opt))
  names(w) <- colnames(Returns)
  if(percentage) w <- w * 100
  new("PortSol", weights = w, opt = list(opt), type = "Global Minimum Variance", call = call)
}
