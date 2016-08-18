##
## Most Diversified Portfolio
##
PMD <- function(Returns, percentage = TRUE, optctrl = ctrl(),...){
  if(is.null(dim(Returns))){
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()  
  V <- cov(Returns, ...)
  C <- cov2cor(V)
  N <- ncol(Returns)
  ## QP
  ## Budget
  A <- matrix(rep(1, N), nrow = 1)
  b <- 1
  ## Nonnegativity constraint (inequality)
  nno1 <- nnoc(G = -diag(N), h = rep(0, N))
  ## Call to cccp
  opt <- cccp(P = C, q = rep(0, N), A = A, b = b, cList = list(nno1), optctrl = optctrl)
  ## Recovering weights for assets
  w <- drop(getx(opt)) / sqrt(diag(V))  
  names(w) <- colnames(Returns)
  wnorm <- w / sum(w)
  if(percentage) wnorm <- wnorm * 100
  new("PortSol", weights = wnorm, opt = list(opt), type = "Most Diversifified", call = call)
}
