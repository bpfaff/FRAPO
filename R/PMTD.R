##
## Minimum Tail dependence Portfolio
##
PMTD <- function(Returns, method = c("EmpTC", "EVT"), k = NULL,
                 percentage = TRUE, optctrl = ctrl(), ...){
  if (is.null(dim(Returns))) {
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()
  V <- tdc(x = Returns, method = method, k = k, ...)
  N <- ncol(Returns)
  ## QP
  ## Budget
  A <- matrix(rep(1, N), nrow = 1)
  b <- 1
  ## Nonnegativity constraint (inequality)
  nno1 <- nnoc(G = -diag(N), h = rep(0, N))
  ## Call to cccp
  opt <- cccp(P = 2 * V, q = rep(0, N), A = A, b = b, cList = list(nno1), optctrl = optctrl)
  w <- drop(getx(opt))
  ## re-scaling weights by assets' sd 
  sd <- sqrt(diag(cov(Returns)))
  w <- w / sd
  w <- w / sum(w)
  names(w) <- colnames(Returns)
  if(percentage) w <- w * 100
  new("PortSol", weights = w, opt = list(opt), type = "Minimum Tail Dependent", call = call)
}
