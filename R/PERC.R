##
## Equal Risk Contribution Portfolio
##
PERC <- function (Sigma, par = NULL, percentage = TRUE, optctrl = ctrl(), ...){
  if(!isSymmetric(Sigma)){
    stop("Matrix provided for Sigma is not symmetric.\n")
  }
  N <- ncol(Sigma)
  mrc <- rep(1/N, N)
  if(is.null(par)){
    par <- mrc
  } else {
    if(length(par) != N){
      stop("Length of 'par' not comformable with dimension of 'Sigma'.\n")
    }
  }
  call <- match.call()
  ## calling rp() from cccp
  opt <- rp(x0 = par, P = Sigma, mrc = mrc, optctrl = optctrl)
  w <- drop(getx(opt))
  w <- w / sum(w)
  if(percentage) w <- w * 100
  if(is.null(dimnames(Sigma))){
    names(w) <- paste("Asset", 1:N, sep = "")
    } else {
      names(w) <- colnames(Sigma)
    }
  new("PortSol", weights = w, opt = list(opt),
      type = "Equal Risk Contribution", call = call)
}
