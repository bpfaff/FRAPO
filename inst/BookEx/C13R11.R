## Defining function for GARCH-forecast
CondVolaFcst <- function(x){
    m <- garchFit(formula = ~ garch(1,1),
                  data = x, trace = FALSE)
    mp <- as.numeric(predict(m, n.ahead = 1))[3]
    mp
}
## Defining functions for entropy pooling
f0 <- function(v, p, Aeq, beq){
    x <- exp(log(p) - 1 - crossprod(Aeq, v))
        x = apply(cbind(x, 0), 1, max)
        L = t(x) %*% (log(x) - log(p) + crossprod(Aeq, v)) -
            crossprod(beq, v)
       -L
}
gf <- function(v, p, Aeq, beq){
    x <- exp(log(p) - 1 - crossprod(Aeq, v))
    beq - Aeq %*% x 
}
ep <- function(x0, Aeq, beq, pprior){
    vopt <- try(optim(par = x0, fn = f0, gr = gf,
              Aeq = Aeq, beq = beq, p = pprior, method = "BFGS"))
    if(class(vopt) == "try-error"){
        return(pprior)
    } else {
        v <- vopt$par
        pbar <- exp(log(pprior) - 1 - crossprod(Aeq, v))
        return(pbar / sum(pbar))
    }
}
