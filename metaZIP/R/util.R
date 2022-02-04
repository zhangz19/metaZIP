

#++++++++++++ Zhen Zhang (zhangquake1@outlook.com), 2022-Jan
require(DEoptim)
require(gaussquad)
nquad <- 21
rule <- ghermite.h.quadrature.rules(nquad, mu=0, normalized=F)[[nquad]]
require(statmod)
gauss.15 <- gauss.quad(15,"hermite")

# eps <- .Machine$double.eps
# eps <- 0.001

metaZIPs <- function(y, x, ni, method=c('POIF','POIR','ZIPF','ZIPR'), para0=NULL, conf.level=.95, eps=.Machine$double.eps){ 
  #meta-ZIP model with single-arm trials on 2 groups (x=0,1): e.g. EUS-FNA in the manuscript
  #each study has y events out of ni.
  n <- length(ni)  #number of studies
  
  # negative log-likelihood functions
  nll_POIF <- function(para){  #para=c(alph, bet, tau)
    alph <- para[['alph']];  bet <- para[['bet']];  tau <- para[['tau']]
    # alph <- para[1];  bet <- para[2];  tau <- para[3]
    ntx <- ni*exp(tau*x)
    val <- n*(alph*log(bet)-lgamma(alph)) + sum( lgamma(y+alph) - (y+alph)*log(bet+ntx) + y*log(ntx) - lgamma(y+1))  #log-likelihood
    return(-val)
  }
  
  nll_POIR <- function(para){  #para=c(alph, bet, mu, sig)
    # alph <- para[['alph']];  bet <- para[['bet']];  mu <- para[['mu']];  sig <- para[['sig']]
    # alph <- para[1];  bet <- para[2];  mu <- para[3];  sig <- para[4]
    alph <- exp(para[['alph']]);  bet <- exp(para[['bet']]);  mu <- para[['mu']];  sig <- exp(para[['sig']])
    lik <- rep(NA, n)
    for(i in 1:n){
      if(x[i]==1){
        # fn <- function(t0){
        #   ntx <- ni[i]*exp(t0*x[i])
        #   ll <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) - (y[i]+alph)*log(bet+ntx) +
        #     y[i]*log(ntx) - lgamma(y[i]+1)
        #   ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
        #   return(exp(ll))
        # }
        # lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        # # lik[i] <-  integrate(fn, lower=-Inf, upper=Inf)$value
        # # pt <- seq(-10,  10, len=1e4); sum(fn(pt)*diff(pt)[1])
        
        #++++ alternative
        # k <- gamma(y[i]+alph)/factorial(y[i])
        # lik[i] <- sum( ((bet+ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2))^(-y[i]-alph))*gauss.15$weights )* exp(0.5*( (mu/sig+y[i]*sig)^2-(mu/sig)^2 )) *(ni[i]^y[i])*(bet^alph)*k/(gamma(alph)*sqrt(pi))
        
        # sum( ((bet+ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2))^(-y[i]-alph))*gauss.15$weights )
        # is equivalent to the following: 
        # fn <- function(t0){
        #   ll <- exp(-t0^2 - (y[i]+alph)*log(bet+ni[i]*exp(sqrt(2)*sig*t0+mu+y[i]*sig^2)) )
        #   return((ll))
        # }
        # ghermite.h.quadrature(fn, rule, weighted=FALSE)
        
        #++++ alternative: this seems to be more numerically stable comparing ghermite.h.quadrature() when using optim
        k <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) + y[i]*log(ni[i]) - lgamma(y[i]+1) - 0.5*log(pi) + 0.5*( (mu/sig+y[i]*sig)^2-(mu/sig)^2 )
        lik[i] <- sum(exp( k + log(gauss.15$weights) - (y[i]+alph)*log(bet+ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2)) ))
        
        if(lik[i]==0){  #this may happen numerically when y[i] is large
          fn <- function(t0){
            ntx <- ni[i]*exp(t0*x[i])
            ll <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) - (y[i]+alph)*log(bet+ntx) +
              y[i]*log(ntx) - lgamma(y[i]+1)
            ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
            return(exp(ll))
          }
          lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        }
        
      }
      if(x[i]==0){
        ntx <- ni[i]
        lik[i] <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) - (y[i]+alph)*log(bet+ntx) +
          y[i]*log(ntx) - lgamma(y[i]+1)
        lik[i] <- exp(lik[i])
      }
    }
    # inds <- which(lik==0); if(length(inds)) lik <- lik[-inds]#lik[inds] <- eps
    return(-sum(log(lik)))
  }
  
  nll_ZIPF <- function(para){ 
    alph <- para[['alph']];  bet <- para[['bet']];  tau <- para[['tau']];  c <- para[['c']];  d <- para[['d']]
    ntx <- ni*exp(tau*x)
    i0 <- which(y==0)
    fac <- (1+ntx[i0]/bet)^(-alph)
    n1 <- length(  i1 <- which(y!=0)  )
    val <- sum(log( d/(c+d)*(1 - fac) + fac )) + 
      n1*log(c/(c+d)) + n1*alph*log(bet) - n1*lgamma(alph) + sum( 
        lgamma(y[i1]+alph) - (y[i1]+alph)*log(bet + ntx[i1]) + y[i1]*log(ntx[i1]) - lgamma(y[i1]+1)  )
    return(-val)
  }
  
  nll_ZIPR <- function(para){ 
    alph <- para[['alph']];  bet <- para[['bet']];  mu <- para[['mu']];  sig <- para[['sig']];  c <- para[['c']];  d <- para[['d']]
    ll <- 0
    n0 <- length(i0 <- which(y==0 & x==1))
    if(n0>0){
      h1i <- rep(NA, n0)
      for(i in 1:n0){
        # integrate(function(t) exp(-t^2)/(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu))^alph, lower=-5, upper=5)$value
        h1i[i] <- ghermite.h.quadrature(function(t) exp(-t^2- alph*log(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu))), rule, weighted=FALSE)
        # pt <- seq(-5, 5, len=1e5)
        # h1i[i] <- sum( exp(-pt^2)/(bet + ni[i0[i]]*exp(sqrt(2)*sig*pt+mu))^alph * diff(pt)[1] )
      }
      ll <- ll + sum(log(  d/(c+d)*(1-bet^alph*h1i/sqrt(pi)) + bet^alph*h1i/sqrt(pi)  ))
    }
    n0 <- length(i0 <- which(y!=0 & x==1))
    if(n0>0){
      h2i <- rep(NA, n0)
      for(i in 1:n0){
        # h2i[i] <- ghermite.h.quadrature(function(t) exp(-t^2)/( bet + ni[i0[i]]*
        # exp(sqrt(2)*sig*t+mu+y[i0[i]]*sig^2) )^(y[i0[i]]+alph), rule, weighted=FALSE)
        h2i[i] <- ghermite.h.quadrature(function(t) exp( lgamma(y[i0[i]]+alph) - lgamma(y[i0[i]]+1) + y[i0[i]]*log(ni[i0[i]]) - 0.5*(mu^2/sig^2 - (mu/sig + y[i0[i]]*sig)^2) -t^2- (y[i0[i]]+alph)*log(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu+y[i0[i]]*sig^2))  ), rule, weighted=FALSE)
        # aa <- function(t){
        #   exp(-t^2 - (y[i0[i]]+alph)*log( bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu+y[i0[i]]*sig^2) ))
        # }
        # ghermite.h.quadrature(aa, rule, weighted=FALSE)
        
        # h2i[i] <- integrate(function(t) exp(-t^2)/(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu+y[i0[i]]*sig^2))^(y[i0[i]]+alph), lower=-Inf, upper=Inf)$value
      }
      # ll <- ll + n0*log(c/(c+d)*bet^alph/sqrt(pi)/gamma(alph)) + sum(  lgamma(y[i0]+alph) - lgamma(y[i0]+1) + y[i0]*log(ni[i0]) - 0.5*(mu^2/sig^2 - (mu/sig + y[i0]*sig)^2) + log(h2i)  )
      ll <- ll + n0*log(c/(c+d)*bet^alph/sqrt(pi)/gamma(alph)) + sum(log(h2i))
    }
    n0 <- length(i0 <- which(y==0 & x==0))
    if(n0>0) 
      ll <- ll + sum(log( d/(c+d)*(1-bet^alph/(bet+ni[i0])^alph) + bet^alph/(bet+ni[i0])^alph ))
    n0 <- length(i0 <- which(y!=0 & x==0))
    if(n0>0) 
      ll <- ll + n0*log(c/(c+d)*bet^alph/gamma(alph)) + 
      sum( y[i0]*log(ni[i0]) + lgamma(y[i0]+alph) - lgamma(y[i0]+1) - (y[i0]+alph)*log(bet+ni[i0]) )
    return(-ll)
  }
  
  summaryFit <- function(fit, method){
    if(method %in% c('POIF','ZIPF')){
      logRisk <- fit$par[['tau']]
      se <- sqrt(diag(solve(fit$hessian)))[['tau']]
    }
    if(method %in% c('POIR','ZIPR')){
      logRisk <- fit$par[['mu']]
      se <- fit$par[['sig']] #sqrt(diag(solve(-fit$hessian)))[['mu']]
    }
    mgn <- qnorm(1-(1-conf.level)/2)*se
    bds <- exp(logRisk+c(-1,1)*mgn)
    out <- c('Log Likelihood'=-fit$value, 'Estimate'=logRisk, 'SE'=se, '95% CI lower'=bds[1], 'upper'=bds[2], 'npara'=length(fit$par), 'AIC'=2*length(fit$par)+2*fit$value, 'BIC'=log(n)*length(fit$par)+2*fit$value)
    paraEst <- fit$par
    return(list(out=out, paraEst=paraEst))
  }
  
  nmed <- length(method)
  out <- matrix(NA, nmed, 8);  paraEst <- matrix(NA, nmed, 6)
  
  # fit a base model: 'POIF'
  para <- c(alph=1, bet=1, tau=0) 
  if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']])
  m_POIF <- optim(para, nll_POIF, method='L-BFGS-B', lower=c(eps,eps,-Inf), upper=c(Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #minimize negative-log-likelihood

  if('POIF' %in% method){
    s_POIF <- summaryFit(m_POIF, 'POIF')
    k <- which(method=='POIF');  out[k,] <- s_POIF$out
    paraEst[k, 1:length(s_POIF$paraEst)] <- s_POIF$paraEst
  }
  
  if('POIR' %in% method || 'ZIPR' %in% method){ #use initial value from 'POIF'
    para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], mu=m_POIF$par[['tau']], sig=1)
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']])
    para[c(1,2,4)] <- log(para[c(1,2,4)])
    m_POIR <- optim(para, nll_POIR, hessian=TRUE, control=list(trace=0, maxit=500))
    m_POIR$par[c(1,2,4)] <- exp(m_POIR$par[c(1,2,4)])
    s_POIR <- summaryFit(m_POIR, 'POIR')
    k <- which(method=='POIR')
    if(length(k)){ out[k,] <- s_POIR$out;  paraEst[k, 1:length(s_POIR$paraEst)] <- s_POIR$paraEst } 
  }
  
  if('ZIPF' %in% method){ #use initial value from 'POIF'
    para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], tau=m_POIF$par[['tau']], c=1, d=1)
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']], c=para0[['c']], d=para0[['d']])
    m_ZIPF <- optim(para, nll_ZIPF, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1
    s_ZIPF <- summaryFit(m_ZIPF, 'ZIPF')
    k <- which(method=='ZIPF');  out[k,] <- s_ZIPF$out;  paraEst[k,c(1,2,3,5,6)] <- s_ZIPF$paraEst
  }
  
  if('ZIPR' %in% method){ #use initial value from 'POIR'
    para <- c(alph=m_POIR$par[['alph']], bet=m_POIR$par[['bet']], mu=m_POIR$par[['mu']], sig=m_POIR$par[['sig']], c=1, d=1) 
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']], c=para0[['c']], d=para0[['d']])
    m_ZIPR <- optim(para, nll_ZIPR, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1, 
    s_ZIPR <- summaryFit(m_ZIPR, 'ZIPR')
    k <- which(method=='ZIPR');  out[k,] <- s_ZIPR$out;  paraEst[k,] <- s_ZIPR$paraEst
  }
  
  out <- as.data.frame(out);  row.names(out) <- method
  names(out) <- c('Log Likelihood','Estimate','SE','95% CI lower','upper','npara','AIC','BIC')
  paraEst <- as.data.frame(paraEst);  row.names(paraEst) <- method
  names(paraEst) <- c('alpha','beta','tau(mu)','sigma','c','d')
  
  
  return(list(out=out, paraEst=paraEst))
}


simZIPs <- function(para0, zeroInflated=TRUE, randomTau=TRUE, ni0=NULL, x0=NULL, n=NULL, seeds=123, ratio_t_c=c(0.5, 0.5)){
  #simulate data using parameter para0
  simN <- FALSE;  if(is.null(n)) n <- length(ni0) else simN <- TRUE  #default sample size using originally observed data
  set.seed(seeds)
  alph <- para0[['alpha']];  bet <- para0[['beta']]
  if(randomTau){ mu <- para0[['tau(mu)']];  sig <- para0[['sigma']] } else tau <- para0[['tau(mu)']]
  if(zeroInflated){ c <- para0[['c']];  d <- para0[['d']] } 
  ys <- rep(NA, n)
  xs <- x0;  if(is.null(x0)) xs <- sample(c(0,1), n, replace=TRUE, prob=ratio_t_c)
  nis <- ni0;  if(simN)  nis <- round(runif(n, 20, 100))  #simulate study size
  # if(simN)  nis <- round(sample(ni0, n, replace=TRUE))  #resampling from the observed study size
  pis <- rep(1,n);  if(zeroInflated) pis <- rbeta(n, c, d)
  ai <- rep(NA, n)
  for(i in 1:n){
    xi <- rgamma(1, shape=alph, rate=bet)
    if(randomTau) tauUse <- rnorm(1, mu, sig)  else tauUse <- tau
    lambda <- xi*exp(tauUse*xs[i])
    ai[i] <- rbinom(1,1,pis[i]); if(ai[i]==1) ys[i] <- rpois(1, nis[i]*lambda[1]) else ys[i] <- 0
  }
  return(data.frame(y=ys, x=xs, ni=nis))
}



metaZIPs0 <- function(y, x, ni, method=c('POIF','POIR','ZIPF','ZIPR'), para0=NULL, conf.level=.95, eps=.Machine$double.eps){ 
  #meta-ZIP model with single-arm trials on 2 groups (x=0,1): e.g. EUS-FNA in the manuscript
  #each study has y events out of ni.
  n <- length(ni)  #number of studies
  
  # negative log-likelihood functions
  nll_POIF <- function(para){  #para=c(alph, bet, tau)
    xi <- para[['xi']];  tau <- para[['tau']]
    lambda <- xi*exp(tau*x)
    ll <- sum(dpois(y, ni*lambda, log=TRUE))
    return(-ll)
  }
  
  nll_POIR <- function(para){  #para=c(alph, bet, mu, sig)
    # alph <- para[['alph']];  bet <- para[['bet']];  mu <- para[['mu']];  sig <- para[['sig']]
    # alph <- para[1];  bet <- para[2];  mu <- para[3];  sig <- para[4]
    xi <- exp(para[['xi']]);  mu <- para[['mu']];  sig <- exp(para[['sig']])
    lik <- rep(NA, n)
    for(i in 1:n){
      if(x[i]==1){
        #++++ alternative: this seems to be more numerically stable comparing ghermite.h.quadrature() when using optim
        k <- y[i]*log(ni[i]*xi) - lgamma(y[i]+1) - 0.5*log(pi) + 0.5*( (mu/sig+y[i]*sig)^2-(mu/sig)^2 )
        lik[i] <- sum(exp( k + log(gauss.15$weights) - xi*ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2) ))
      }
      if(x[i]==0) lik[i] <- dpois(y[i], ni[i]*xi)
    }
    return(-sum(log(lik)))
  }
  
  nll_ZIPF <- function(para){ 
    # xi <- para[['xi']];  tau <- para[['tau']];  c <- para[['c']];  d <- para[['d']]
    xi <- exp(para[['xi']]);  tau <- para[['tau']];  c <- exp(para[['c']]);  d <- exp(para[['d']])
    i0 <- which(y==0)
    n1 <- length(i1 <- which(y!=0))
    ll <- n1*log(c/(c+d)) + 
      sum( dpois(y[i1], ni[i1]*xi*exp(tau*x[i1]), log=TRUE) ) + 
    sum(log( d/(c+d) + c/(c+d)*dpois(y[i0], ni[i0]*xi*exp(tau*x[i0])) ))
    return(-ll)
  }
  
  nll_ZIPR <- function(para){ 
    # xi <- para[['xi']];  mu <- para[['mu']];  sig <- para[['sig']];  c <- para[['c']];  d <- para[['d']]
    xi <- exp(para[['xi']]);  mu <- para[['mu']];  sig <- exp(para[['sig']]);  c <- exp(para[['c']]);  d <- exp(para[['d']])
    i0 <- which(y==0)
    n1 <- length(i1 <- which(y!=0))
    ll <- 0
    for(i in 1:n){
      # h <- ghermite.h.quadrature(function(tau) dpois(y[i], ni[i]*xi*exp(tau*x[i]))*dnorm(tau, mu, sig), rule, weighted=FALSE)
      if(y[i]!=0){
        if(x[i]==1){
          k <- y[i]*log(ni[i]*xi) - lgamma(y[i]+1) - 0.5*log(pi) + 0.5*( (mu/sig+y[i]*sig)^2-(mu/sig)^2 )
          h <- sum(exp( k + log(gauss.15$weights) - xi*ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2) ))
          if(h==0){
            fn <- function(t0){
              lambda <- xi*exp(t0)
              ll <- -ni[i]*lambda + yi[i]*log(ni[i]*lambda) - lgamma(yi[i]+1) + dnorm(t0, mu, sig, log=TRUE)
              return(exp(ll))
            }
            h <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
          }
        }
        if(x[i]==0) h <- dpois(y[i], ni[i]*xi)
        ll <- ll + log(c/(c+d)) + log(h)
      }
      if(y[i]==0){
        if(x[i]==1){
          k <- - 0.5*log(pi)
          h <- sum(exp( k + log(gauss.15$weights) - xi*ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu) ))
          if(h==0){
            fn <- function(t0){
              lambda <- xi*exp(t0)
              ll <- -ni[i]*lambda + dnorm(t0, mu, sig, log=TRUE)
              return(exp(ll))
            }
            h <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
          }
        }
        if(x[i]==0) h <- dpois(0, ni[i]*xi)
        ll <- ll + log(d/(c+d) + c/(c+d)*h)
      }
    }
    return(-ll)
  }
  
  summaryFit <- function(fit, method){
    if(method %in% c('POIF','ZIPF')){
      logRisk <- fit$par[['tau']]
      se <- sqrt(diag(solve(fit$hessian)))[['tau']]
    }
    if(method %in% c('POIR','ZIPR')){
      logRisk <- fit$par[['mu']]
      se <- fit$par[['sig']] #sqrt(diag(solve(-fit$hessian)))[['mu']]
    }
    mgn <- qnorm(1-(1-conf.level)/2)*se
    bds <- exp(logRisk+c(-1,1)*mgn)
    out <- c('Log Likelihood'=-fit$value, 'Estimate'=logRisk, 'SE'=se, '95% CI lower'=bds[1], 'upper'=bds[2], 'npara'=length(fit$par), 'AIC'=2*length(fit$par)+2*fit$value, 'BIC'=log(n)*length(fit$par)+2*fit$value)
    paraEst <- fit$par
    return(list(out=out, paraEst=paraEst))
  }
  
  nmed <- length(method)
  out <- matrix(NA, nmed, 8);  paraEst <- matrix(NA, nmed, 5)
  
  # fit a base model: 'POIF'
  para <- c(xi=1, tau=0) 
  # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']])
  m_POIF <- optim(para, nll_POIF, method='L-BFGS-B', lower=c(eps,-Inf), upper=c(Inf,Inf), hessian=TRUE, control=list(trace=0)) #minimize negative-log-likelihood
  
  if('POIF' %in% method){
    s_POIF <- summaryFit(m_POIF, 'POIF')
    k <- which(method=='POIF');  out[k,] <- s_POIF$out
    paraEst[k, 1:length(s_POIF$paraEst)] <- s_POIF$paraEst
  }
  
  if('POIR' %in% method || 'ZIPR' %in% method){ #use initial value from 'POIF'
    para <- c(xi=m_POIF$par[['xi']], mu=m_POIF$par[['tau']], sig=1)
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']])
    para[c(1,3)] <- log(para[c(1,3)])
    m_POIR <- optim(para, nll_POIR, hessian=TRUE, control=list(trace=0, maxit=500))
    m_POIR$par[c(1,3)] <- exp(m_POIR$par[c(1,3)])
    s_POIR <- summaryFit(m_POIR, 'POIR')
    k <- which(method=='POIR')
    if(length(k)){ out[k,] <- s_POIR$out;  paraEst[k, 1:length(s_POIR$paraEst)] <- s_POIR$paraEst } 
  }
  
  if('ZIPF' %in% method){ #use initial value from 'POIF'
    # para <- c(xi=m_POIF$par[['xi']], tau=m_POIF$par[['tau']], c=1, d=1)
    # # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']], c=para0[['c']], d=para0[['d']])
    # m_ZIPF <- optim(para, nll_ZIPF, method='L-BFGS-B', lower=c(eps,-Inf,eps,eps), upper=c(Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1
    
    para <- c(xi=log(m_POIF$par[['xi']]), tau=m_POIF$par[['tau']], c=log(1), d=log(1))
    m_ZIPF <- optim(para, nll_ZIPF, hessian=TRUE)
    m_ZIPF$par[c(1,3:4)] <- exp(m_ZIPF$par[c(1,3:4)])
    
    s_ZIPF <- summaryFit(m_ZIPF, 'ZIPF')
    k <- which(method=='ZIPF');  out[k,] <- s_ZIPF$out;  paraEst[k,c(1,2,4,5)] <- s_ZIPF$paraEst
  }
  
  if('ZIPR' %in% method){ #use initial value from 'POIR'
    # para <- c(xi=m_POIR$par[['xi']], mu=m_POIR$par[['mu']], sig=m_POIR$par[['sig']], c=1, d=1) 
    # para <- c(xi=m_ZIPF$par[['xi']], mu=m_ZIPF$par[['tau']], sig=0.02, c=m_ZIPF$par[['c']], d=m_ZIPF$par[['d']])
    # # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']], c=para0[['c']], d=para0[['d']])
    # m_ZIPR <- optim(para, nll_ZIPR, method='L-BFGS-B', lower=c(eps,-Inf,eps,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=5)) #maximize log-likelihood fnscale=-1, 
    
    para <- c(xi=log(m_ZIPF$par[['xi']]), mu=m_ZIPF$par[['tau']], sig=log(0.02), c=log(m_ZIPF$par[['c']]), d=log(m_ZIPF$par[['d']]))
    m_ZIPR <- optim(para, nll_ZIPR, hessian=TRUE)
    m_ZIPR$par[c(1,3:5)] <- exp(m_ZIPR$par[c(1,3:5)])
    
    s_ZIPR <- summaryFit(m_ZIPR, 'ZIPR')
    k <- which(method=='ZIPR');  out[k,] <- s_ZIPR$out;  paraEst[k,] <- s_ZIPR$paraEst
  }
  
  out <- as.data.frame(out);  row.names(out) <- method
  names(out) <- c('Log Likelihood','Estimate','SE','95% CI lower','upper','npara','AIC','BIC')
  paraEst <- as.data.frame(paraEst);  row.names(paraEst) <- method
  names(paraEst) <- c('xi','tau(mu)','sigma','c','d')
  
  
  return(list(out=out, paraEst=paraEst))
}



metaZIPr <- function(y1, n1, y0, n0, method=c('POIF','POIR','ZIPF','ZIPR'), para0=NULL, conf.level=.95, eps=.Machine$double.eps){ 
  #meta-ZIP model with randomized trials on 2 groups (x=0,1): e.g. ROS in the manuscript
  #each study has y1 events out of n1 (treatment), and y0 events out of n0 (control). 
  n <- length(y1)  #number of studies
  
  # negative log-likelihood functions
  nll_POIF <- function(para){  #para=c(alph, bet, tau)
    alph <- para[['alph']];  bet <- para[['bet']];  tau <- para[['tau']]
    ntx <- n1*exp(tau)
    val <- n*(alph*log(bet)-lgamma(alph)) + sum( lgamma(y1+y0+alph) - (y1+y0+alph)*log(bet+ntx+n0) + y1*log(ntx) - lgamma(y1+1) + y0*log(n0) - lgamma(y0+1))  #log-likelihood
    return(-val)
  }
  
  nll_POIR <- function(para){  #para=c(alph, bet, mu, lsig)
    alph <- para[['alph']];  bet <- para[['bet']];  mu <- para[['mu']];  sig <- para[['sig']]
    lik <- rep(NA, n)
    # const <- lgamma(y1+y0+alph) - lgamma(y1+1) - lgamma(y0+1) + y1*log(n1) + y0*log(n0) + alph*log(bet) - lgamma(alph) - 0.5*log(pi) + 0.5*((mu/sig+y1*sig)^2-(mu/sig)^2)
    for(i in 1:n){
      fn <- function(t0){
        ntx <- n1[i]*exp(t0)
        ll <- alph*log(bet) - lgamma(alph) + lgamma(y1[i]+y0[i]+alph) - (y1[i]+y0[i]+alph)*log(bet+ntx+n0[i]) + y1[i]*log(ntx) - lgamma(y1[i]+1) + y0[i]*log(n0[i]) - lgamma(y0[i]+1)
        ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
        return(exp(ll))
      }
      lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
      if(is.nan(lik[i]) || lik[i]==0){  #this may happen numerically in extreme case
        const <- lgamma(y1[i]+y0[i]+alph) - lgamma(y1[i]+1) - lgamma(y0[i]+1) + y1[i]*log(n1[i]) + y0[i]*log(n0[i]) + alph*log(bet) - lgamma(alph) - 0.5*log(pi) + 0.5*((mu/sig+y1[i]*sig)^2-(mu/sig)^2)
        lik[i] <- sum(exp( const + log(gauss.15$weights) - (y0[i]+y1[i]+alph)*log(bet+n1[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y1[i]*sig^2)+n0[i]) ))
      }
    }
    return(-sum(log(lik)))
  }
  
  nll_ZIPF <- function(para){ 
    alph <- para[['alph']];  bet <- para[['bet']];  tau <- para[['tau']];  c <- para[['c']];  d <- para[['d']]
    ll <- 0
    m <- length(i <- which(y1>0 & y0>0))
    if(m>0){
      ntx <- n1[i]*exp(tau)
      ll <- ll + 2*m*log(c/(c+d)) + m*alph*log(bet) - m*lgamma(alph) + sum(
        lgamma(y1[i]+y0[i]+alph) - (y1[i]+y0[i]+alph)*log(bet+ntx+n0[i]) + y1[i]*log(ntx) - lgamma(y1[i]+1) + y0[i]*log(n0[i]) - lgamma(y0[i]+1)
      )
    }
    m <- length(i <- which(y1==0 & y0==0))
    if(m>0){
      ntx <- n1[i]*exp(tau)
      ll <- ll - 2*m*log(c+d) + sum(log(
        d^2 + (bet/(bet+n0[i]))^alph*c*d + (bet/(bet+ntx))^alph*c*d + 
          (bet/(bet+ntx+n0[i]))^alph*c^2
      ))
    }
    m <- length(i <- which(y1==0 & y0>0))
    if(m>0){
      ntx <- n1[i]*exp(tau)
      # ll <- ll + m*log(c/(c+d)^2) + m*alph*log(bet) - m*lgamma(alph) + sum(
      #   lgamma(y0[i]+alph) + y0[i]*log(n0[i]) - lgamma(y0[i]+1) + 
      #     log(d/(bet+n0[i])^(y0[i]+alph) + c/(bet+ntx+n0[i])^(y0[i]+alph))
      # )
      ll <- ll + m*log(c/(c+d)^2) + m*alph*log(bet) - m*lgamma(alph) + sum(
        lgamma(y0[i]+alph) + y0[i]*log(n0[i]) - lgamma(y0[i]+1) - (y0[i]+alph)*log(bet+n0[i]) + 
          log(d + c/(1+ntx/(bet+n0[i]))^(y0[i]+alph))
      )
      
    }
    m <- length(i <- which(y1>0 & y0==0))
    if(m>0){
      ntx <- n1[i]*exp(tau)
      # ll <- ll + m*log(c/(c+d)^2) + m*alph*log(bet) - m*lgamma(alph) + sum(
      #   lgamma(y1[i]+alph) + y1[i]*log(ntx) - lgamma(y1[i]+1) + 
      #     log(d/(bet+ntx)^(y1[i]+alph) + c/(bet+ntx+n0[i])^(y1[i]+alph))
      # )
      ll <- ll + m*log(c/(c+d)^2) + m*alph*log(bet) - m*lgamma(alph) + sum(
        lgamma(y1[i]+alph) + y1[i]*log(ntx) - lgamma(y1[i]+1) - (y1[i]+alph)*log(bet+ntx) + 
          log(d + c/(1+n0[i]/(bet+ntx))^(y1[i]+alph))
      )
    }
    return(-ll)
  }
  
  nll_ZIPR <- function(para){ 
    alph <- para[['alph']];  bet <- para[['bet']];  mu <- para[['mu']];  sig <- para[['sig']];  c <- para[['c']];  d <- para[['d']]
    lik <- rep(NA, n)
    
    m <- length(i0 <- which(y1>0 & y0>0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        fn <- function(t0){
          ntx <- n1[i]*exp(t0)
          ll <- 2*log(c/(c+d)) + alph*log(bet) - lgamma(alph) + lgamma(y1[i]+y0[i]+alph) - (y1[i]+y0[i]+alph)*log(bet+ntx+n0[i]) + y1[i]*log(ntx) - lgamma(y1[i]+1) + y0[i]*log(n0[i]) - lgamma(y0[i]+1)
          ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
          return(exp(ll))
        }
        lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
      }
    }
    
    m <- length(i0 <- which(y1==0 & y0==0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        L0 <- d^2 + c*d/(1+n0[i]/bet)^alph
        fn <- function(t0){
          ntx <- n1[i]*exp(t0)
          L <- c*d/(1+ntx/bet)^alph + c^2/(1+(n0[i]+ntx)/bet)^alph
          L <- L/(sqrt(2*pi)*sig)*exp(- 0.5*((t0-mu)/sig)^2)
          return(L)
        }
        lik[i] <- ( L0 + ghermite.h.quadrature(fn, rule, weighted=FALSE) )/(c+d)^2
      }
    }
    
    m <- length(i0 <- which(y1==0 & y0>0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        const <- log(c/(c+d)^2) + alph*log(bet) - lgamma(alph) +  lgamma(y0[i]+alph) + y0[i]*log(n0[i]) - lgamma(y0[i]+1)
        
        fn1 <- function(t0){
          ll <- const - (y0[i]+alph)*log(bet+n0[i])
          ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
          return(d*exp(ll))
        }
        fn2 <- function(t0){
          ntx <- n1[i]*exp(t0)
          ll <- const - (y0[i]+alph)*log(bet+n0[i]+ntx)
          ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
          return(c*exp(ll))
        }
        lik[i] <- ghermite.h.quadrature(fn1, rule, weighted=FALSE) +
          ghermite.h.quadrature(fn2, rule, weighted=FALSE)
      }
    }
    
    m <- length(i0 <- which(y1>0 & y0==0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        fn <- function(t0){
          ntx <- n1[i]*exp(t0)
          ll <- log(c/(c+d)^2) + alph*log(bet) - lgamma(alph) +  lgamma(y1[i]+alph) - (y1[i]+alph)*log(bet+ntx) + y1[i]*log(ntx) - lgamma(y1[i]+1) + log(d + c/(1+n0[i]/(bet+ntx))^(y1[i]+alph))
          ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
          return(exp(ll))
        }
        lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        
      }
    }
    return(-sum(log(lik)))
  }
  
  summaryFit <- function(fit, method){
    if(method %in% c('POIF','ZIPF')){
      logRisk <- fit$par[['tau']]
      se <- sqrt(diag(solve(fit$hessian)))[['tau']]
    }
    if(method %in% c('POIR','ZIPR')){
      logRisk <- fit$par[['mu']]
      se <- sqrt(diag(solve(fit$hessian)))[['mu']]
      # se <- fit$par[['sig']] #sqrt(diag(solve(-fit$hessian)))[['mu']]
    }
    mgn <- qnorm(1-(1-conf.level)/2)*se
    bds <- exp(logRisk+c(-1,1)*mgn)
    out <- c('Log Likelihood'=-fit$value, 'Estimate'=logRisk, 'SE'=se, '95% CI lower'=bds[1], 'upper'=bds[2], 'npara'=length(fit$par), 'AIC'=2*length(fit$par)+2*fit$value, 'BIC'=log(n)*length(fit$par)+2*fit$value)
    paraEst <- fit$par
    return(list(out=out, paraEst=paraEst))
  }
  
  nmed <- length(method)
  out <- matrix(NA, nmed, 8);  paraEst <- matrix(NA, nmed, 6)
  
  # fit a base model: 'POIF'
  para <- c(alph=1, bet=1, tau=0) 
  if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']])
  m_POIF <- optim(para, nll_POIF, method='L-BFGS-B', lower=c(eps,eps,-Inf), upper=c(Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #minimize negative-log-likelihood
  
  if('POIF' %in% method){
    s_POIF <- summaryFit(m_POIF, 'POIF')
    k <- which(method=='POIF');  out[k,] <- s_POIF$out
    paraEst[k, 1:length(s_POIF$paraEst)] <- s_POIF$paraEst
  }
  
  if('POIR' %in% method || 'ZIPR' %in% method){ #use initial value from 'POIF'
    para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], mu=m_POIF$par[['tau']], sig=1e-4)
    # para <- c(alph=1.440586519511588, bet=384.857870764441486, mu=0.284836966884400, sig=0.000000005240601)
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']])
    # para[c(1,2,4)] <- log(para[c(1,2,4)])
    # m_POIR <- optim(para, nll_POIR, hessian=TRUE, control=list(trace=0)) #, maxit=500
    m_POIR <- optim(para, nll_POIR, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps), upper=c(Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0))
    # m_POIR$par[c(1,2,4)] <- exp(m_POIR$par[c(1,2,4)])
    s_POIR <- summaryFit(m_POIR, 'POIR')
    k <- which(method=='POIR')
    if(length(k)){ out[k,] <- s_POIR$out;  paraEst[k, 1:length(s_POIR$paraEst)] <- s_POIR$paraEst } 
  }
  
  if('ZIPF' %in% method){ #use initial value from 'POIF'
    para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], tau=m_POIF$par[['tau']], c=1, d=1) 
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']], c=para0[['c']], d=para0[['d']])
    m_ZIPF <- optim(para, nll_ZIPF, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1
    s_ZIPF <- summaryFit(m_ZIPF, 'ZIPF')
    k <- which(method=='ZIPF');  out[k,] <- s_ZIPF$out;  paraEst[k,c(1,2,3,5,6)] <- s_ZIPF$paraEst
  }
  
  if('ZIPR' %in% method){ #use initial value from 'POIR'
    # para <- c(alph=m_POIR$par[['alph']], bet=m_POIR$par[['bet']], mu=m_POIR$par[['mu']], sig=m_POIR$par[['sig']], c=1, d=0.01) 
    para <- c(alph=m_POIR$par[['alph']], bet=m_POIR$par[['bet']], mu=m_POIR$par[['mu']], sig=1e-2, c=1, d=1) 
    # inds <- c(1,2,4,5,6)
    # inds <- c(1,2,4)
    # para[inds] <- log(para[inds])
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']], c=para0[['c']], d=para0[['d']])
    m_ZIPR <- optim(para, nll_ZIPR, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1,
    # m_ZIPR <- optim(para, nll_ZIPR, hessian=TRUE, control=list(trace=0, maxit=5000))  #, maxit=5000
    # m_ZIPR$par[inds] <- exp(m_ZIPR$par[inds])
    s_ZIPR <- summaryFit(m_ZIPR, 'ZIPR')
    k <- which(method=='ZIPR');  out[k,] <- s_ZIPR$out;  paraEst[k,] <- s_ZIPR$paraEst
  }
  
  out <- as.data.frame(out);  row.names(out) <- method
  names(out) <- c('Log Likelihood','Estimate','SE','95% CI lower','upper','npara','AIC','BIC')
  paraEst <- as.data.frame(paraEst);  row.names(paraEst) <- method
  names(paraEst) <- c('alpha','beta','tau(mu)','sigma','c','d')
  return(list(out=out, paraEst=paraEst))
  
}


simZIPr <- function(para0, zeroInflated=TRUE, randomTau=TRUE, ni0=NULL, x0=NULL, n=NULL, seeds=123, ratio_t_c=c(0.5, 0.5)){
  #simulate data using parameter para0
  simN <- FALSE;  if(is.null(n)) n <- length(ni0) else simN <- TRUE  #default sample size using originally observed data
  set.seed(seeds)
  alph <- para0[['alpha']];  bet <- para0[['beta']]
  if(randomTau){ mu <- para0[['tau(mu)']];  sig <- para0[['sigma']] } else tau <- para0[['tau(mu)']]
  if(zeroInflated){ c <- para0[['c']];  d <- para0[['d']] } 
  ys <- rep(NA, n)
  xs <- x0;  if(is.null(x0)) xs <- sample(c(0,1), n, replace=TRUE, prob=ratio_t_c)
  nis <- ni0;  if(simN)  nis <- round(runif(n, 20, 100))  #simulate study size
  # if(simN)  nis <- round(sample(ni0, n, replace=TRUE))  #resampling from the observed study size
  pis <- rep(1,n);  ai <- rep(NA, n)
  if(zeroInflated) pis <- rbeta(n, c, d);  ai <- rep(NA, n)
  for(i in 1:n){
    xi <- rgamma(1, shape=alph, rate=bet)
    if(randomTau) tauUse <- rnorm(1, mu, sig)  else tauUse <- tau
    lambda <- xi*exp(tauUse*xs[i])
    ai[i] <- rbinom(1,1,pis[i]); if(ai[i]==1) ys[i] <- rpois(1, nis[i]*lambda[1]) else ys[i] <- 0
  }
  return(data.frame(y=ys, x=xs, ni=nis))
}



metaZIPr0 <- function(y1, n1, y0, n0, method=c('POIF','POIR','ZIPF','ZIPR'), para0=NULL, conf.level=.95, eps=.Machine$double.eps){ 
  #meta-ZIP model with randomized trials on 2 groups (x=0,1): e.g. ROS in the manuscript
  #each study has y1 events out of n1 (treatment), and y0 events out of n0 (control).
  #fixed baseline xi instead of xi~Gamma(alph, bet)
  n <- length(y1)  #number of studies
  
  # negative log-likelihood functions
  nll_POIF <- function(para){  #para=c(alph, bet, tau)
    xi <- para[['xi']];  tau <- para[['tau']]
    ntx <- n1*exp(tau)
    # val <- n*(alph*log(bet)-lgamma(alph)) + sum( lgamma(y1+y0+alph) - (y1+y0+alph)*log(bet+ntx+n0) + y1*log(ntx) - lgamma(y1+1) + y0*log(n0) - lgamma(y0+1))  #log-likelihood
    val <- sum( y1*log(n1*xi) - lgamma(y1+1) + y0*log(n0*xi) - lgamma(y0+1) + tau*y1 - xi*(n0 + n1*exp(tau)))  #log-likelihood
    # val <- sum(-n1*xi*exp(tau) + y1*log(n1*xi) + y1*tau - lgamma(y1+1) - n0*xi + y0*log(n0*xi) - lgamma(y0+1))
    # val <- sum(dpois(y1, n1*xi*exp(tau), log=TRUE) + dpois(y0, n0*xi, log=TRUE))
    return(-val)
  }
  
  nll_POIR <- function(para){  #para=c(alph, bet, mu, lsig)
    xi <- para[['xi']];  mu <- para[['mu']];  sig <- para[['sig']]
    lik <- rep(NA, n)
    # const <- lgamma(y1+y0+alph) - lgamma(y1+1) - lgamma(y0+1) + y1*log(n1) + y0*log(n0) + alph*log(bet) - lgamma(alph) - 0.5*log(pi) + 0.5*((mu/sig+y1*sig)^2-(mu/sig)^2)
    for(i in 1:n){
      fn <- function(t0){
        ll <- y1[i]*log(n1[i]*xi) - lgamma(y1[i]+1) + y0[i]*log(n0[i]*xi) - lgamma(y0[i]+1) - xi*(n0[i]+n1[i]*exp(t0)) + t0*y1[i]
        ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
        return(exp(ll))
      }
      lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
      # if(is.nan(lik[i]) || lik[i]==0){  #this may happen numerically in extreme case
      #   const <- lgamma(y1[i]+y0[i]+alph) - lgamma(y1[i]+1) - lgamma(y0[i]+1) + y1[i]*log(n1[i]) + y0[i]*log(n0[i]) + alph*log(bet) - lgamma(alph) - 0.5*log(pi) + 0.5*((mu/sig+y1[i]*sig)^2-(mu/sig)^2)
      #   lik[i] <- sum(exp( const + log(gauss.15$weights) - (y0[i]+y1[i]+alph)*log(bet+n1[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y1[i]*sig^2)+n0[i]) ))
      # }
    }
    return(-sum(log(lik)))
  }
  
  
  
  nll_ZIPF <- function(para){ 
    # xi <- para[['xi']];  tau <- para[['tau']];  c <- para[['c']];  d <- para[['d']]
    xi <- exp(para[['xi']]);  tau <- para[['tau']];  c <- exp(para[['c']]);  d <- exp(para[['d']])
    ll <- 0
    m <- length(i <- which(y1>0 & y0>0))
    if(m>0){
      lambda0 <- n0[i]*xi
      lambda1 <- n1[i]*xi*exp(tau)
      ll <- ll + m*log((c+1)*c/(c+d+1)/(c+d)) + sum(
        -lambda1 + y1[i]*log(lambda1) - lgamma(y1[i]+1) - lambda0 + y0[i]*log(lambda0) - lgamma(y0[i]+1)
      )
    }
    m <- length(i <- which(y1==0 & y0==0))
    if(m>0){
      p0 <- exp(-n0[i]*xi)
      p1 <- exp(-n1[i]*xi*exp(tau))
      ll <- ll - m*log((c+d+1)*(c+d)) + sum(log( (d+1)*d + c*d*(p0+p1) + (c+1)*c*p0*p1 ))
    }
    m <- length(i <- which(y1==0 & y0>0))
    if(m>0){
      p1 <- exp(-n1[i]*xi*exp(tau))
      ll <- ll + m*log(c/(c+d+1)/(c+d)) + sum(-n0[i]*xi + y0[i]*log(n0[i]*xi) - lgamma(y0[i]+1) + log(d + (c+1)*p1))
    }
    m <- length(i <- which(y1>0 & y0==0))
    if(m>0){
      lamnda1 <- n1[i]*xi*exp(tau)
      logP1 <- -lamnda1 + y1[i]*log(lamnda1) - lgamma(y1[i]+1)
      ll <- ll + m*log(c/(c+d+1)/(c+d)) + sum(logP1 - lgamma(y0[i]+1) + log(d + (c+1)*p0))
    }
    return(-ll)
  }
  
  
  
  nll_ZIPR <- function(para){ 
    # xi <- para[['xi']];  mu <- para[['mu']];  sig <- para[['sig']];  c <- para[['c']];  d <- para[['d']]
    xi <- exp(para[['xi']]);  mu <- para[['mu']];  sig <- exp(para[['sig']]);  c <- exp(para[['c']]);  d <- exp(para[['d']])
    # lik <- rep(NA, n)
    ll <- 0
    m <- length(i0 <- which(y1>0 & y0>0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        #   fn <- function(t0){
        #     lambda0 <- n0[i]*xi
        #     lambda1 <- n1[i]*xi*exp(t0)
        #     ll <- log((c+1)*c/(c+d+1)/(c+d)) + y1[i]*log(lambda1) - lambda1 - lgamma(y1[i]+1) + y0[i]*log(lambda0) - lgamma(y0[i]+1) - lambda0
        #     ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
        #     return(exp(ll))
        #   }
        #   lik[i] <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        lambda0 <- n0[i]*xi
        k <- -lambda0 + y0[i]*log(lambda0) - lgamma(y0[i]+1) + y1[i]*log(n1[i]*xi) - lgamma(y1[i]+1) - 0.5*log(pi) + 0.5*((mu/sig+y1[i]*sig)^2-(mu/sig)^2)
        p01 <- sum(exp( k + log(gauss.15$weights) - n1[i]*xi*exp(sqrt(2)*sig*gauss.15$nodes+mu+sig^2*y1[i]) ))
        ll <- ll + log((c+1)*c/(c+d+1)/(c+d)) + log(p01)
      }
    }
    
    m <- length(i0 <- which(y1==0 & y0==0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        p0 <- exp(-n0[i]*xi)
        k <- - 0.5*log(pi)
        p1 <- sum(exp( k + log(gauss.15$weights) - n1[i]*xi*exp(sqrt(2)*sig*gauss.15$nodes+mu) ))
        ll <- ll - log((c+d+1)*(c+d)) + log( (d+1)*d + c*d*(p0+p1) + (c+1)*c*p0*p1 )
      }
    }
    
    m <- length(i0 <- which(y1==0 & y0>0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        p1 <- sum(exp( - 0.5*log(pi) + log(gauss.15$weights) - ni[i]*xi*exp(sqrt(2)*sig*gauss.15$nodes+mu) ))
        if(p1==0){
          fn <- function(t0){
            lambda <- xi*exp(t0)
            ll <- -ni[i]*lambda + dnorm(t0, mu, sig, log=TRUE)
            return(exp(ll))
          }
          p1 <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        }
        ll <- ll + log(c/(c+d+1)/(c+d)) - n0[i]*xi + y0[i]*log(n0[i]*xi) - lgamma(y0[i]+1) + log(d + (c+1)*p1)
      }
    }
    
    m <- length(i0 <- which(y1>0 & y0==0))
    if(m>0){
      for(j in 1:m){
        i <- i0[j]
        p0 <- exp(-n0[i]*xi)
        p1 <- sum(exp(  y1[i]*log(n1[i]*xi) - lgamma(y1[i]+1) - 0.5*log(pi) + 0.5*((mu/sig+y1[i]*sig)^2-(mu/sig)^2) + log(gauss.15$weights) - ni[i]*xi*exp(sqrt(2)*sig*gauss.15$nodes+mu+sig^2*y1[i]) ))
        # if(p1==0){
        #   fn <- function(t0){
        #     lambda <- xi*exp(t0)
        #     ll <- -ni[i]*lambda + dnorm(t0, mu, sig, log=TRUE)
        #     return(exp(ll))
        #   }
        #   p1 <- ghermite.h.quadrature(fn, rule, weighted=FALSE)
        # }
        ll <- ll + log(c/(c+d+1)/(c+d)) + log(d*p1 + (c+1)*p0*p1)
      }
    }
    return(-ll)
  }
  
  summaryFit <- function(fit, method){
    if(method %in% c('POIF','ZIPF')){
      logRisk <- fit$par[['tau']]
      se <- sqrt(diag(solve(fit$hessian)))[['tau']]
    }
    if(method %in% c('POIR','ZIPR')){
      logRisk <- fit$par[['mu']]
      se <- sqrt(diag(solve(fit$hessian)))[['mu']]
      # se <- fit$par[['sig']] #sqrt(diag(solve(-fit$hessian)))[['mu']]
    }
    mgn <- qnorm(1-(1-conf.level)/2)*se
    bds <- exp(logRisk+c(-1,1)*mgn)
    out <- c('Log Likelihood'=-fit$value, 'Estimate'=logRisk, 'SE'=se, '95% CI lower'=bds[1], 'upper'=bds[2], 'npara'=length(fit$par), 'AIC'=2*length(fit$par)+2*fit$value, 'BIC'=log(n)*length(fit$par)+2*fit$value)
    paraEst <- fit$par
    return(list(out=out, paraEst=paraEst))
  }
  
  nmed <- length(method)
  out <- matrix(NA, nmed, 8);  paraEst <- matrix(NA, nmed, 5)
  
  # fit a base model: 'POIF'
  para <- c(xi=0.003745472, tau=0.25) 
  # if(!is.null(para0)) para <- c(alph=para0[['xi']], tau=para0[['tau']])
  m_POIF <- optim(para, nll_POIF, method='L-BFGS-B', lower=c(eps,-Inf), upper=c(Inf,Inf), hessian=TRUE, control=list(trace=0)) #minimize negative-log-likelihood
  
  if('POIF' %in% method){
    s_POIF <- summaryFit(m_POIF, 'POIF')
    k <- which(method=='POIF');  out[k,] <- s_POIF$out
    paraEst[k, 1:length(s_POIF$paraEst)] <- s_POIF$paraEst
  }
  
  if('POIR' %in% method || 'ZIPR' %in% method){ #use initial value from 'POIF'
    para <- c(xi=m_POIF$par[['xi']], mu=m_POIF$par[['tau']], sig=1)
    # para <- c(alph=1.440586519511588, bet=384.857870764441486, mu=0.284836966884400, sig=0.000000005240601)
    # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']])
    # para[c(1,2,4)] <- log(para[c(1,2,4)])
    # m_POIR <- optim(para, nll_POIR, hessian=TRUE, control=list(trace=0)) #, maxit=500
    m_POIR <- optim(para, nll_POIR, method='L-BFGS-B', lower=c(eps,-Inf,eps), upper=c(Inf,Inf,Inf), hessian=TRUE, control=list(trace=0))
    # m_POIR$par[c(1,2,4)] <- exp(m_POIR$par[c(1,2,4)])
    s_POIR <- summaryFit(m_POIR, 'POIR')
    k <- which(method=='POIR')
    if(length(k)){ out[k,] <- s_POIR$out;  paraEst[k, 1:length(s_POIR$paraEst)] <- s_POIR$paraEst } 
  }
  
  if('ZIPF' %in% method){ #use initial value from 'POIF'
    para <- c(xi=m_POIF$par[['xi']], tau=m_POIF$par[['tau']], c=1, d=1) 
    # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']], c=para0[['c']], d=para0[['d']])
    para[c(1,3,4)] <- log(para[c(1,3,4)])
    # m_ZIPF <- optim(para, nll_ZIPF, method='L-BFGS-B', lower=c(eps,-Inf,eps,eps), upper=c(Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1
    m_ZIPF <- optim(para, nll_ZIPF, hessian=TRUE)
    m_ZIPF$par[c(1,3,4)] <- exp(m_ZIPF$par[c(1,3,4)])
    s_ZIPF <- summaryFit(m_ZIPF, 'ZIPF')
    k <- which(method=='ZIPF');  out[k,] <- s_ZIPF$out;  paraEst[k,c(1,2,4,5)] <- s_ZIPF$paraEst
  }
  
  if('ZIPR' %in% method){ #use initial value from 'POIR'
    # para <- c(alph=m_POIR$par[['alph']], bet=m_POIR$par[['bet']], mu=m_POIR$par[['mu']], sig=m_POIR$par[['sig']], c=1, d=0.01) 
    para <- c(xi=m_POIR$par[['xi']], mu=m_POIR$par[['mu']], sig=0.02, c=1, d=1) 
    inds <- c(1,3,4,5)
    # inds <- c(1,2,4)
    para[inds] <- log(para[inds])
    # if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']], c=para0[['c']], d=para0[['d']])
    # m_ZIPR <- optim(para, nll_ZIPR, method='L-BFGS-B', lower=c(eps,-Inf,eps,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0)) #maximize log-likelihood fnscale=-1,
    m_ZIPR <- optim(para, nll_ZIPR, hessian=TRUE, control=list(trace=0, maxit=5000))  #, maxit=5000
    m_ZIPR$par[inds] <- exp(m_ZIPR$par[inds])
    s_ZIPR <- summaryFit(m_ZIPR, 'ZIPR')
    k <- which(method=='ZIPR');  out[k,] <- s_ZIPR$out;  paraEst[k,] <- s_ZIPR$paraEst
  }
  
  out <- as.data.frame(out);  row.names(out) <- method
  names(out) <- c('Log Likelihood','Estimate','SE','95% CI lower','upper','npara','AIC','BIC')
  paraEst <- as.data.frame(paraEst);  row.names(paraEst) <- method
  names(paraEst) <- c('xi','tau(mu)','sigma','c','d')
  return(list(out=out, paraEst=paraEst))
  
}


# not run


