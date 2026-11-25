# Historical single-arm calculation from Zhen Zhang's original metaZIPs().
# Likelihoods, starting values, optimizers, and quadrature orders are retained.
# The wrapper below labels historical spread/limits and parameter counts.
# A failed inverse Hessian returns NA rather than substituting another spread.

.metaZIPs_legacy_fit <- function(y, x, ni, method=c('POIF','POIR','ZIPF','ZIPR'), para0=NULL, conf.level=.95, eps=.Machine$double.eps){
    rule <- gaussquad::ghermite.h.quadrature.rules(21, mu = 0, normalized = FALSE)[[21]]
    gauss.15 <- statmod::gauss.quad(15, "hermite")
    n <- length(ni)

    nll_POIF <- function(para){
        alph <- para[['alph']];  bet <- para[['bet']];  tau <- para[['tau']]
        ntx <- ni*exp(tau*x)
        val <- n*(alph*log(bet)-lgamma(alph)) + sum( lgamma(y+alph) - (y+alph)*log(bet+ntx) + y*log(ntx) - lgamma(y+1))
        return(-val)
    }

    nll_POIR <- function(para){
        alph <- exp(para[['alph']]);  bet <- exp(para[['bet']]);  mu <- para[['mu']];  sig <- exp(para[['sig']])
        lik <- rep(NA, n)
        for(i in 1:n){
            if(x[i]==1){



                k <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) + y[i]*log(ni[i]) - lgamma(y[i]+1) - 0.5*log(pi) + 0.5*( (mu/sig+y[i]*sig)^2-(mu/sig)^2 )
                lik[i] <- sum(exp( k + log(gauss.15$weights) - (y[i]+alph)*log(bet+ni[i]*exp(sqrt(2)*sig*gauss.15$nodes+mu+y[i]*sig^2)) ))

                if(lik[i]==0){
                    fn <- function(t0){
                        ntx <- ni[i]*exp(t0*x[i])
                        ll <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) - (y[i]+alph)*log(bet+ntx) +
                            y[i]*log(ntx) - lgamma(y[i]+1)
                        ll <- ll - 0.5*log(2*pi) - log(sig) - 0.5*((t0-mu)/sig)^2
                        return(exp(ll))
                    }
                    lik[i] <- gaussquad::ghermite.h.quadrature(fn, rule, weighted=FALSE)
                }

            }
            if(x[i]==0){
                ntx <- ni[i]
                lik[i] <- alph*log(bet) + lgamma(y[i]+alph) - lgamma(alph) - (y[i]+alph)*log(bet+ntx) +
                    y[i]*log(ntx) - lgamma(y[i]+1)
                lik[i] <- exp(lik[i])
            }
        }
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
                h1i[i] <- gaussquad::ghermite.h.quadrature(function(t) exp(-t^2- alph*log(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu))), rule, weighted=FALSE)
            }
            ll <- ll + sum(log(  d/(c+d)*(1-bet^alph*h1i/sqrt(pi)) + bet^alph*h1i/sqrt(pi)  ))
        }
        n0 <- length(i0 <- which(y!=0 & x==1))
        if(n0>0){
            h2i <- rep(NA, n0)
            for(i in 1:n0){
                h2i[i] <- gaussquad::ghermite.h.quadrature(function(t) exp( lgamma(y[i0[i]]+alph) - lgamma(y[i0[i]]+1) + y[i0[i]]*log(ni[i0[i]]) - 0.5*(mu^2/sig^2 - (mu/sig + y[i0[i]]*sig)^2) -t^2- (y[i0[i]]+alph)*log(bet + ni[i0[i]]*exp(sqrt(2)*sig*t+mu+y[i0[i]]*sig^2))  ), rule, weighted=FALSE)

            }
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
            se <- tryCatch(sqrt(diag(solve(fit$hessian)))[['tau']],
                           error = function(e) NA_real_)
        }
        if(method %in% c('POIR','ZIPR')){
            logRisk <- fit$par[['mu']]
            se <- fit$par[['sig']]
        }
        mgn <- stats::qnorm(1-(1-conf.level)/2)*se
        bds <- exp(logRisk+c(-1,1)*mgn)
        out <- c('Log Likelihood'=-fit$value, 'Estimate'=logRisk, 'SE'=se, '95% CI lower'=bds[1], 'upper'=bds[2], 'npara'=length(fit$par), 'AIC'=2*length(fit$par)+2*fit$value, 'BIC'=log(n)*length(fit$par)+2*fit$value)
        paraEst <- fit$par
        return(list(out = out, paraEst = paraEst))
    }

    nmed <- length(method)
    out <- matrix(NA, nmed, 8);  paraEst <- matrix(NA, nmed, 6)

    para <- c(alph=1, bet=1, tau=0)
    if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']])
    m_POIF <- stats::optim(para, nll_POIF, method='L-BFGS-B', lower=c(eps,eps,-Inf), upper=c(Inf,Inf,Inf), hessian=TRUE, control=list(trace=0))

    if('POIF' %in% method){
        s_POIF <- summaryFit(m_POIF, 'POIF')
        k <- which(method=='POIF');  out[k,] <- s_POIF$out
        paraEst[k, 1:length(s_POIF$paraEst)] <- s_POIF$paraEst
    }

    if('POIR' %in% method || 'ZIPR' %in% method){
        para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], mu=m_POIF$par[['tau']], sig=1)
        if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']])
        para[c(1,2,4)] <- log(para[c(1,2,4)])
        m_POIR <- stats::optim(para, nll_POIR, hessian=TRUE, control=list(trace=0, maxit=500))
        m_POIR$par[c(1,2,4)] <- exp(m_POIR$par[c(1,2,4)])
        s_POIR <- summaryFit(m_POIR, 'POIR')
        k <- which(method=='POIR')
        if(length(k)){ out[k,] <- s_POIR$out;  paraEst[k, 1:length(s_POIR$paraEst)] <- s_POIR$paraEst }
    }

    if('ZIPF' %in% method){
        para <- c(alph=m_POIF$par[['alph']], bet=m_POIF$par[['bet']], tau=m_POIF$par[['tau']], c=1, d=1)
        if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], tau=para0[['tau']], c=para0[['c']], d=para0[['d']])
        m_ZIPF <- stats::optim(para, nll_ZIPF, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0))
        s_ZIPF <- summaryFit(m_ZIPF, 'ZIPF')
        k <- which(method=='ZIPF');  out[k,] <- s_ZIPF$out;  paraEst[k,c(1,2,3,5,6)] <- s_ZIPF$paraEst
    }

    if('ZIPR' %in% method){
        para <- c(alph=m_POIR$par[['alph']], bet=m_POIR$par[['bet']], mu=m_POIR$par[['mu']], sig=m_POIR$par[['sig']], c=1, d=1)
        if(!is.null(para0)) para <- c(alph=para0[['alph']], bet=para0[['bet']], mu=para0[['mu']], sig=para0[['sig']], c=para0[['c']], d=para0[['d']])
        m_ZIPR <- stats::optim(para, nll_ZIPR, method='L-BFGS-B', lower=c(eps,eps,-Inf,eps,eps,eps), upper=c(Inf,Inf,Inf,Inf,Inf,Inf), hessian=TRUE, control=list(trace=0))
        s_ZIPR <- summaryFit(m_ZIPR, 'ZIPR')
        k <- which(method=='ZIPR');  out[k,] <- s_ZIPR$out;  paraEst[k,] <- s_ZIPR$paraEst
    }

    out <- as.data.frame(out);  row.names(out) <- method
    names(out) <- c('Log Likelihood','Estimate','SE','95% CI lower','upper','npara','AIC','BIC')
    paraEst <- as.data.frame(paraEst);  row.names(paraEst) <- method
    names(paraEst) <- c('alpha','beta','tau(mu)','sigma','c','d')


    fits <- list(POIF = m_POIF)
    if (any(c("POIR", "ZIPR") %in% method)) fits$POIR <- m_POIR
    if ("ZIPF" %in% method) fits$ZIPF <- m_ZIPF
    if ("ZIPR" %in% method) fits$ZIPR <- m_ZIPR
    return(list(out = out, paraEst = paraEst, fits = fits[method]))
}

.legacy_start <- function(para0, eps) {
  if (is.null(para0)) return(NULL)
  if (is.list(para0)) para0 <- unlist(para0, use.names = TRUE)
  allowed <- c("alpha", "alph", "beta", "bet", "tau", "mu", "tau(mu)",
               "sigma", "sig", "q", "c", "d")
  if (!is.numeric(para0) || !length(para0) || is.null(names(para0)) ||
      anyNA(names(para0)) || anyDuplicated(names(para0)) ||
      any(!names(para0) %in% allowed) || any(!is.finite(para0)))
    stop("para0 must contain finite, named legacy starting values; see ?metaZIPs.",
         call. = FALSE)
  pick <- function(keys, default) {
    k <- keys[keys %in% names(para0)]
    if (length(k)) unname(para0[[k[[1L]]]]) else default
  }
  p <- c(alph = pick(c("alpha", "alph"), 1),
         bet = pick(c("beta", "bet"), 1),
         tau = pick(c("tau", "mu", "tau(mu)"), 0),
         mu = pick(c("mu", "tau(mu)", "tau"), 0),
         sig = pick(c("sigma", "sig"), 1),
         c = pick("c", 1), d = pick("d", 1))
  if ("q" %in% names(para0)) {
    q <- para0[["q"]]
    if (q <= 0 || q >= 1)
      stop("Legacy starting q must be strictly between 0 and 1.", call. = FALSE)
    p[c("c", "d")] <- c(q, 1 - q)
  }
  if (any(p[c("alph", "bet", "sig", "c", "d")] < eps))
    stop("Legacy alpha, beta, sigma, c and d starts must be at least eps (> 0).",
         call. = FALSE)
  p
}

.metaZIPs_legacy <- function(y, x, ni, method, para0, conf.level, eps,
                             control, call) {
  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!method %in% names(.model_indices)) || anyDuplicated(method))
    stop("method must contain unique values from POIF, POIR, ZIPF, ZIPR.",
         call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) ||
      eps <= 0 || eps >= 1)
    stop("Legacy eps must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  ctl <- .metaZIP_control(control)
  if (length(setdiff(names(control), "warn")))
    stop("The legacy engine retains fixed numerical settings; ",
         "control supports only warn.", call. = FALSE)
  if (any(vapply(c(0, 1), function(g) sum(y[x == g]) == 0, logical(1))))
    stop("A finite log-risk-ratio fit requires at least one event in each group overall.",
         call. = FALSE)

  result <- .metaZIPs_legacy_fit(y, x, ni, method,
                                 .legacy_start(para0, eps), conf.level, eps)
  names(result$out) <- c("Log Likelihood", "Estimate", "Legacy spread",
                         "Legacy lower", "Legacy upper", "Legacy npara",
                         "Legacy AIC", "Legacy BIC")
  spread_type <- ifelse(method %in% c("POIR", "ZIPR"),
                         "Heterogeneity SD", "Inverse-Hessian SE")
  result$out$RR <- exp(result$out$Estimate)
  result$out[["Spread type"]] <- spread_type
  optimizer_code <- vapply(result$fits, function(f) f$convergence, integer(1))
  messages <- vapply(result$fits, function(f)
    if (is.null(f$message)) "" else f$message, character(1))
  diagnostics <- data.frame(converged = optimizer_code == 0L,
                            optimizer_code = optimizer_code,
                            spread_type = spread_type,
                            spread_available = is.finite(result$out[["Legacy spread"]]),
                            message = messages, row.names = method,
                            stringsAsFactors = FALSE)
  if (ctl$warn && any(!diagnostics$converged | !diagnostics$spread_available))
    warning("Some historical fits did not converge or have unavailable spreads; ",
            "inspect $diagnostics.", call. = FALSE)
  result$diagnostics <- diagnostics
  result$design <- "single"
  result$nstudies <- length(y)
  result$data <- list(y = y, x = x, ni = ni, design = "single")
  result$conf.level <- conf.level
  result$control <- list(warn = ctl$warn, nquad_POIR = 15L, nquad_ZIPR = 21L)
  result$engine <- "legacy"
  result$call <- call
  result$estimand <- "mean log risk ratio; historical spread is not uniformly its SE"
  class(result) <- "metaZIP_legacy"
  result
}
