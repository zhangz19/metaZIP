

#+++++++++++++++++++++++++++++ ROS data: multiple randomized studies on the same drug
rm(list=ls())
source('util.R')

# real data
dat <- read.table("midata.txt", header=T, quote="\"") #48 obs.
n <- nrow(dat)  #total number of studies
n1 <- dat$Rosiglitazone.SampSize  #size of treatment group (Rosiglitazone)
n0 <- dat$Placebo.SampSize  #size of placebo group (Rosiglitazone)
y1 <- dat$Rosiglitazone.MI  #number of mi cases in treatment group (Rosiglitazone)
y0 <- dat$Placebo.MI  #number of mi cases in placebo group

# dat <- read.csv("RosiglitazoneData.csv") #56 obs.
# n <- nrow(dat)  #total number of studies
# n1 <- dat$n_rosi  #size of treatment group (Rosiglitazone)
# n0 <- dat$n_ctl   #size of placebo group (Rosiglitazone)
# y1 <- dat$MI_rosi  #number of mi cases in treatment group (Rosiglitazone)
# y0 <- dat$MI_ctl  #number of mi cases in placebo group

metaZIPr(y1=y1, n1=n1, y0=y0, n0=n0, method=c('POIF','POIR','ZIPF','ZIPR'), eps=0.01)

print(fit <- metaZIPr(y1=y1, n1=n1, y0=y0, n0=n0, eps=0.01))
u <- fit$out[['Log Likelihood']]

# check results using Binomial-logit regression
require(lme4)
y <- c(y0,y1);  ni <- c(n0,n1);  study <- c(1:length(y0), 1:length(y1));  x <- c(rep(0,length(y0)), rep(1,length(y1)))
summary(gm1 <- glmer(cbind(y, ni - y) ~ (1 | study) + x, family=binomial(link = "logit"), na.action=na.omit))
logLik(gm1)
AIC(gm1) #142.261
-2*logLik(gm1)+2*(1+length(fixef(gm1)))  #add one variation parameter for study

metaZIPr0(y1=y1, n1=n1, y0=y0, n0=n0, method=c('POIF','POIR','ZIPF','ZIPR'))







simZIPr <- function(para0, zeroInflated=TRUE, randomTau=TRUE, n1s=NULL, n0s=NULL, n=NULL, seeds=123){
  #simulate data using parameter para0
  simN <- FALSE;  if(is.null(n)) n <- length(n1s) else simN <- TRUE  #default sample size using originally observed data
  set.seed(seeds)
  alph <- para0[['alpha']];  bet <- para0[['beta']]
  if(randomTau){ mu <- para0[['tau(mu)']];  sig <- para0[['sigma']] } else tau <- para0[['tau(mu)']]
  if(zeroInflated){ c <- para0[['c']];  d <- para0[['d']] } 
  n1 <- n1s;  n0 <- n0s
  if(simN){
    n1 <- n0 <- rep(NA, n)
    n0 <- round(runif(n, 20, 50))
    n1 <- n0*sample(c(1,2,3), n, replace=TRUE, prob=c(0.4,0.4,0.2))
  }
  y1 <- y0 <- rep(NA, n)
  for(i in 1:n){
    if(zeroInflated) pis <- rbeta(2, c, d)  else  pis <- rep(1,2)  #zero-inflated or not
    xi <- rgamma(1, shape=alph, rate=bet)
    if(randomTau) tauUse <- rnorm(1, mu, sig)  else tauUse <- tau
    lambda <- xi*exp(tauUse*c(0,1))
    ai <- rbinom(1,1,pis[1]); if(ai==1) y0[i] <- rpois(1, n0[i]*lambda[1]) else y0[i] <- 0
    ai <- rbinom(1,1,pis[2]); if(ai==1) y1[i] <- rpois(1, n1[i]*lambda[2]) else y1[i] <- 0
  }
  return(data.frame(y1=y1, n1=n1, y0=y0, n0=n0))
}

para0 <- fit$paraEst[4,]  #c(alph=3, bet=30, tau=.3, mu=.3, sig=.1, c=2, d=5)
sdat <- simZIPr(para0=para0, zeroInflated=T, randomTau=T, n1s=n1, n0s=n0, seeds=123)
m <- metaZIPr(y1=sdat$y1, n1=sdat$n1, y0=sdat$y0, n0=sdat$n0, method=c('POIF','POIR','ZIPF','ZIPR'), eps=0.01)


# simulation
set.seed(123)
simN <- TRUE
# para0 <- c(alph=3, bet=30, tau=3, mu=3, sig=0.01, c=2, d=10);  ZIP <- TRUE; re <- FALSE #ZIPF
para0 <- c(alph=3, bet=30, tau=3, mu=3, sig=1, c=10, d=10);  ZIP <- TRUE; re <- TRUE #ZIPR
alph <- para0[['alph']];  bet <- para0[['bet']];  mu <- para0[['mu']]; sig <- para0[['sig']];  c <- para0[['c']];  d <- para0[['d']]
if(simN) n <- 100 #2e3 #large sample for better estimation
y1 <- y0 <- rep(NA, n)
if(simN){
  n1 <- n0 <- rep(NA, n)
  n0 <- round(runif(n, 20, 50))
  n1 <- n0*sample(c(1,2,3), n, replace=TRUE, prob=c(0.4,0.4,0.2))
}
for(i in 1:n){
  if(ZIP) pis <- rbeta(2, c, d)  else  pis <- rep(1,2)  #zero-inflated or not
  xi <- rgamma(1, shape=alph, rate=bet)
  if(re) tau <- rnorm(1, mu, sig) else tau <- mu #random effect or fixed parameter
  lambda <- xi*exp(tau*c(0,1))
  ai <- rbinom(1,1,pis[1]); if(ai==1) y0[i] <- rpois(1, n0[i]*lambda[1]) else y0[i] <- 0
  ai <- rbinom(1,1,pis[2]); if(ai==1) y1[i] <- rpois(1, n1[i]*lambda[2]) else y1[i] <- 0
}

tmp <- metaZIPr(y1=y1, n1=n1, y0=y0, n0=n0, para0=para0, method=c('POIF','POIR','ZIPF','ZIPR'))
print(tmp$out); rbind(tmp$paraEst, para0[-3])





#+++++++++++++++++++++++++++++ EUS-FNA data: multiple studies, each on one of the 2 drugs
rm(list=ls())
source('util.R')

dat <- read.table('EUS-FNA.txt', h=F)  #51 obs.
inds <- unlist(gregexpr('/', dat$V2))
x <- 2 - dat$V1
# x <-  dat$V1-1
y <- as.numeric(substr(dat$V2, 1,    inds-1))  #y cases out of ni
ni <- as.numeric(substr(dat$V2, inds+1, 100L))

# real data
print(fit <- metaZIPs(y=y, x=x, ni=ni))
print(fit <- metaZIPs0(y=y, x=x, ni=ni, method=c('POIF','POIR','ZIPF','ZIPR')))


u <- fit$out[['Log Likelihood']]
# check results using Binomial-logit regression
require(lme4)
study <- 1:length(y)
summary(gm1 <- glmer(cbind(y, ni - y) ~ (1 | study) + x, family=binomial(link = "logit"), na.action=na.omit))
logLik(gm1)
AIC(gm1) #142.261
-2*logLik(gm1)+2*(1+length(fixef(gm1)))  #add one variation parameter for study

# model comparison
nsim <- 300  #more runs than target, since the optimization is numerically unstable
nmodel <- length(models <- c('POIF','POIR','ZIPF','ZIPR'))
modelComp <- numeric()
verbose <- FALSE
for(i in 1:nmodel){
  mat <- matrix(NA, nsim, nmodel)  # 4 models: 'POIF','POIR','ZIPF','ZIPR'
  for(j in 1:nsim){ 
    if(verbose) {cat(j,''); if(!j%%20) cat('\n')}
    if(models[i]=='POIF') {zeroInflated=FALSE; randomTau=FALSE}
    if(models[i]=='POIR') {zeroInflated=FALSE; randomTau=TRUE}
    if(models[i]=='ZIPF') {zeroInflated=TRUE; randomTau=FALSE}
    if(models[i]=='ZIPR') {zeroInflated=TRUE; randomTau=TRUE}
    sdat <- simZIPs(para0=fit$paraEst[i,], zeroInflated=zeroInflated, randomTau=randomTau, ni0=ni, x0=x, seeds=j)
    if(TRUE){ #if re-estimate the parameters, would be using a form of parametric boostrap (Efron & Tibshirani, 1993)
      tmp <- try(metaZIPs(y=sdat$y, x=sdat$x, ni=sdat$ni, method=models)$out[['Log Likelihood']], silent=TRUE)
      if(class(tmp)!='try-error') mat[j,] <- tmp
    }
  }
  a <- na.omit(mat); if(nrow(a)>200) a <- a[1:200,]
  # calculate the Mahalanobis distances for general model comparisons
  ui <- colMeans(a); Si <- cov(a)
  D2 <- t(u-ui)%*%solve(Si)%*%(u-ui)
  modelComp <- rbind(modelComp,  data.frame(nsim=nrow(a), D2=D2,  pval=1-pf(D2/4, 4, nrow(a)-1)))
}
row.names(modelComp) <- models
print(modelComp)


#simulation to check code
set.seed(1234)
simN <- TRUE  #if FALSE, will use the sample size in the real data
para0 <- c(alph=5, bet=200, tau=1, mu=1, sig=.4, c=5, d=5);  ZIP <- TRUE;  re <- TRUE
#if too extreme, e.g. c=0.2, d=10, estimated alph, bet for the Poisson part may deviate from the truth
alph <- para0[['alph']];  bet <- para0[['bet']];  mu <- para0[['mu']];  sig <- para0[['sig']]; c <- para0[['c']];  d <- para0[['d']]
if(simN) n <- 1e3  else n <- nrow(read.table('EUS-FNA.txt', h=F))
# if(simN) n <- 5e2  else n <- nrow(read.table('EUS-FNA.txt', h=F))
ys <- rep(NA, n)
xs <- sample(c(0,1), n, replace=TRUE, prob=c(0.7, 0.3))
nis <- ni;  if(simN)  nis <- round(runif(n, 20, 100))
# if(simN)  nis <- round(sample(ni, n, replace=TRUE))
pis <- rep(1,n);  ai <- rep(NA, n)
if(ZIP) pis <- rbeta(n, c, d);  ai <- rep(NA, n)
for(i in 1:n){
  xi <- rgamma(1, shape=alph, rate=bet)
  if(re) tau <- rnorm(1, mu, sig)  else tau <- mu
  lambda <- xi*exp(tau*xs[i])
  ai[i] <- rbinom(1,1,pis[i]); if(ai[i]==1) ys[i] <- rpois(1, nis[i]*lambda[1]) else ys[i] <- 0
}
all(ys<nis)
range(ys)
tmp <- metaZIPs(y=ys, x=xs, ni=nis, method=c('POIF','POIR','ZIPF','ZIPR'), para0=para0)
print(tmp$out); rbind(tmp$paraEst, para0[-3])



# not run


