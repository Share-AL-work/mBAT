library(survey)
library(plyr)
library(ggplot2)
library(mvtnorm)
library(mgcv)

setwd("")

## Settings
n = 1000
beta1 = 1
beta2vec = c(-0.1, -0.3, -0.5, -0.8, -1)
LDvec = c(seq(0,0.9,0.1), 0.99)
m = 2
vare = 100

## Theoretical power
p = 0.05
res = c()
for (i in 1:length(LDvec)) {
  LD = LDvec[i]
  R = matrix(c(1,LD,LD,1), nrow=2)
  eig = eigen(R, symmetric=TRUE)
  U = eig$vectors
  D = diag(eig$values)
  RinvSqrt = U %*% diag(sqrt(1/diag(D))) %*% t(U)
  for (j in 1:length(beta2vec)) {
    beta2 = beta2vec[j]
    varg = beta1^2 + beta2^2 + 2*LD*beta1*beta2
    vary = varg + vare
    q2 = varg/vary
    ncp = n*q2/(1-q2)
    cutoff = qchisq(p, m, ncp = 0, lower.tail = FALSE)
    power_mbat = pchisq(cutoff, m, ncp = ncp, lower.tail = FALSE)

    b1 = beta1 + LD*beta2
    b2 = beta2 + LD*beta1
    q2SNP1 = b1^2/vary
    q2SNP2 = b2^2/vary
    ncp1 = n*q2SNP1/(1-q2SNP1)
    ncp2 = n*q2SNP2/(1-q2SNP2)
    mu1 = sqrt(n/vary/(1-q2SNP1))*b1
    mu2 = sqrt(n/vary/(1-q2SNP2))*b2
    mu = c(mu1, mu2)
    ncp_transf = (t(U) %*% RinvSqrt %*% mu)^2
    cutoff = qgamma(p, shape = 1/(1+LD^2), scale = 2*(1+LD^2), lower.tail = FALSE)
    power_fastbat = psum.chisq(cutoff, lb = eig$values, df = rep(1,m), nc = ncp_transf, lower.tail = FALSE)

    res = rbind(res, data.frame(LD=LD, beta2=beta2, Method=c("mbat","fastbat"), Power=c(power_mbat,power_fastbat)))
  }
}

res.theo = res


## Empirical power
nrep = 1000
res = c()
for (i in 1:length(LDvec)) {
  LD = LDvec[i]
  R = matrix(c(1,LD,LD,1), nrow=2)
  eig = eigen(R, symmetric=TRUE)
  for (j in 1:length(beta2vec)) {
    beta2 = beta2vec[j]
    for (k in 1:nrep) {
      X = rmvnorm(n, mean=c(0,0), sigma=R)
      beta = c(beta1, beta2)
      g = X%*%beta
      y = g + rnorm(n, 0, 10)
      fit = apply(X, 2, function(x){summary(lm(y~x))$coefficients[2,1:2]})
      bhat = fit[1,]
      se = fit[2,]  
      z = bhat/se
      
      # fastbat
      P_fastbat = pchisqsum(sum(z^2), df=rep(1,m), a=eig$values, method="sad", lower=FALSE)
      
      # mbat
      UPz = crossprod(eig$vectors, z)
      chisq = crossprod(UPz, 1/eig$values * UPz) 
      P_mbat = pchisq(chisq, df=m, lower.tail=FALSE)
      
      # multiple regression
      #P_lm = anova(lm(y ~ X))$P[1]    
      
      res = rbind(res, data.frame(LD=LDvec[i], beta2=beta2, Method=c("fastbat", "mbat"), P=c(P_fastbat, P_mbat)))
    }
  }
}

res.emp = res
dd.emp = ddply(res.emp, .(LD, beta2, Method), summarise, Power = sum(P<0.05)/length(P))

res.theo$Scen = "Theoretical"
dd.emp$Scen = "Empirical"

res = rbind(res.theo, dd.emp)
res$LD = as.factor(res$LD)
res$beta2 = factor(res$beta2, levels=beta2vec, labels=sapply(beta2vec, function(x){paste("beta2 = ",x,sep="")}))
res$Scen = factor(res$Scen, levels=c("Theoretical","Empirical"))

#save(res, file="theoretical_power.RData")


## plot
p = ggplot(res, aes(LD, Power, col=Method, group=interaction(Method,Scen), shape=Scen, linetype=Scen)) + 
  geom_point() + geom_line() + 
  facet_wrap(~beta2, nrow=1) + ylim(c(0,1)) +
  theme(legend.title= element_blank())
#ggsave("theoretical_power.pdf", p, width=15, height=4)
