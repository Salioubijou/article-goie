library(stats4)
library(splines)
library(VGAM)
library("AdequacyModel")
library(moments)


# ************* Exponential********
# cdf of GIE-E distribution
cdf_GIE_E <- function(par,x){
  lambda <- par[1]
  gamma <- par[2]
  a <- par[3]
  G <-  pexp(x,a)
  cdf <- exp(-lambda*(1-G)*(1+gamma*G)/G)
  return(cdf)
}
# pdf of GIE-E distribution
pdf_GIE_E <- function(par,x){
  lambda <- par[1]
  gamma <- par[2]
  a <- par[3]
  g <- dexp(x,a)
  G <-  pexp(x,a)
  pdf <- lambda*g*(1+gamma*G^2)*exp(-lambda*(1-G)*(1+gamma*G)/G)/G^2
  return(pdf)
}
# log-Likelihood function of GIE_E distribution
lgie_E<-function(par,x){ 
  pdf<-rep(0,0);L <- length(x)
  for(i in 1:L)
    pdf[i] <- par[1]*dexp(x[i],par[3])*(1+par[2]*pexp(x[i],par[3])^2)*exp(-par[1]*(1-pexp(x[i],
                  par[3]))*(1+par[2]*pexp(x[i],par[3]))/pexp(x[i],par[3]))/pexp(x[i],par[3])^2
  return(sum(log(pdf)))
}
# *********************************************************
# cdf of GIE-W distribution
cdf_GIE_W <- function(par, x){
  lambda <- par[1]
  gamma <- par[2]
  a <- par[3]
  b <- par[4]
  G <-  pweibull(x,a,b)
  cdf <- exp(-lambda*(1-G)*(1+gamma*G)/G)
  return(cdf)
}
# pdf of GIE-W distribution
pdf_GIE_W <- function(par, x){
  lambda <- par[1]
  gamma <- par[2]
  a <- par[3]
  b <- par[4]
  g <- dweibull(x,a,b)
  G <-  pweibull(x,a,b)
  pdf <- (lambda*g*(1+gamma*G^2)*cdf_GIE_W(par, x))/G^2
  return(pdf)
}

# log-Likelihood function of GIE_W distribution
lgie_W<-function(par,x){ 
  pdf<-rep(0,0);L <- length(x)
  for(i in 1:L)
    pdf[i] <- par[1]*dweibull(x[i],par[3],par[4])*(1+par[2]*pweibull(x[i],par[3],par[4])^2
          )*exp(-par[1]*(1-pweibull(x[i],par[3],par[4]))*(1+par[2]*pweibull(x[i],par[3],
              par[4]))/pweibull(x[i],par[3],par[4]))/pweibull(x[i],par[3],par[4])^2
  return(sum(log(pdf)))
}
# *******************************************************************
# cdf of GIE-L distribution
cdf_GIE_L <- function(par,x){
  lambda <- par[1]; gamma <- par[2]
  theta <- par[3]
  G <-  plind(x,theta)
  cdf <- exp(-lambda*(1-G)*(1+gamma*G)/G)
  return(cdf)
}
# pdf of GIE-L distribution
pdf_GIE_L <- function(par, x){
  lambda <- par[1]; gamma <- par[2]
  theta <- par[3]
  g <- dlind(x,theta)
  G <-  plind(x,theta)
  pdf <- lambda*g*(1+gamma*G^2)*exp(-lambda*(1-G)*(1+gamma*G)/G)/G^2
  return(pdf)
}
# log-Likelihood function of GIE_L distribution
lgie_L<-function(par,x){ 
  pdf<-rep(0,0);L <- length(x)
  for(i in 1:L)
    pdf[i] <- par[1]*dlind(x[i],par[3])*(1+par[2]*plind(x[i],par[3])^2)*exp(-par[1]*(1-plind(x[i],
            par[3]))*(1+par[2]*plind(x[i],par[3]))/plind(x[i],par[3]))/plind(x[i],par[3])^2
  return(sum(log(pdf)))
}

# ******************** First data set *************************************

# OK Fatigue time of 101 6061-T6 aluminum coupons pdf 2931326
x2 <- c(70,90,96,97,99,100,103,104,104,105,107,108,108,108,
       109,109,112,112,113,114,114,114,116,119,120,120,120,
       121,121,123,124,124,124,124,124,128,128,129,129,130,
       130,131,131,131,131,131,132,132,132,133,134,134,134,
       134,134,136,136,137,138,138,138,139,139,141,141,142,
       142,142,142,142,142,144,144,145,146,148,148,149,151,
       151,152,155,156,157,157,157,157,158,159,162,163,163,
       164,166,166,168,170,174,196,212)

res2 <- goodness.fit(pdf=pdf_GIE_L, cdf=cdf_GIE_L,
                     starts = c(3,.1,1), data = x2, method = "N",
                     domain = c(0,Inf), mle = NULL)
res12 <- goodness.fit(pdf=pdf_GIE_W, cdf=cdf_GIE_W,
                      starts = c(.2,7,.11,1), data = x2, method = "N",
                      domain = c(0,Inf), mle = NULL)
res32 <- goodness.fit(pdf=pdf_GIE_E, cdf=cdf_GIE_E,
                      starts = c(2,3.2,1), data = x2, method = "N",
                      domain = c(0,Inf), mle = NULL)
        
summary(x2)
skewness(x2)
kurtosis(x2)

xp2 <- x2

# loglikelihood
lgie_E(res2$mle, x2)
lgie_W(res12$mle, x2)
lgie_L(res32$mle, x2)

# Test Likelihood Ratio (LR)
# Hypothises H0: gamma =0 versus H1: gamma ≠ 0
res121 <- goodness.fit(pdf=pdf_GIE_W, cdf=cdf_GIE_W,
                       starts = c(.2,0,.11,1), data = x2, method = "N",
                       domain = c(0,Inf), mle = NULL)
w2 <- 2*(lgie_W(res12$mle, x2)-lgie_W(res121$mle, x2))
pvalueLR2 =1 - pchisq(w2,1)
# 0.0002471858 < 0.001 then, we reject H0. 
# Fin test (LR)

# profile Loglikelihood when lambda is variable
# la represent lambda vector
la2 <- seq(.1,7,length.out = 50);la1 <- seq(.1,19,length.out = 50)
la3 <- seq(.1,6,length.out = 50);la4 <- seq(27,500,length.out = 50)
G2 <- G12 <- G32 <- G42 <- rep(0,0)
for (i in 1:50){
  G2[i] <- lgie_W(c(la2[i],res12$mle[-1]),x2)
  G12[i] <- lgie_W(c(res12$mle[1],la1[i],res12$mle[3:4]),x2)
  G32[i] <- lgie_W(c(res12$mle[-3:-4],la3[i],res12$mle[4]),x2)
  G42[i] <- lgie_W(c(res12$mle[-4],la4[i]),x2)
}

par(mfrow=c(1,1))
hist(x2,probability = T, ylim = c(0,.018), main = "Fitted pdf", ylab = "pdf",
     xlab = "x", col = "gray90")
box()
grid()
# ************* Estimated pdf *****************************
curve(pdf_GIE_L(x, par = res2$mle), col=1, add = T, lty=1,lwd=2)
curve(pdf_GIE_W(x, par = res12$mle), col=2, add = T, lty=2,lwd=2)
curve(pdf_GIE_E(x, par = res32$mle), col=4, add = T, lty=4,lwd=2)
legend("topleft",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1,cex = .8)

#************************ Estimated cdf ***********************

x2 <- seq(0,1,length.out=200)
plot(ecdf(xp2), lty=1, lwd=2, do.points=FALSE, verticals=TRUE, 
     ylim=c(0,1.0), ylab="cdf", main="Fitted cdf", col.01line="white",pch=4,cex=2)
grid()
curve(cdf_GIE_L(x, par = res2$mle), col=1, add = T, lty=1,lwd=2)
curve(cdf_GIE_W(x, par = res12$mle), col=2, add = T, lty=2,lwd=2)
curve(cdf_GIE_E(x, par = res32$mle), col=4, add = T, lty=4,lwd=2)
legend("topleft",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1, cex = .8)

par(mfrow=c(1,4))
plot(G2,type="l",col="red", ylab="profile Loglikelihood",xlab=expression(lambda),lwd=2)
grid()
plot(G12,type="l",col="#007FFF", ylab="profile Loglikelihood",xlab=expression(gamma),lwd=2)
grid()
plot(G32,type="l",col="#16B84E", ylab="profile Loglikelihood",xlab="a",lwd=2)
grid()
plot(G42,type="l",col="#2E006C", ylab="profile Loglikelihood",xlab="b",lwd=2)
grid()
# ******************* End first data set **************************

# **************** Second data set ********************************
x3 <- c(0.55, 0.93, 1.25, 1.36, 1.49, 1.52, 1.58, 1.61, 1.64, 1.68, 1.73, 1.81, 2, 0.74, 1.04, 1.27, 1.39,
        1.49, 1.53, 1.59, 1.61, 1.66, 1.68, 1.76, 1.82, 2.01, 0.77, 1.11, 1.28, 1.42, 1.5, 1.54, 1.6, 1.62,
        1.66, 1.69, 1.76, 1.84, 2.24, 0.81, 1.13, 1.29, 1.48, 1.5, 1.55, 1.61, 1.62, 1.66, 1.7, 1.77, 1.84,
        0.84, 1.24, 1.3, 1.48, 1.51, 1.55, 1.61, 1.63, 1.67, 1.7, 1.78, 1.89)

summary(x3)
skewness(x3)
kurtosis(x3)

xp3 <- x3

res3 <- goodness.fit(pdf=pdf_GIE_L, cdf=cdf_GIE_L,
                     starts = c(3,1,1), data = x3, method = "N",
                     domain = c(0,Inf), mle = NULL)
res13 <- goodness.fit(pdf=pdf_GIE_W, cdf=cdf_GIE_W,
                      starts = c(.2,7,.11,1), data = x3, method = "N",
                      domain = c(0,Inf), mle = NULL)
res33 <- goodness.fit(pdf=pdf_GIE_E, cdf=cdf_GIE_E,
                      starts = c(2,3.2,1), data = x3, method = "N",
                      domain = c(0,Inf), mle = NULL)

# loglikelihood
lgie_E(res3$mle, x3)
lgie_W(res13$mle, x3)
lgie_L(res33$mle, x3)

la32 <- seq(.001,.1,length.out = 50);la31 <- seq(75,178,length.out = 50)
la33 <- seq(2,7,length.out = 50); la34 <- seq(.2,15,length.out = 50); 
G32 <- G31 <- G33 <- G34 <- rep(0,0)
for (i in 1:50){
  G32[i] <- lgie_W(c(la32[i],res13$mle[-1]),x3)
  G31[i] <- lgie_W(c(res13$mle[1],la31[i],res13$mle[3:4]),x3)
  G33[i] <- lgie_W(c(res13$mle[c(-3,-4)],la33[i],res13$mle[4]),x3)
  G34[i] <- lgie_W(c(res13$mle[-4],la34[i]),x3)
}

par(mfrow=c(1,1))
hist(x3,probability = T, ylim = c(0,1.8), main = "Fitted pdf", ylab = "pdf",
     xlab = "x", col = "gray90")
box()
grid()
# ************* Estimated pdf *****************************
curve(pdf_GIE_L(x, par = res3$mle), col=1, add = T, lty=1,lwd=2)
curve(pdf_GIE_W(x, par = res13$mle), col=2, add = T, lty=2,lwd=2)
curve(pdf_GIE_E(x, par = res33$mle), col=4, add = T, lty=4,lwd=2)
legend("topleft",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1,cex = .8)

#************************ Estimated cdf ***********************

x3 <- seq(0,1,length.out=200)
plot(ecdf(xp3), lty=1, lwd=2, do.points=FALSE, verticals=TRUE, 
     ylim=c(0.0,1.0), ylab="cdf", main="Fitted cdf", col.01line="white",pch=4,cex=2)
grid()
curve(cdf_GIE_L(x, par = res3$mle), col=1, add = T, lty=1,lwd=2)
curve(cdf_GIE_W(x, par = res13$mle), col=2, add = T, lty=2,lwd=2)
curve(cdf_GIE_E(x, par = res33$mle), col=4, add = T, lty=4,lwd=2)
legend("topleft",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1,cex = .8)

par(mfrow=c(1,4))
plot(G32,type="l",col="red", ylab="profile Loglikelihood",xlab=expression(lambda),lwd=2)
grid()
plot(G31,type="l",col="#007FFF", ylab="profile Loglikelihood",xlab=expression(gamma),lwd=2)
grid()
plot(G33,type="l",col="#16B84E", ylab="profile Loglikelihood",xlab="a",lwd=2)
grid()
plot(G34,type="l",col="#2E006C", ylab="profile Loglikelihood",xlab="b",lwd=2)
grid()

# *************************** End seconde data set ******************

# ********************************* Third data set *************************
# x4 <- c(1.6, 2.0, 2.6, 3.0, 3.5, 3.9, 4.5, 4.6, 4.8, 5.0, 5.1, 5.3, 5.4, 5.6, 5.8, 6.0, 6.0, 6.1, 6.3, 6.5,
#         6.5, 6.7, 7.0, 7.1, 7.3, 7.3, 7.3, 7.7, 7.7, 7.8, 7.9, 8.0, 8.1, 8.3, 8.4, 8.4, 8.5, 8.7, 8.8, 9.0)

x4 <- c(1.901, 2.132, 2.203, 2.228, 2.257, 2.350, 2.361, 2.396, 2.397, 2.445, 2.454, 2.474, 2.518, 2.522,
        2.525, 2.532, 2.575, 2.614, 2.616, 2.618, 2.624, 2.659, 2.675, 2.738, 2.740, 2.856, 2.917, 2.928,
        2.937, 2.937, 2.977, 2.996, 3.030, 3.125, 3.139, 3.145, 3.220, 3.223, 3.235, 3.243, 3.264, 3.272, 
        3.294, 3.332, 3.346, 3.377, 3.408, 3.435, 3.493, 3.501, 3.537, 3.554, 3.562, 3.628, 3.852, 3.871, 
        3.886, 3.971, 4.024, 4.027, 4.225, 4.395, 5.020)
summary(x4)
skewness(x4)
kurtosis(x4)

xp4 <- x4
res4 <- goodness.fit(pdf=pdf_GIE_L, cdf=cdf_GIE_L,
                     starts = c(3,2,1), data = x4, method = "N",
                     domain = c(0,Inf), mle = NULL)
res14 <- goodness.fit(pdf=pdf_GIE_W, cdf=cdf_GIE_W,
                      starts = c(.2,7,.11,1), data = x4, method = "N",
                      domain = c(0,Inf), mle = NULL)
res34 <- goodness.fit(pdf=pdf_GIE_E, cdf=cdf_GIE_E,
                      starts = c(2,.2,1), data = x4, method = "N",
                      domain = c(0,Inf), mle = NULL)
# loglikelihood
lgie_E(res4$mle, x4)
lgie_W(res14$mle, x4)
lgie_L(res34$mle, x4)

la42 <- seq(3,5,length.out = 50);la41 <- seq(.01,70,length.out = 50)
la43 <- seq(1,1.8,length.out = 50);la44 <- seq(.1,15,length.out = 50);
G42 <- G41 <- G43 <- G44 <- rep(0,0)
for (i in 1:50){
  G42[i] <- lgie_W(c(la42[i],res14$mle[-1]),x4)
  G41[i] <- lgie_W(c(res14$mle[1],la41[i],res14$mle[3:4]),x4)
  G43[i] <- lgie_W(c(res14$mle[c(-3,-4)],la43[i],res14$mle[4]),x4)
  G44[i] <- lgie_W(c(res14$mle[-4],la44[i]),x4)
}

par(mfrow=c(1,1))
hist(x4,probability = T, ylim = c(0,.75), main = "Fitted pdf", ylab = "pdf",
     xlab = "x", col = "gray90")
box()
grid()
# ************* Estimated pdf *****************************
curve(pdf_GIE_L(x, par = res4$mle), col=1, add = T, lty=1,lwd=2)
curve(pdf_GIE_W(x, par = res14$mle), col=2, add = T, lty=2,lwd=2)
curve(pdf_GIE_E(x, par = res34$mle), col=4, add = T, lty=4,lwd=2)
legend("topright",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1,cex = .8)

#************************ Estimated cdf ***********************

x4 <- seq(0,1,length.out=200)
plot(ecdf(xp4), lty=1, lwd=2, do.points=FALSE, verticals=TRUE, 
     ylim=c(0.0,1.0), ylab="cdf", main="Fitted cdf", col.01line="white",pch=4,cex=2)
grid()
curve(cdf_GIE_L(x, par = res4$mle), col=1, add = T, lty=1,lwd=2)
curve(cdf_GIE_W(x, par = res14$mle), col=2, add = T, lty=2,lwd=2)
curve(cdf_GIE_E(x, par = res34$mle), col=4, add = T, lty=4,lwd=2)
legend("topleft",c("GIE-L","GIE-W","GIE-E"),lwd=c(2,2,2),
       lty=c(1,2,4),col=c(1,2,4), bg="gray90", box.lwd = .1,cex = .8)

par(mfrow=c(1,4))
plot(G42,type="l",col="red", ylab="profile Loglikelihood",xlab=expression(lambda),lwd=2)
grid()
plot(G41,type="l",col="#007FFF", ylab="profile Loglikelihood",xlab=expression(gamma),lwd=2)
grid()
plot(G43,type="l",col="#16B84E", ylab="profile Loglikelihood",xlab="a",lwd=2)
grid()
plot(G44,type="l",col="#2E006C", ylab="profile Loglikelihood",xlab="b",lwd=2)
grid()

# ********************************* End data set 3 ***************************

