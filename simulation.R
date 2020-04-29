library(AdequacyModel)
library(stats4)
library(splines)
library(VGAM)

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

qgie <- function(u, lambda, gamma, a, b){
  if(gamma==0)
    return(qweibull(lambda/(lambda-log(u)), a, b))
  return(qweibull((log(u)+lambda*(gamma-1)+((log(u)+
          lambda*(gamma-1))**2+4*gamma*lambda**2)**0.5)/(2*lambda*gamma),a, b))
}

lambda <- gam <- a <- b <- NULL
slambda <- sgam <- sa <- sb <- NULL
llambda <- llgam <- la <- lb <- ulambda <- ugam <- ua <- ub <-  NULL
ALlambda <- ALgam <- ALa <-ALb <-  NULL
cplambda <- cpgam <- cpa<- cpb <- NULL
biaslambda <- biasgam <- biasa <-biasb <-  NULL
mselambda <- msegam <- msea <-mseb <-  NULL

N <- 350
a1 <- 2; gam1 <- 1; theta1 <- 1; b1 <- .5
set.seed(15000)
for( k in 1:221) {
  for( j in 1:N) {
    echan=NULL
    if(k <= 50)
      echan <- qgie(runif(155,0,1),a1,gam1,theta1,b1)
    if(k>50 && k <=100)
      echan <- qgie(runif(655,0,1),a1,gam1,theta1,b1)
    if(k>100 && k <=150)
      echan <- qgie(runif(870,0,1),a1,gam1,theta1,b1)
    else 
      echan <- qgie(runif(1000,0,1),a1,gam1,theta1,b1)
    res_gf=goodness.fit(pdf=pdf_GIE_W, cdf=cdf_GIE_W, starts = c(a1,gam1,theta1,b1), 
              data = echan,method="N", domain=c(0,Inf), mle = NULL)
  
    lambda[j]=res_gf$mle[1]
    gam[j]=res_gf$mle[2] 
    a[j]=res_gf$mle[3]
    b[j]=res_gf$mle[4]
    
    slambda[j]=res_gf$Erro[1]
    sgam[j]=res_gf$Erro[2]
    sa[j]=res_gf$Erro[3]
    sb[j]=res_gf$Erro[4]

    llambda[j]=lambda[j]-1.96*slambda[j]
    ulambda[j]=lambda[j]+1.96*slambda[j]
    llgam[j]=gam[j]-1.96*sgam[j]
    ugam[j]=gam[j]+1.96*sgam[j]
    la[j]=a[j]-1.96*sa[j]
    ua[j]=a[j]+1.96*sa[j]
    lb[j]=b[j]-1.96*sb[j]
    ub[j]=b[j]+1.96*sb[j]
}
  
  ALlambda[k]=3.92*mean(slambda)
  ALgam[k]=3.92*mean(sgam)
  ALa[k]=3.92*mean(sa)
  ALb[k]=3.92*mean(sb)

  cplambda[k]=mean((llambda) < a1 & (ulambda) > a1)
  cpgam[k]=mean((llgam) < gam1 & (ugam) > gam1)
  cpa[k]=mean((la) < theta1 & (ua) >theta1)
  cpa[k]=mean((lb) < b1 & (ub) >theta1)
  
  biaslambda[k]=mean(lambda-a1) 
  biasgam[k]=mean(gam-gam1) 
  biasa[k]=mean(a-theta1)
  biasb[k]=mean(b-b1)
  
  mselambda[k]=mean((lambda-a1)^2)
  msegam[k]=mean((gam-gam1)^2)
  msea[k]=mean((a-theta1)^2)
  mseb[k]=mean((b-b1)^2)
}


# bias plot
plot(ALa, type = "l", col="red", main = expression(paste(theta)), xlab = "n",
     xlim=c(0,130), ylab="AL")
grid()
plot(ALlambda, type = "l", col="blue", main = expression(paste(lambda)), xlab = "n",
     xlim=c(0,130), ylab="AL")
grid()
plot(ALgam, type = "l", col="green", main = expression(paste(gamma)), xlab = "n",
     xlim=c(0,130), ylab="AL")
grid()

par(mfrow=c(1,1))
plot(cpa, type = "l", xlab = "n", xlim=c(0,130), col="blue",ylim=c(.8,1), ylab="cp")
points(cplambda, type = "l", col=2)
points(cpgam, type = "l", col=3)
abline(0.95,0, col="#EE1010")
legend("bottomright",c(expression(lambda),expression(gamma),expression(theta)),
       lwd=c(4,4,4),col=c(4,2,3), bg="gray90", box.lwd = .1)
grid()

par(mfrow=c(1,3))
plot(biasa, type = "l", col="red", main = expression(paste(theta)), xlab = "n",
     xlim=c(0,130), ylab="Estimated Bias")
grid()
plot(biaslambda, type = "l", col="blue", main = expression(paste(lambda)),
     xlab = "n", xlim=c(0,130), ylab="Estimated Bias")
grid()
plot(biasgam, type = "l", col="green", main = expression(paste(gamma)),
     xlab = "n", xlim=c(0,130), ylab="Estimated Bias")
grid()

# mse plot
plot(msea, type = "l", col="red", main = expression(paste(theta)),
     xlab = "n", xlim=c(0,130), ylab="Estimated MSE")
grid()
plot(mselambda, type = "l", col="blue", main = expression(paste(lambda)),
     xlab = "n", xlim=c(0,130), ylab="Estimated MSE")
grid()
plot(msegam, type = "l", col="green", main = expression(paste(gamma)),
     xlab = "n", xlim=c(0,130), ylab="Estimated MSE")
grid()
#

