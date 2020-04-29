library(stats4)
library(splines)
library(VGAM)
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
# hrf of GIE-E distribution
hrf_GIE_E <- function(par,x){
  lambda <- par[1]
  gamma <- par[2]
  a <- par[3]
  g <- dexp(x,a)
  G <-  pexp(x,a)
  hrf <- (lambda*g*(1+gamma*G^2)*exp(-lambda*(1-G)*(1+gamma*G)/G)/
            G^2)/(1-exp(-lambda*(1-G)*(1+gamma*G)/G))
  return(hrf)
}

curve(pdf_GIE_E(c(3,-1,1),x), ylab = "density-E", xlim = c(0,5), ylim=c(0,1.1),lwd=2)
curve(pdf_GIE_E(c(5,-.1,1),x), add = T, col=2, lty=2,lwd=2)
curve(pdf_GIE_E(c(1,.3,1),x), add = T, col=3, lty=3,lwd=2)
curve(pdf_GIE_E(c(.5,1, 1),x), add = T, col=4, lty=4,lwd=2)
curve(pdf_GIE_E(c(.4,1.5,.7),x), add = T, col=5, lty=5,lwd=2)
legend("topright", legend = c(expression(lambda==3~~gamma==-1~~theta==1),
                              expression(lambda==5~~gamma==-0.1~~theta==1),
                              expression(lambda==1~~gamma==0.3~~theta==1),
                              expression(lambda==0.5~~gamma==1~~theta==1),
                              expression(lambda==0.4~~gamma==1.5~~theta==0.7)),
       lwd=c(2,2,2,2,2) ,lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()

curve(hrf_GIE_E(c(1, .2,1),x), ylab = "hazard rate-E", xlim = c(0,4.2), ylim=c(0,1.1),lwd=2)
curve(hrf_GIE_E(c(1.5,.8,1),x), add = T, col=2, lty=2,lwd=2)
curve(hrf_GIE_E(c(5,-.1,1),x), add = T, col=3, lty=3,lwd=2)
curve(hrf_GIE_E(c(.4,2,1),x), add = T, col=4, lty=4,lwd=2)
curve(hrf_GIE_E(c(1,2,1),x), add = T, col=5, lty=5,lwd=2)
legend("bottomright", legend = c(expression(lambda==1~~gamma==0.2~~theta==.5),
                              expression(lambda==1.5~~gamma==0.8~~theta==1),
                              expression(lambda==5~~gamma==-0.1~~theta==1), 
                              expression(lambda==0.4~~gamma==2~~theta==1),
                              expression(lambda==1~~gamma==2~~theta==1)),
       lwd=c(2,2,2,2,2) ,lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()
# ************************* Fin exponential ******************

# **************************** Weibull *****************************
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

# hrf of GIE-W distribution
hrf_GIE_W <- function(par, x){
  pdf_GIE_W(par, x)/(1-cdf_GIE_W(par, x))
}
curve(pdf_GIE_W(c(1,4,2,1.5),x), ylab = "density-W", xlim = c(0,7), ylim=c(0,1.1),lwd=2)
curve(pdf_GIE_W(c(2,4,.8,.7),x), add = T, col=2, lty=2,lwd=2)
curve(pdf_GIE_W(c(1.5,3,2,1),x), add = T, col=3, lty=3,lwd=2)
curve(pdf_GIE_W(c(7,1,1.4,1.2),x), add = T, col=4, lty=4,lwd=2)
curve(pdf_GIE_W(c(1,3,1,1),x), add = T, col=5, lty=5,lwd=2)
legend("topright", legend = c(expression(lambda==1~~gamma==4~~a==2~~b==1.5),
                              expression(lambda==2~~gamma==4~~a==0.8~~b==0.7),
                              expression(lambda==1.5~~gamma==3~~a==0.5~~b==1),
                              expression(lambda==7~~gamma==1~~a==1.4~~b==1.2),
                              expression(lambda==1~~gamma==3~~a==1~~b==1)),
       lwd=c(2,2,2,2,2),lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()
curve(hrf_GIE_W(c(1,-.5,2,1.5),x), ylab = "hazard rate-W", ylim = c(0,4.5), xlim = c(0,3),lwd=2)
curve(hrf_GIE_W(c(2,0,.3,.1),x), add = T, col=2, lty=2,lwd=2)
curve(hrf_GIE_W(c(4,-.3,1,.6),x), add = T, col=3, lty=3,lwd=2)
curve(hrf_GIE_W(c(7,-1,1,.7),x), add = T, col=4, lty=4,lwd=2)
curve(hrf_GIE_W(c(.5,.03,2,2),x), add = T, col=5, lty=5,lwd=2)
legend("topright", legend = c(expression(lambda==1~~gamma==-.5~~a==02~~b==1.5),
                              expression(lambda==2~~gamma==-0~~a==.3~~b==.1),
                              expression(lambda==4~~gamma==-0.3~~a==1~~b==0.6),
                              expression(lambda==7~~gamma==-1~~a==1~~b==0.7),
                              expression(lambda==.5~~gamma==-0.03~~a==2~~b==2)),
       lwd=c(2,2,2,2,2),lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()
# ************************ Fin Weibull **************************
# **************************** Lindley *****************************
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
# hrf of GIE-L distribution
hrf_GIE_L <- function(par,x){
  lambda <- par[1]; gamma <- par[2]
  theta <- par[3]
  g <- (theta^2)*(1+x)*exp(-theta*x)/(1+theta)
  G <-  1-((1+theta+theta*x)*exp(-theta*x))/(1+theta)
  hrf <- (lambda*g*(1+gamma*G^2)*exp(-lambda*(1-G)*(1+gamma*G)/G)/G^2)/
    (1-exp(-lambda*(1-G)*(1+gamma*G)/G))
  return(hrf)
}
# Traçons la courbe de la pdf_GIE_L
curve(pdf_GIE_L(c(2,.5,2),x), ylab = "density-L", xlim = c(0,7), ylim=c(0,1.3),lwd=2)
curve(pdf_GIE_L(c(7,3,4),x), add = T, col=2, lty=2,lwd=2)
curve(pdf_GIE_L(c(5,-1,.7),x), add = T, col=3, lty=3,lwd=2)
curve(pdf_GIE_L(c(2.7,.3, 1),x), add = T, col=4, lty=4,lwd=2)
curve(pdf_GIE_L(c(3,-1, 1),x), add = T, col=5, lty=5,lwd=2)
legend("topright", legend = c(expression(lambda==2~~gamma==0.5~~theta==2),
                              expression(lambda==7~~gamma==3~~theta==4),
                              expression(lambda==5~~gamma==-1~~theta==0.7),
                              expression(lambda==2.7~~gamma==0.3~~theta==1),
                              expression(lambda==3~~gamma==-1~~theta==1)),
       lwd=c(2,2,2,2,2),lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()
# Traçons la courbe de la hrf_GIE_L
curve(hrf_GIE_L(c(11,2,4.5),x), ylab = "hazard rate-L", xlim = c(0,13), ylim=c(0,7.5),lwd=2)
curve(hrf_GIE_L(c(12,4,3),x), add = T, col=2, lty=2,lwd=2)
curve(hrf_GIE_L(c(15,3,4),x), add = T, col=3, lty=3,lwd=2)
curve(hrf_GIE_L(c(1,5,5),x), add = T, col=4, lty=4,lwd=2)
curve(hrf_GIE_L(c(10,5,3.4),x), add = T, col=5, lty=5,lwd=2)
legend("topright", legend = c(expression(lambda==11~~gamma==2~~theta==4.5),
                              expression(lambda==12~~gamma==4~~theta==3),
                              expression(lambda==15~~gamma==3~~theta==4),
                              expression(lambda==1~~gamma==5~~theta==5),
                              expression(lambda==10~~gamma==5~~theta==3.4)),
       lwd=c(2,2,2,2,2),lty=1:5, col=1:5, cex = .8,bg="gray90", box.lwd = .1)
grid()

# ************************ Fin Lindley **************************

# plot of skewness and kurtosis of GIE-weibull distribution
# qgie quantile function of GIE-G
qgie <- function(u, lambda, gamma, a, b){
  if(gamma==0)
    return(qweibull(lambda/(lambda-log(u)), a, b))
  return(qweibull((log(u)+lambda*(gamma-1)+((log(u)+
      lambda*(gamma-1))**2+4*gamma*lambda**2)**0.5)/(2*lambda*gamma),a, b))
}

vlambda <- seq(1,3,length.out = 5)
vgamma <- seq(-1,2,length.out = 5)
# Skewness function
skew1 <- function(a, b){
  res <- matrix(data = NA, ncol = 5, nrow = 5)
  for (j in 1:5) {
    for (i in 1:5) {
      lambda <- vlambda[i]
      gamma <- vgamma[j]
      res[i,j] <- (qgie(.25, lambda, gamma, a, b)+qgie(.75, lambda, gamma, a, b)
                 -2*qgie(.5, lambda, gamma, a, b))/
        (qgie(.75, lambda, gamma, a, b)-qgie(.25, lambda, gamma, a, b))
    }
  }
  return(res)
}
skewness1 <- skew1(1,1)
library(plotly)
plot_ly(x = vlambda, y = vgamma, z = ~skewness1, type = "surface") %>%layout(
  title = "Plot of GOIE-W skewness",
  scene = list(xaxis = list(title = paste(expression(lambda))),
    yaxis = list(title = paste(expression(gamma))),
    zaxis = list(title = "skewness")))

# Kurtosis function
kurto <- function(a, b){
  res <- matrix(data = NA, ncol = 5, nrow = 5)
  for (j in 1:5) {
    for (i in 1:5) {
      lambda <- vlambda[i]
      gamma <- vgamma[j]
      res[i,j] <- (qgie(7/8, lambda, gamma, a, b)-qgie(5/8, lambda, gamma, a, b)+
                   qgie(3/8, lambda, gamma, a, b)-qgie(1/8, lambda, gamma, a, b))/
        (qgie(6/8, lambda, gamma, a, b)-qgie(2/8, lambda, gamma, a, b))
    }
  }
  return(res)
}
kurtosis1 <- kurto(1,1)
plot_ly(x = vlambda, y = vgamma, z = ~kurtosis1, type = "surface") %>%layout(
  title = "Plot of GOIE-W kurtosis",
  scene = list(xaxis = list(title = paste(expression(lambda))),
    yaxis = list(title = paste(expression(gamma))),
    zaxis = list(title = "kurtosis")))

