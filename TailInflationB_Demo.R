source('TailInflationB.R')

### Example 1: Piecewise linear theta ###

# Generate data from P such that
# P(dx) = e^theta(x) P0(dx):
n       <- 400
alpha   <- 1
beta    <- 1
# Preliminary choice of tau and theta for rate beta = 1: 
t.true  <- c(0, 2.0, 4.0, 6.0)
th.true <- c(0, 0, 1, 2.2, 0.8)
# Normalization of theta:
th.true <- LocalNormalize.2B(t.true, th.true, alpha)
# Rescaling to actual beta:
t.true  <- t.true/beta
th.true <- c(th.true[1:length(t.true)],
             th.true[1+length(t.true)]*beta)

# Generate data
x <- Simulate.2B(n, t.true, th.true, alpha, beta)
w <- rep(1/n,n)

# Estimation of tail inflation function with
# graphical display of all intermediate steps
# (make sure the console window is in the foreground
#  and graphics window beside but visible):
res <- TailInflation.B0(x,w,alpha,beta)

# Same computation without graphics:
res <- TailInflation.B(x,w,alpha,beta)
system.time(res <- TailInflation.B(x,w,alpha,beta))

t.hat  <- res$tau
th.hat <- res$theta
tt     <- seq(0,max(x),length.out=501)
tt     <- unique(sort(c(tt,t.hat)))

# Depicting directional derivatives to check
# optimality:
x0      <- x*beta
t.hat0  <- t.hat*beta
th.hat0 <- c(th.hat[1:length(t.hat)],
             th.hat[1+length(t.hat)]/beta)
tt0     <- tt*beta
H.tt0 <- LocalDirDeriv2.2B(tt0,x0,w,t.hat0,th.hat0,alpha)
# par(cex=1.5,mai=c(1.0,1.1,0.1,0.05),mgp=c(2.3,1,0))
plot(tt0,H.tt0[,1],type='l',lwd=2,
     xlab=expression(italic(t)),
     ylab=expression(italic(h(t))))
abline(h=0,col='red')
abline(v=t.hat0,col='blue')
rug(x0)
# The function h should be <= 0 with equality at t.hat...


# Display true (green) and estimated (black)
# tail inflation function:
th.true.tt <- Evaluate.2B(tt,t.true,th.true)
D.tt <- LocalLinearSplines1.2B(t.hat,tt)
th.hat.tt <- D.tt %*% th.hat
plot(tt,th.true.tt,type='l',lwd=2,col='green',
	xlab=expression(italic(x)),
	ylab=expression(italic(theta(x))),
	ylim=range(c(th.true.tt,th.hat.tt)))
abline(v=t.hat,col='blue')
lines(tt,th.true.tt,lty=5)
lines(tt,th.hat.tt,lwd=2)
rug(x)

# Display true (green), estimated (black) and
# reference density (magenta):
p0.tt     <- dgamma(tt,alpha,beta)
p.true.tt <- exp(th.true.tt)*p0.tt
p.hat.tt  <- exp(th.hat.tt)*p0.tt
plot(tt,p0.tt,type='l',lwd=2,col='magenta',
     xlab=expression(italic(x)),
     ylab=expression(italic(p(x))),
     ylim=c(0,max(p0.tt[2],p.true.tt[2],p.hat.tt[2])))
lines(tt,p.true.tt,col='green',lwd=2)
lines(tt,p.hat.tt,lwd=2)
rug(x)


### Example 2: Scale mixtures of chi_1^2 r.v.s ###

n <- 1000
w <- rep(1/n,n)
alpha <- 1/2
beta  <- 1/2

# Number of observations with alternative scale:
k <- 100
# Alternative scale:
S1 <- 2
S2 <- 1.4

# set.seed(71220)
x <- rchisq(n,df=1)
x[1:k] <- S1*x[1:k]
x[(k+1):(2*k)] <- S2*x[(k+1):(2*k)]
x <- sort(x)

# Estimation of tail inflation function with
# graphical display of all intermediate steps
# (make sure the console window is in the foreground
#  and graphics window beside but visible):

res <- TailInflation.B0(x,w,alpha,beta)

# Same computation without graphics:
res <- TailInflation.B(x,w,alpha,beta)
system.time(res <- TailInflation.B(x,w,alpha,beta))

t.hat <- res$tau
th.hat <- res$theta
tt <- seq(0,max(x),length.out=501)
tt <- unique(sort(c(tt,t.hat)))

TLR <- n*res$LL
TLR

# Depicting directional derivatives to check
# optimality:
x0      <- x*beta
t.hat0  <- t.hat*beta
th.hat0 <- c(th.hat[1:length(t.hat)], 
             th.hat[1+length(t.hat)]/beta)
tt0     <- tt*beta
H.tt0 <- LocalDirDeriv2.2B(tt0,x0,w,t.hat0,th.hat0,alpha)
# par(cex=1.5,mai=c(1.0,1.1,0.1,0.05),mgp=c(2.3,1,0))
plot(tt,H.tt0[,1],type='l',lwd=2,
     xlab=expression(italic(t)),
     ylab=expression(italic(h(t))))
abline(h=0,col='red')
abline(v=t.hat,col='blue')
rug(x)
# The function h should be <= 0 with equality at t.hat...

# Display true (green) and estimated (black)
# tail inflation function:
th.true.tt <- log(n - 2*k +
		k*sqrt(1/S1)*exp(0.5*(1 - 1/S1)*tt) +
		k*sqrt(1/S2)*exp(0.5*(1 - 1/S2)*tt)) - log(n)
D.tt <- LocalLinearSplines1.2B(t.hat,tt)
th.hat.tt <- D.tt %*% th.hat
plot(tt,th.true.tt,type='l',lwd=2,col='green',
     xlab=expression(italic(x)),
     ylab=expression(italic(theta(x))),
     ylim=range(c(th.true.tt,th.hat.tt)))
abline(v=t.hat,col='blue')
lines(tt,th.hat.tt,lwd=2)
lines(tt,th.true.tt,lty=5)
rug(x)

# Display true (green), estimated (black) and
# reference density (magenta):
p0.tt     <- dgamma(tt,alpha,rate=beta)
p.true.tt <- ((n-2*k)*dgamma(tt,alpha,rate=beta) +
              k*dgamma(tt,alpha,rate=beta/S1) +
              k*dgamma(tt,alpha,rate=beta/S2))/n
p.hat.tt  <- exp(th.hat.tt)*p0.tt
plot(tt,p0.tt,type='l',lwd=2,col='magenta',
     xlab=expression(italic(x)),
     ylab=expression(italic(p(x))),
     ylim=c(0,max(p0.tt[2],p.true.tt[2],p.hat.tt[2])))
lines(tt,p.true.tt,col='green',lwd=2)
lines(tt,p.true.tt,lty=5)
lines(tt,p.hat.tt,lwd=2)
rug(x)
