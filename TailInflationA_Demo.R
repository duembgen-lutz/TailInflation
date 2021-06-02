source('TailInflationA.R')

# Simulate heterogeneous Gaussian data:
n <- 400
w <- rep(1/n,n)
mu <- 1.5
sigma <- 1.1
k <- floor(sqrt(n))
x <- rnorm(n)
x[1:k] <- sigma*x[1:k] + mu
x <- sort(x)


# Estimation of tail inflation function with
# graphical display of all intermediate steps
# (make sure the console window is in the foreground
#  and graphics window beside but visible):

res <- TailInflation.A0(x)

# Same computation without graphics:
res <- TailInflation.A(x)
system.time(res <- TailInflation.A(x))


# Likelihood ratio test statistic for the null hypothesis
# that all observations are standard Gaussian:
TLR <- n*res$LL
TLR

tau <- res$tau
theta <- res$theta

tt <- seq(min(-3.5,mu-3.5*sigma),max(3.5,mu+3.5*sigma),length.out=1001)
tt <- unique(sort(c(x,tau,tt)))

# Check optimality of fit by plotting directional
# derivatives for extra kinks:
# par(cex=1.5,mai=c(1.0,1.1,0.1,0.05),mgp=c(2.3,1,0))
H.tt <- LocalDirDeriv2.2A(tt,x,w,tau,theta)
plot(tt,H.tt[,1],type='l',lwd=2,
	xlab=expression(italic(t)),
	ylab=expression(italic(h(t))))
abline(h=0,col='red')
abline(v=tau,col='blue')
rug(x)

# Display true (green) and estimated (black)
# tail inflation function:
th0.tt <- log(1 - k/n + (k/n)*exp(mu*tt - mu^2/2))
D.tt <- LocalLinearSplines1.2A(tau,tt)
thh.tt <- D.tt %*% theta
plot(tt,th0.tt,type='l',lwd=2,col='green',
	xlab=expression(italic(x)),
	ylab=expression(italic(theta(x))),
	ylim=c(-0.3,4)
	#ylim=range(c(th0.tt,thh.tt))
	)
abline(v=tau,col='blue')
lines(tt,th0.tt,lty=5)
lines(tt,thh.tt,lwd=2)
rug(x)

# Display true (green), estimated (black) and
# reference density (magenta):
p0.tt <- dnorm(tt)
f0.tt <- exp(th0.tt)*p0.tt
fh.tt <- exp(thh.tt)*p0.tt
plot(tt,p0.tt,type='l',lwd=2,col='magenta',
	xlab=expression(italic(x)),
	ylab=expression(italic(p(x))),
	ylim=c(0,max(c(p0.tt,f0.tt,fh.tt))))
abline(h=0)
lines(tt,f0.tt,col='green',lwd=3)
lines(tt,f0.tt,lty=5)
lines(tt,fh.tt,lwd=2)
rug(x)


# Likelihood ratio test statistic for the null hypothesis
# that all observations are standard Gaussian:
TLR <- n*res$LL
TLR

# Monte-Carlo p-value for the null hypotheses that
# the observations are i.i.d. standard Gaussian:

mcsim <- 999
Tv <- rep(0,mcsim+1)
Tv[1] <- TLR
for (s in 1:mcsim){
	Tv[s+1] <- n*TailInflation.A(sort(rnorm(n)))$LL
}
pv <- mean(Tv >= Tv[1])
pv
