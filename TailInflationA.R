#######################################################
### Active Set Algorithm for the Computation of the ###
### NPMLE of a Log-Convex Density Ratio             ###
#######################################################
#
# Lutz Duembgen, Alexandre Moesching, Christof Straehl
# December 2020
#
# Source:
# L. Duembgen, A. Moesching, C. Straehl (2021).
# Active Set Algorithms for Estimating Shape-Constrained
# Density Ratios.
# Computational Statistics and Data Analysis 163.
# arXiv:1808.09340


### Main programs ----

TailInflation.A0 <- function(X,W=rep(1,length(X)),
	xmin=1.1*min(X)-0.1*max(X),
	xmax=1.1*max(X)-0.1*min(X),
	delta1=10^(-10)/length(X),
	delta2=10^(-3)/length(X))
# Computation of a log-convex density ratio with respect to
# the standard Gaussian distribution.
# Graphical display of all intermediate steps...
{
	# Preparations of raw data:
	n <- length(X)
	if (sum(X[1:(n-1)] >= X[2:n]) == 0){
		x <- X
		w <- W/sum(W)
	}else{
		tmp <- LocalPrepareData(X,W)
		x <- tmp$x
		w <- tmp$w
		n <- tmp$n
	}
	# Gaussian MLE theta (for variance 1):
	mu <- sum(w*x)
	LL <- sum(w*(mu*x - mu^2/2))
	tmp <- Optimality1.2A(x,w,mu)
	if (length(tmp$tnew) > 0){
		ht0 <- max(tmp$htnew)
		t0 <- tmp$tnew[which.max(tmp$htnew)]
	}else{
		ht0 <- 0
	}
		
	message <- readline("First plot:")
	plot(c(xmin,xmax),rep(0,2),type='l',
		lwd=3,col='green',
		xlab=expression(italic(x)),
		ylab=expression(italic(theta(x))),
		xlim=range(x),
		ylim=range(mu*c(xmin,xmax) - mu^2/2),
		main=paste('Starting point: ',
			'LL = ',round(LL,digits=9),sep=''))
	lines(c(xmin,xmax),mu*c(xmin,xmax)-mu^2/2,lwd=2)
	rug(x)
	
	if (ht0 < delta2){
		print(paste('Parametric fit sufficient: N(',
			round(mu,4),', 1)',sep=''))
		return(list(tau=NULL,theta=mu,x=x,w=w,m=0,
			LL=LL,nr.local.search=0,nr.Newton=0))
	}

	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0
	
	tau <- NULL
	theta <- mu
	
	while (max(ht0) > delta2){
		nr.local.search <- nr.local.search + 1
		# Store current value of LL:
		LL.old <- LL
		# Add new knots
		tmp <- LocalNewKnots.2A(tau,theta,t0,ht0,delta2)
		tau <- tmp$knots
		theta <- tmp$theta
		isnewknot <- tmp$isnewknot
		# Graphical display:
		message <- readline("New knots:")
		xtx <- c(xmin,tau,xmax)
		Dtmp <- LocalLinearSplines1.2A(tau,xtx)
		thetaxtx <- Dtmp %*% theta
		plot(c(xmin,xmax),rep(0,2),type='l',
			lwd=3,col='green',
			xlab=expression(italic(x)),
			ylab=expression(italic(theta(x))),
			xlim=range(x),
			ylim=range(thetaxtx),
			main=paste('LL =',round(LL,digits=9)))
		lines(xtx,thetaxtx,lwd=2)
		points(tau,thetaxtx[2:(length(tau)+1)])
		abline(v=tau[isnewknot],col='blue')
		rug(x)
		
		wt <- LocalLinearSplines2.2A(tau,x,w)
		proposal <- LocalNewton.2A(wt,tau,theta,delta1)
		nr.Newton <- nr.Newton + 1
		while (proposal$dirderiv > delta1){
			theta.new <- proposal$theta.new
			# Graphical display:
			message <- readline("New proposal:")
			thetaxtx.new <- Dtmp %*% theta.new
			plot(c(xmin,xmax),rep(0,2),type='l',
				lwd=3,col='green',
				xlab=expression(italic(x)),
				ylab=expression(italic(theta(x))),
				xlim=range(x),
				ylim=range(thetaxtx),
				main=paste('LL =',
					round(LL,digits=9)))
			lines(xtx,thetaxtx,lwd=2)
			lines(xtx,thetaxtx.new,lwd=2,col='blue')
			points(tau,thetaxtx.new[2:(length(tau)+1)],
				col='blue')
			abline(v=tau[isnewknot],col='blue')
			rug(x)
			# 2nd step size correction:
			corr <- LocalStepSize2.2A(tau,theta,theta.new,isnewknot)
			theta <- corr$theta
			if (length(corr$tau) < length(tau)){
				LL.ref <- -Inf
				# Update of knots, wt and isnewknot:
				tau <- corr$tau
				isnewknot <- corr$isnewknot
				wt <- LocalLinearSplines2.2A(tau,x,w)
				# Graphical display:
				message <- readline("2nd step size correction:")
				xtx <- c(xmin,tau,xmax)
				Dtmp <- LocalLinearSplines1.2A(tau,xtx)
				thetaxtx <- Dtmp %*% theta
				lines(xtx,thetaxtx,lwd=2,col='magenta')
			}else{
				LL.ref <- LL
			}
			# Normalization:
			theta <- LocalNormalize.2A(tau,theta)
			LL <- sum(wt*theta)
			if (LL > LL.ref){
				proposal <- LocalNewton.2A(wt,tau,theta,delta1)
				nr.Newton <- nr.Newton + 1
				# For graphical display:
				thetaxtx <- Dtmp %*% theta
			}else{
				proposal$dirderiv <- 0	
			}
		}

		# Graphical display:
		message <- readline("Local optimum:")
		thetaxtx     <- Dtmp %*% theta
		plot(c(xmin,xmax),rep(0,2),type='l',
			lwd=3,col='green',
			xlab=expression(italic(x)),
			ylab=expression(italic(theta(x))),
			xlim=range(x),
			ylim=range(thetaxtx),
			main=paste('LL =',
				round(LL,digits=9)))
		lines(xtx,thetaxtx,lwd=2)
		points(tau,thetaxtx[2:(length(tau)+1)])
		rug(x)
		
		# Check for global optimality:
		if (LL > LL.old){
			tmp <- Optimality2.2A(x,w,tau,theta)
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
			kk <- ht0 >= max(delta2, max(ht0)*0.0001)
			t0 <- t0[kk]
			if (length(t0) > 0){
				ht0 <- ht0[kk]
			} else {
				ht0 <- 0
			}
		}else{
			ht0 <- 0
		}
	}
	
	# Graphical display:
	message <- readline("Final estimate:")
	thetaxtx <- Dtmp %*% theta
	plot(c(xmin,xmax),rep(0,2),type='l',
		lwd=3,col='green',
		xlab=expression(italic(x)),
		ylab=expression(italic(theta(x))),
		xlim=range(x),
		ylim=range(thetaxtx),
		main=paste('LL =',
			round(LL,digits=9)))
	lines(xtx,thetaxtx,lwd=2)
	points(tau,thetaxtx[2:(length(tau)+1)])
	rug(x)
	
	return(list(tau=tau,theta=theta,x=x,w=w,m=length(tau),
		LL=LL,
		nr.local.search=nr.local.search,nr.Newton=nr.Newton))
}

TailInflation.A <- function(X,W=rep(1,length(X)),
	delta1=10^(-10)/length(X),
	delta2=10^(-3)/length(X))
# Computation of a log-convex density ratio with respect to
# the standard Gaussian distribution.
# Pure computation, no graphical displays...
{
	# Preparations of raw data:
	n <- length(X)
	if (sum(X[1:(n-1)] >= X[2:n]) == 0){
		x <- X
		w <- W/sum(W)
	}else{
		tmp <- LocalPrepareData(X,W)
		x <- tmp$x
		w <- tmp$w
		n <- tmp$n
	}
	# Gaussian MLE theta (for variance 1):
	mu <- sum(w*x)
	LL <- sum(w*(mu*x - mu^2/2))
	tmp <- Optimality1.2A(x,w,mu)
	if (length(tmp$tnew) > 0){
		ht0 <- max(tmp$htnew)
		t0 <- tmp$tnew[which.max(tmp$htnew)]
	}else{
		ht0 <- 0
	}
	
	if (ht0 < delta2){
		return(list(tau=NULL,theta=mu,x=x,w=w,m=0,
			LL=LL,nr.local.search=0,nr.Newton=0))
	}
	
	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0

	tau <- NULL
	theta <- mu
	
	while (max(ht0) > delta2){
		nr.local.search <- nr.local.search + 1
		# Store current value of LL:
		LL.old <- LL
		# Add new knots
		tmp <- LocalNewKnots.2A(tau,theta,t0,ht0,delta2)
		tau <- tmp$knots
		theta <- tmp$theta
		isnewknot <- tmp$isnewknot
		
		wt <- LocalLinearSplines2.2A(tau,x,w)
		proposal <- LocalNewton.2A(wt,tau,theta,delta1)
		nr.Newton <- nr.Newton + 1
		while (proposal$dirderiv > delta1){
			theta.new <- proposal$theta.new
			# 2nd step size correction:
			corr <- LocalStepSize2.2A(tau,theta,theta.new,isnewknot)
			theta <- corr$theta
			if (length(corr$tau) < length(tau)){
				LL.ref <- -Inf
				# Update of knots, wt and isnewknot:
				tau <- corr$tau
				isnewknot <- corr$isnewknot
				wt <- LocalLinearSplines2.2A(tau,x,w)
			}else{
				LL.ref <- LL
			}
			# Normalization:
			theta <- LocalNormalize.2A(tau,theta)
			LL <- sum(wt*theta)
			if (LL > LL.ref){
				proposal <- LocalNewton.2A(wt,tau,theta,delta1)
				nr.Newton <- nr.Newton + 1
			}else{
				proposal$dirderiv <- 0
			}
		}
		
		# Check for global optimality:
		if (LL > LL.old){
			tmp <- Optimality2.2A(x,w,tau,theta)
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
			kk <- (ht0 >= max(delta2, max(ht0)*0.0001))
			t0 <- t0[kk]
			if (length(t0) > 0){
				ht0 <- ht0[kk]
			}else{
				ht0 <- 0
			}
		}else{
			ht0 <- 0
		}
	}
	return(list(tau=tau,theta=theta,x=x,w=w,m=length(tau),
		LL=LL,
		nr.local.search=nr.local.search,nr.Newton=nr.Newton))
}

### Auxiliary programs 4 ----
### Checking for new potential knots 

LocalNewKnots.2A <- function(tau,theta,t0,ht0,delta=0)
# Deactivation of constraints:
# We look for new knots on the left of the left-most knot, on the right of the
# right-most knot, and in between any consecutive two knots.
{
	if (is.null(tau)){
		k <- which.max(ht0)
		return(list(knots=t0[k],
			theta=c(theta,theta*t0[k] - theta^2/2,theta),
			isnewknot=TRUE))
	}
	
	m <- length(tau)
	knots.new <- rep(NA,2*m+1)
	knots.new[2*(1:m)] <- tau
	theta.new <- rep(NA,2*m+3)
	theta.new[2*(0:(m+1))+1] <- theta
	isnewknot <- rep(FALSE,2*m+1)
	
	# (i) To the left of the left-most knot
	II <- which(t0 < tau[1])
	if (length(II) > 0 && max(ht0[II]) > delta) {
		k <- which.max(ht0[II])
		i <- II[k]
		knots.new[1] <- t0[i]
		isnewknot[1] <- TRUE
		theta.new[2] <- theta[2] + (t0[i] - tau[1])*theta[1]
	}
	# (ii) To the right of the right-most knot
	II <- which(t0 > tau[m])
	if (length(II) > 0 && max(ht0[II]) > delta) {
		k <- which.max(ht0[II])
		i <- II[k]
		knots.new[2*m+1] <- t0[i]
		isnewknot[2*m+1] <- TRUE
		theta.new[2*m+2] <- theta[m+1] +
			(t0[i] - tau[m])*theta[m+2]
	}
	# (iii) In between consecutive knots tau[j] and tau[j+1]
	if (m > 1) {
		for (j in 1:(m-1)) {
			II <- which(t0 > tau[j] & t0 < tau[j+1])
			if (length(II) > 0 && max(ht0[II]) > delta) {
				k <- which.max(ht0[II])
				i <- II[k]
				knots.new[2*j+1] <- t0[i]
				isnewknot[2*j+1] <- TRUE
				theta.new[2*j+2] <- ((tau[j+1] - t0[i])*theta[j+1] +
						(t0[i] - tau[j])*theta[j+2])/
					(tau[j+1] - tau[j])
			}
		}
	}
	
	tmp <- !is.na(knots.new)
	knots.new <- knots.new[tmp]
	isnewknot <- isnewknot[tmp]
	theta.new <- theta.new[!is.na(theta.new)]
	return(list(knots=knots.new,
		theta=theta.new,
		isnewknot=isnewknot))
}


LocalDirDeriv1.2A <- function(t,x,w,mu=sum(x*w))
# This procedure computes the directional derivatives
#    h(t[i]) := DL(theta,V_{t[i]})
# for 1 <= i <= length(t) in case of
#    theta(x) = mu*x - mu^2/2 .
# In addition, it returns the right-sided derivatives
#    h'(t[i]+) .
# It is used to illustrate and check the methods.
# The estimation procedure uses only Optimality1.2A().
{
	z <- c(x,t)
	wz <- c(w,rep(0,length(t)))
	z.is.t <- c(rep(FALSE,length(x)),rep(TRUE,length(t)))
	ord <- order(z)
	z <- z[ord]
	wz <- wz[ord]
	z.is.t <- z.is.t[ord]
	nz <- length(z)
	h1t <- cumsum(wz) - pnorm(z - mu)
	h0t <- z*h1t - cumsum(z*wz) + mu*pnorm(z-mu) - dnorm(z-mu)
	return(cbind("h(t)"=h0t[z.is.t],"h'(t+)"=h1t[z.is.t]))
}

Optimality1.2A <- function(x,w,mu=sum(w*x))
# Determine knots t0 with locally maximal directional
# derivatives DL(theta,V_{t0}) in case of
# theta(x) = mu*x - mu^2/2.
{
	n <- length(x)
	csw <- cumsum(w[1:(n-1)])
	t0 <- mu + qnorm(csw)
	ht0 <- rep(NA,n-1)
	ht0a <- cumsum(x[n:2]*w[n:2])[(n-1):1]
	ht0b <- cumsum(w[n:2])[(n-1):1]
	tmp <- (x[1:(n-1)] < t0 & t0 < x[2:n])
	if (any(tmp)){
		ht0[tmp] <- ht0a[tmp] -
			ht0b[tmp]*t0[tmp] -
			(mu - t0[tmp])*pnorm(mu - t0[tmp]) -
			dnorm(mu - t0[tmp])
	}
	return(list(tnew=t0[tmp],htnew=ht0[tmp]))
}


LocalDirDeriv2.2A <- function(t,x,w,tau,theta)
# This procedure computes the directional derivatives
#    h(t[i]) := DL((tau,theta),V_{t[i],(tau,theta)})
# and the right-sided derivatives
#    h'(t[i]+)
# for 1 <= i <= length(t).
# It is used to illustrate and check the methods.
# The estimation procedure uses only Optimality2.2A().
{
	m <- length(tau)
	z <- c(x,tau,t)
	wz <- c(w,rep(0,m),rep(0,length(t)))
	z.is.t <- c(rep(FALSE,length(x)+m),rep(TRUE,length(t)))
	ord <- order(z)
	z <- z[ord]
	wz <- wz[ord]
	z.is.t <- z.is.t[ord]
	nz <- length(z)
	h0t <- rep(0,nz)
	h1t <- rep(0,nz)
	Phatz <- rep(0,nz)
	Pthetaz <- rep(0,nz)
	Yhatz <- rep(0,nz)
	thetaz <- rep(0,nz)

	# Points z[i] < tau[1]:
	ib <- 1
	while (z[ib+1] <= tau[1]){
		ib <- ib+1
	}
		# Now, z[ib] = tau[1]
	thetaz[1:ib] <- theta[2] + (z[1:ib] - tau[1])*theta[1]
	Pthetaz[1:ib] <- K0.2A(thetaz[1:ib],-theta[1],-z[1:ib])
	Phatz[1:ib] <- cumsum(wz[1:ib])
	Yhatz[(ib-1):1] <- cumsum(wz[ib:2]*(tau[1] - z[ib:2]))
	h1t[1:ib] <- Phatz[1:ib] - Pthetaz[1:ib]
	h0t[1:ib] <- -Yhatz[1:ib] + (z[1:ib] - tau[1])*
		(h1t[1:ib] - J10.2A(thetaz[1:ib],theta[2],z[1:ib],tau[1]))
	# Points z[i] in [tau[1], tau[m]):
	if (m > 1){
		for (j in 1:(m-1)){
			ia <- ib
			while (z[ib+1] <= tau[j+1]){
				ib <- ib+1
			}
				# Now, z[ia] = tau[j] and z[ib] = tau[j+1]
			dj <- tau[j+1] - tau[j]
			Hj <- sum(wz[ia:ib]*(tau[j+1] - z[ia:ib]))/dj -
				J10.2A(theta[j+1],theta[j+2],tau[j],tau[j+1])
			thetaz[ia:ib] <- ((tau[j+1] - z[ia:ib])*theta[j+1] +
					(z[ia:ib] - tau[j])*theta[j+2])/dj
			Pthetaz[ia:ib] <- J00.2A(theta[j+1],thetaz[ia:ib],
				tau[j],z[ia:ib])
			Phatz[ia:ib] <- cumsum(wz[ia:ib])
			Yhatz[ia:ib] <- cumsum(wz[ia:ib]*(z[ia:ib] - tau[j]))
			h1t[ia:ib] <- Phatz[ia:ib] - Pthetaz[ia:ib] - Hj
			h0t[ia:ib] <- -Yhatz[ia:ib] +
				(z[ia:ib] - tau[j])*(h1t[ia:ib] +
					J01.2A(theta[j+1],thetaz[ia:ib],tau[j],z[ia:ib]))
		}
	}
	# Points z[i] >= tau[m]:
	ia <- ib
		# Now, z[ia] = tau[m]
	thetaz[ia:nz] <- theta[m+1] + (z[ia:nz] - tau[m])*theta[m+2]
	Pthetaz[ia:nz] <- K0.2A(thetaz[ia:nz],theta[m+2],z[ia:nz])
	Phatz[nz:ia] <- c(0,cumsum(wz[nz:(ia+1)]))
	Yhatz[ia:nz] <- cumsum(wz[ia:nz]*(z[ia:nz] - tau[m]))
	h1t[ia:nz] <- Pthetaz[ia:nz] - Phatz[ia:nz]
	h0t[ia:nz] <- -Yhatz[ia:nz] + (z[ia:nz] - tau[m])*
			(h1t[ia:nz] +
				J01.2A(theta[m+1],thetaz[ia:nz],tau[m],z[ia:nz]))
	return(cbind("h(t)"=h0t[z.is.t],"h'(t+)"=h1t[z.is.t]))
}

Optimality2.2A <- function(x,w,tau,theta)
# Determine knots t0 with locally maximal directional
# derivatives DL(theta,V_{t0,theta}) in case of a
# non-affine current fit given by (tau,theta).
{
	# Determine vector t of boundary points of open intervals
	# containing no data points or knots but potentially a new
	# knot. In addition, store the empirical weights of these
	# points:
	t <- c(x,tau)
	wt <- c(w,rep(0,length(tau)))
	ord <- order(t)
	t <- t[ord]
	wt <- wt[ord]
	nt <- length(t)
	# Vector of potential new knots and vector of
	# corresponding directional derivatives:
	tnew <- rep(NA,nt-1)
	htnew <- rep(NA,nt-1)

	# Auxiliary vectors:
	h1tp <- rep(NA,nt-1)
	h1tm <- rep(NA,nt-1)
		# Vectors of one-sided derivatives
		# h1tp[i] = h_theta'(t[i]+) and
		# h1tm[i] = h_theta'(t[i+1]-)
	thetat <- rep(NA,nt)
		# theta, evaluated at t
	thetastar <- rep(NA,nt-1)
		# theta, evaluated at potential new knots
	Phatt <- rep(0,nt)
		# Empirical probabilities of certain intervals
	Pthetat <- rep(0,nt)
		# Estimated probabilities of certain intervals
	Yhatt <- rep(NA,nt)
		# Further empirical integrals

	m <- length(tau)

	# Search for potential new knots to the left of tau[1]: 
	ib <- 1
	while (t[ib+1] <= tau[1]){
		ib <- ib+1
	}
		# Now, we consider potential new knots between
		# t[1] = x[1] and t[ib] = tau[1].
	Phatt[1:(ib-1)] <- cumsum(wt[1:(ib-1)])
		# Phatt[i] = empirical prob. of (-Inf,t[i]]
	thetat[1:ib] <- theta[2] + (t[1:ib] - tau[1])*theta[1]
		# thetat[i] = theta(t[i])
	Pthetat[1:ib] <- K0.2A(thetat[1:ib],-theta[1],-t[1:ib])
		# Pthetat[i] = estimated prob. of (-Inf,t[i]]
	Yhatt[(ib-1):1] <- cumsum(wt[ib:2]*(tau[1]-t[ib:2]))
		# Yhatt[i] = empirical integral of (tau[1] - x)
		# over (t(i),tau[1]]
	h1tp[1:(ib-1)] <- Phatt[1:(ib-1)] - Pthetat[1:(ib-1)]
	h1tm[1:(ib-1)] <- Phatt[1:(ib-1)] - Pthetat[2:ib]
	JJ <- which(h1tp[1:(ib-1)] > 0 & h1tm[1:(ib-1)] < 0)
		# Which intervals (t[i],t[i+1]) contain
		# a local maximum of h_theta?
	tnew[JJ] <- InverseK0.2A(theta[2],theta[1],tau[1],Phatt[JJ],-1)
	failures <- (tnew[JJ] <= t[JJ] | tnew[JJ] >= t[JJ+1])
	if (any(failures)){
		tnew[JJ[failures]] <- NA
		JJ <- JJ[!failures]
	}
	thetastar[JJ] <- theta[2] + (tnew[JJ] - tau[1])*theta[1]
	htnew[JJ] <- -Yhatt[JJ] + (tau[1] - tnew[JJ])*
			J10.2A(thetastar[JJ],theta[2],tnew[JJ],tau[1])

	if (m > 1){
		for (j in 1:(m-1)){
			# Search for potential new knots between tau[j] and tau[j+1]:
			ia <- ib
			while (t[ib+1] <= tau[j+1]){
				ib <- ib+1
			}
				# Now, t[ia] = tau[j] and t[ib] = tau[j+1].
			Phatt[ia:(ib-1)] <- cumsum(wt[ia:(ib-1)])
				# Phatt[i] = empirical prob. of (tau[j],t[i]]
			dj <- tau[j+1] - tau[j]
			Hj <- sum(wt[ia:ib]*(tau[j+1] - t[ia:ib]))/dj -
				J10.2A(theta[j+1],theta[j+2],tau[j],tau[j+1])
				# Hj = empirical minus estimated integral of
				# j_{10}(x; tau[j],tau[j+1]) over (tau[j],tau[j+1]]
			thetat[ia:ib] <- ((tau[j+1] - t[ia:ib])*theta[j+1] +
					(t[ia:ib] - tau[j])*theta[j+2])/dj
				# thetat[i] = theta(t[i])
			Pthetat[ia:ib] <- J00.2A(theta[j+1],thetat[ia:ib],
				tau[j],t[ia:ib])
				# Pthetat[i] = estimated prob. of (tau[j],t[i]]
			Yhatt[ia:ib] <- cumsum(wt[ia:ib]*(t[ia:ib] - tau[j]))
				# Yhatt[i] = empirical integral of (x - tau[j])
				# over (tau[j],t[i]]
			h1tp[ia:(ib-1)] <- Phatt[ia:(ib-1)] -
				Pthetat[ia:(ib-1)] - Hj
			h1tm[ia:(ib-1)] <- Phatt[ia:(ib-1)] -
				Pthetat[(ia+1):ib] - Hj
			
			JJ <- ia-1 + which(h1tp[ia:(ib-1)] > 0 & h1tm[ia:(ib-1)] < 0)
				# Which intervals (t[i],t[i+1]) contain
				# a local maximum of h_theta?
			theta1j <- (theta[j+2] - theta[j+1])/dj
			tnew[JJ] <- InverseJ00.2A(theta[j+1],theta1j,
				tau[j],Phatt[JJ]-Hj)
			failures <- (tnew[JJ] <= t[JJ] | tnew[JJ] >= t[JJ+1])
			if (any(failures)){
				tnew[JJ[failures]] <- NA
				JJ <- JJ[!failures]
			}
			thetastar[JJ] <- ((tau[j+1] - tnew[JJ])*theta[j+1] +
				(tnew[JJ] - tau[j])*theta[j+2])/dj
			htnew[JJ] <- -Yhatt[JJ] + (tnew[JJ] - tau[j])*
				J01.2A(theta[j+1],thetastar[JJ],tau[j],tnew[JJ])
		}
	}
	
	# Search for potential new knots to the right of tau[m]: 
	ia <- ib
	ib <- nt
		# Now, t[ia] = tau[m] and t[ib] = x[n].
	Phatt[(ib-1):ia] <- cumsum(wt[ib:(ia+1)])
		# Phatt[i] = empirical prob. of (t[i],Inf)
	thetat[ia:ib] <- theta[m+1] + (t[ia:ib] - tau[m])*theta[m+2]
		# thetat[i] = theta(t[i])
	Pthetat[ia:ib] <- K0.2A(thetat[ia:ib],theta[m+2],t[ia:ib])
		# Pthetat[i] = estimated prob. of (t[i],Inf)
	Yhatt[ia:ib] <- cumsum(wt[ia:ib]*(t[ia:ib] - tau[m]))
		# Yhatt[i] = empirical integral of (x - tau[m])
		# over (tau[m],t[i]]
	h1tp[ia:(ib-1)] <- Pthetat[ia:(ib-1)] - Phatt[ia:(ib-1)]
	h1tm[ia:(ib-1)] <- Pthetat[(ia+1):ib] - Phatt[ia:(ib-1)]
	JJ <- ia-1 + which(h1tp[ia:(ib-1)] > 0 & h1tm[ia:(ib-1)] < 0)
		# Which intervals (t[i],t[i+1]) contain
		# a local maximum of h_theta?
	tnew[JJ] <- InverseK0.2A(theta[m+1],theta[m+2],
		tau[m],Phatt[JJ],+1)
	failures <- (tnew[JJ] <= t[JJ] | tnew[JJ] >= t[JJ+1])
	if (any(failures)){
		tnew[JJ[failures]] <- NA
		JJ <- JJ[!failures]
	}
	thetastar[JJ] <- theta[m+1] + (tnew[JJ] - tau[m])*theta[m+2]
	htnew[JJ] <- -Yhatt[JJ] + (tnew[JJ] - tau[m])*
			J01.2A(theta[m+1],thetastar[JJ],tau[m],tnew[JJ])

	return(list(tnew=tnew[!is.na(tnew)],htnew=htnew[!is.na(tnew)]))
}


### Auxiliary programs 3 ----
### 2nd step size correction

LocalConvexity.2A <- function(tau,theta)
# Changes of slope,
# all of which should be > 0 if theta defines
# a strictly convex function on {tau[i] : ...}.
{
	m <- length(tau)
	if (m >= 2){
		dtau <- tau[2:m] - tau[1:(m-1)]
		dtheta <- theta[3:(m+1)] - theta[2:m]
		dtdt <- c(theta[1],dtheta/dtau,theta[m+2])
	}
	else{
		dtdt <- c(theta[1],theta[3])
	}
	chofslope <- dtdt[2:(m+1)] - dtdt[1:m]
	return(chofslope)
}

LocalStepSize2.2A <- function(tau,theta,theta.new,
                              isnewknot=rep(FALSE, length(tau)))
# Replace theta with
#    (1 - t0)*theta + t0*theta.new ,
# where t0 is the largest number in [0,1] such that
# the new theta defines still a convex function.
# Return new tuples tau and theta, possibly shorter
# than the original ones.
{
	m <- length(tau)
	chofsl.new <- LocalConvexity.2A(tau,theta.new)
	# 1st case, should not happen!
	if (m == 1 && chofsl.new <= 0){
	  # In this case, the largest t0 such that phi = (1-t0)*theta + t0*theta.new
	  # is feasible is such that phi is a straight line. This should not happen
	  # at this stage, since we excluded this option in the beginning of the
	  # algorithm.
	  stop("m == 1 and chofsl.new <= 0 should not happen at this stage.")
	}
	# 2nd easy (in fact, ideal) case:
	if (min(chofsl.new) > 0){
	  return(list(tau=tau,theta=theta.new))
	}
	# All other cases with min(chofsl.new) <= 0 lead to
	# removal of one knot / activation of one constraint:
	#	
	# (1) If at least one of the new knots is responsible for a negative change of
	# slope in theta.new, this means that the constraint at the worst knot has to
	# be activated, resulting in the deletion of that knot. Then theta.new is 
	# discarded and will have to be recomputed with a Newton step.
	if (any(isnewknot) && min(chofsl.new[isnewknot]) <= 0) {
	  tmp  <- which(chofsl.new <= 0 & isnewknot)
	  j0 <- which.min(chofsl.new[tmp])
	  j0 <- tmp[j0]
	  theta <- theta[-(j0+1)]
	  tau <- tau[-j0]
	  isnewknot <- isnewknot[-j0]
	  return(list(tau=tau,theta=theta,isnewknot=isnewknot))
	}
	# (2) At this stage, none of the newly added knots is responsible for a 
	# negative change of slope in theta.new. Nevertheless, at least one of the old
	# constraints has to be activated and the corresponding knot has to be
	# deleted. We therefore look for the largest t0 such that
	# (1 - t0)*theta + t0*theta.new is feasible.
	chofsl <- LocalConvexity.2A(tau,theta)
	tmp  <- which(chofsl.new <= 0)
	t0 <- chofsl[tmp]/(chofsl[tmp] - chofsl.new[tmp])
	j0 <- which.min(t0)
	t0 <- t0[j0]
	j0 <- tmp[j0]
	theta <- (1 - t0)*theta + t0*theta.new
	tau <- tau[-j0]
	isnewknot <- isnewknot[-j0]
	theta <- theta[-(j0+1)]
	m <- m-1
	# Now make sure that phi is even strictly concave at each knot:
	chofsl <- LocalConvexity.2A(tau,theta)
	while (m > 1 && min(chofsl) <= 0){
	  j0 <- which.min(chofsl)
	  tau <- tau[-j0]
	  isnewknot <- isnewknot[-j0]
	  theta <- theta[-(j0+1)]
	  m <- m-1
	  chofsl <- LocalConvexity.2A(tau,theta)
	}
	if (m == 1 && chofsl <= 0){
	  theta[c(1,3)] <- mean(theta[c(1,3)])
	}
	return(list(tau=tau,theta=theta,isnewknot=isnewknot))
}

### Auxiliary programs 2 ----
### Normalization, log-likelihood plus derivatives
### and Newton step with 1st step size correction

LocalNormalize.2A <- function(tau,theta)
# Normalize theta such that it defines
# a log-probability density
{
	m <- length(tau)
	integral <- K0.2A(theta[2],-theta[1],-tau[1]) +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1)
	{
		integral <- integral +
			sum(J00.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m]))
	}
	theta.new <- theta
	theta.new[2:(m+1)] <- theta.new[2:(m+1)] - log(integral)
	return(theta.new)
}

LocalLL.2A <- function(wt,tau,theta)
# Log-likelihood
{
	m <- length(tau)
	integral <- K0.2A(theta[2],-theta[1],-tau[1]) +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1)
	{
		integral <- integral +
			sum(J00.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m]))
	}
	LL <- sum(wt*theta) + 1 - integral
	return(LL)
}

LocalLL2.2A <- function(wt,tau,theta)
# Log-likelihood plus its gradient vector
# and minus Hessian matrix
{
	m <- length(tau)
	integral <- K0.2A(theta[2],-theta[1],-tau[1]) +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		integral <- integral +
			sum(J00.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m]))
	}
	LL <- sum(wt*theta) + 1 - integral
	GLL <- wt
	GLL[1] <- GLL[1] + K1.2A(theta[2], -theta[1], -tau[1])
	GLL[2] <- GLL[2] - K0.2A(theta[2], -theta[1], -tau[1])
	GLL[m+1] <- GLL[m+1] - K0.2A(theta[m+1], theta[m+2], tau[m])
	GLL[m+2] <- GLL[m+2] - K1.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		GLL[2:m] <- GLL[2:m] -
			J10.2A(theta[2:m],theta[3:(m+1)],tau[1:(m-1)],tau[2:m])
		GLL[3:(m+1)] <- GLL[3:(m+1)] -
			J01.2A(theta[2:m],theta[3:(m+1)],tau[1:(m-1)],tau[2:m])
	}
	MHLL <- matrix(0,m+2,m+2)
	MHLL[1,1] <- K2.2A(theta[2], - theta[1], -tau[1])
	MHLL[1,2] <- - K1.2A(theta[2], - theta[1], -tau[1])
	MHLL[2,1] <- MHLL[1,2]
	MHLL[2,2] <- K0.2A(theta[2], -theta[1], -tau[1])
	MHLL[m+1,m+1] <- MHLL[m+1,m+1] +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	MHLL[m+1,m+2] <- K1.2A(theta[m+1], theta[m+2], tau[m])
	MHLL[m+2,m+1] <- MHLL[m+1,m+2]
	MHLL[m+2,m+2] <- K2.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		A20 <- cbind(2:m,2:m)
		MHLL[A20] <- MHLL[A20] +
			J20.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m])
		A02 <- cbind(3:(m+1),3:(m+1))
		MHLL[A02] <- MHLL[A02] +
			J02.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m])
		tmp11 <- J11.2A(theta[2:m],theta[3:(m+1)],
			tau[1:(m-1)],tau[2:m])
		MHLL[cbind(2:m,3:(m+1))] <- tmp11
		MHLL[cbind(3:(m+1),2:m)] <- tmp11
	}
	return(list(LL=LL,GLL=GLL,MHLL=MHLL))
}

LocalNewton.2A <- function(wt,tau,theta,
	delta0=10^(-11))
# Newton step with (1st) step size correction
{
	m <- length(tau)
	integral <- K0.2A(theta[2],-theta[1],-tau[1]) +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		integral <- integral +
			sum(J00.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m]))
	}
	LL <- sum(wt*theta) + 1 - integral
	GLL <- wt
	GLL[1] <- GLL[1] + K1.2A(theta[2], -theta[1], -tau[1])
	GLL[2] <- GLL[2] - K0.2A(theta[2], -theta[1], -tau[1])
	GLL[m+1] <- GLL[m+1] - K0.2A(theta[m+1], theta[m+2], tau[m])
	GLL[m+2] <- GLL[m+2] - K1.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		GLL[2:m] <- GLL[2:m] -
			J10.2A(theta[2:m],theta[3:(m+1)],tau[1:(m-1)],tau[2:m])
		GLL[3:(m+1)] <- GLL[3:(m+1)] -
			J01.2A(theta[2:m],theta[3:(m+1)],tau[1:(m-1)],tau[2:m])
	}
	MHLL <- matrix(0,m+2,m+2)
	MHLL[1,1] <- K2.2A(theta[2], - theta[1], -tau[1])
	MHLL[1,2] <- - K1.2A(theta[2], - theta[1], -tau[1])
	MHLL[2,1] <- MHLL[1,2]
	MHLL[2,2] <- K0.2A(theta[2], -theta[1], -tau[1])
	MHLL[m+1,m+1] <- MHLL[m+1,m+1] +
		K0.2A(theta[m+1], theta[m+2], tau[m])
	MHLL[m+1,m+2] <- K1.2A(theta[m+1], theta[m+2], tau[m])
	MHLL[m+2,m+1] <- MHLL[m+1,m+2]
	MHLL[m+2,m+2] <- K2.2A(theta[m+1], theta[m+2], tau[m])
	if (m > 1){
		A20 <- cbind(2:m,2:m)
		MHLL[A20] <- MHLL[A20] +
			J20.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m])
		A02 <- cbind(3:(m+1),3:(m+1))
		MHLL[A02] <- MHLL[A02] +
			J02.2A(theta[2:m],theta[3:(m+1)],
				tau[1:(m-1)],tau[2:m])
		tmp11 <- J11.2A(theta[2:m],theta[3:(m+1)],
			tau[1:(m-1)],tau[2:m])
		MHLL[cbind(2:m,3:(m+1))] <- tmp11
		MHLL[cbind(3:(m+1),2:m)] <- tmp11
	}
	dtheta <- qr.solve(MHLL,GLL)
	dirderiv <- sum(dtheta * GLL)
	theta.new <- theta + dtheta
	LL.new <- LocalLL.2A(wt,tau,theta.new)
	while (LL.new < LL + dirderiv/3 && dirderiv > delta0){
		theta.new <- (theta + theta.new)/2
		dirderiv <- dirderiv/2
		LL.new <- LocalLL.2A(wt,tau,theta.new)
	}
	return(list(theta.new=theta.new,dirderiv=dirderiv,
		LL=LL,LL.new=LL.new))
}


### Auxiliary programs 1 ----
### Linear splines and
### Basic functions J(.,.), K(.,.) plus 1st and 2nd
### partial derivatives

LocalLinearSplines1.2A <- function(tau,x=tau)
# For a vector tau with m >= 1 strictly increasing
# components, this function computes a design matrix Dx
# for the set of piecewise linear functions f on the real
# line with kinks only on {tau[j] : 1 <= j <= m}, evaluated
# at x. The basis corresponds to the parameter vector theta
# with components
#    theta[1]   = f'(tau[1] -) ,
#    theta[j+1] = f(tau[j]), 1 <= j <= m ,
#    theta[m+2] = f'(tau[m] +) .
{
	m <- length(tau)
	n <- length(x)
	Dx <- matrix(0,n,m+2)
	tmp <- (x <= tau[1])
	Dx[tmp,1] <- x[tmp] - tau[1]
	Dx[tmp,2] <- 1
	if (m > 1){
		dtau <- tau[2:m] - tau[1:(m-1)]
		for (j in 1:(m-1)){
			tmp <- (x > tau[j] & x <= tau[j+1])
			Dx[tmp,j+1] <- (tau[j+1] - x[tmp])/dtau[j]
			Dx[tmp,j+2] <- (x[tmp] - tau[j])/dtau[j]
		}
	}
	tmp <- (x > tau[m])
	Dx[tmp,m+1] <- 1
	Dx[tmp,m+2] <- x[tmp] - tau[m]
	return(Dx)
}


LocalLinearSplines2.2A <- function(tau,
	x=tau,w=rep(1/length(x),length(x)))
# For linear splines with given knots
#   tau[1] < tau[2] < ... < tau[m]
# and arbitrary vectors x and w of equal length with w >= 0,
# this procedure returns a weight vector wt, such that for
# the basis functions f_1, f_2, ..., f_m, f_{m+1}, f_{m+2}
# described/computed in LocalLinearSplines1.2A(),
#    wt[j] = sum_{i=1}^n w[i]*f_j(x[i]) .
{
	Dx <- LocalLinearSplines1.2A(tau,x)
	wt <- colSums(Dx * w)
	return(wt)
}


J00.2A <- function(x,y,a=0,b=1)
# J00.2A(x,y,a,b)
#    =  integral_a^b
#          exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <- (y - x)/(b - a + (a == b))
	at <- a - yt
	bt <- b - yt
	mt <- (a + b)/2 - yt
	dt <- (b - a)/2
	pdiff <- pnorm(bt) - pnorm(at)
	tmp <- (at > 0)
	pdiff[tmp] <- pnorm(-at[tmp]) - pnorm(-bt[tmp])
	J00xy <- exp(xt + yt^2/2 +
		pmax(log(pdiff),-mt^2/2 + log(2*pnorm(dt) - 1)))
	return(J00xy)
}

J10.2A <- function(x,y,a=0,b=1)
# J10.2A(x,y,a,b)
#    =  integral_a^b
#          (1-v(z))*exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <- (y - x)/(b - a + (a == b))
	at <- a - yt
	bt <- b - yt
	pdiff <- pnorm(bt) - pnorm(at)
	tmp <- (at > 0)
	pdiff[tmp] <- pnorm(-at[tmp]) - pnorm(-bt[tmp])
	J10xy <- exp(xt + yt^2/2)*
		(bt*pdiff + dnorm(bt) - dnorm(at))/
			(b - a + (a == b))
	return(J10xy)
}

J01.2A <- function(x,y,a=0,b=1)
# J01.2A(x,y,a,b)
#    =  integral_a^b
#          v(z)*exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
# Symmetry considerations reveal that
#    J01.2A(x,y,a,b) = J10.2A(y,x,-b,-a)
{
	return(J10.2A(y,x,-b,-a))
}

J20.2A <- function(x,y,a=0,b=1)
# J20.2A(x,y,a,b)
#    =  integral_a^b
#          (1-v(z))^2*exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <- (y - x)/(b - a + (a == b))
	at <- a - yt
	bt <- b - yt
	pdiff <- pnorm(bt) - pnorm(at)
	tmp <- (at > 0)
	pdiff[tmp] <- pnorm(-at[tmp]) - pnorm(-bt[tmp])
	J20xy <- exp(xt + yt^2/2)*
		((1 + bt^2)*pdiff +
		 (at - 2*bt)*dnorm(at) + bt*dnorm(bt))/
			(b - a + (a == b))^2
	return(J20xy)
}

J11.2A <- function(x,y,a=0,b=1)
# J11.2A(x,y,a,b)
#    =  integral_a^b
#          (1-v(z))*v(z)*exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <- (y - x)/(b - a + (a == b))
	at <- a - yt
	bt <- b - yt
	pdiff <- pnorm(bt) - pnorm(at)
	tmp <- (at > 0)
	pdiff[tmp] <- pnorm(-at[tmp]) - pnorm(-bt[tmp])
	J11xy <- exp(xt + yt^2/2)*
		(-(1 + at*bt)*pdiff + bt*dnorm(at) - at*dnorm(bt))/
			(b - a + (a == b))^2
	return(J11xy)
}

J02.2A <- function(x,y,a=0,b=1)
# J02.2A(x,y,a,b)
#    =  integral_a^b
#          v(z)^2*exp((1-v(z))x + v(z)y)*dnorm(z) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <- (y - x)/(b - a + (a == b))
	at <- a - yt
	bt <- b - yt
	pdiff <- pnorm(bt) - pnorm(at)
	tmp <- (at > 0)
	pdiff[tmp] <- pnorm(-at[tmp]) - pnorm(-bt[tmp])
	J02xy <- exp(xt + yt^2/2)*
		((1 + at^2)*pdiff +
		 (2*at - bt)*dnorm(bt) - at*dnorm(at))/
			(b - a + (a == b))^2
	return(J02xy)
}


K0.2A <- function(x,y,a=0)
# K0.2A(x,y,a) = integral_a^Inf exp(x + y*(z-a))*dnorm(z) dz
{
	return(exp(x - a*y + y^2/2)*pnorm(y-a))
}

K1.2A <- function(x,y,a=0)
# K1.2A(x,y,a) = dK0.2A(x,y,a)/dy
{
	return(exp(x - a*y + y^2/2)*
		((y - a)*pnorm(y-a) + dnorm(y - a)))
}

K2.2A <- function(x,y,a=0)
# K2.2A(x,y,a) = d^2K0.2A(x,y,a)/dy^2
{
	return(exp(x - a*y + y^2/2)*
		((1 + (y - a)^2)*pnorm(y-a) +
			(y - a)*dnorm(y - a)))
}


InverseK0.2A <- function(theta0, theta1, tau0, c, s=1)
# Determines the unique solution tau of
#    K0.2A(theta0 + theta1*(tau-tau0),s*theta1,s*tau) = c.
{
	c <- pmax(c,0)
	d <- c*exp(-theta0 + theta1*tau0 - theta1^2/2)
	d <- pmin(d,1)
	return(theta1 - s*qnorm(d))
}

InverseJ00.2A <- function(theta0, theta1, tau0, c)
# Determines the unique solution tau of
#    J00.2A(theta0, theta0 + theta1*(tau-tau0),tau0,tau) = c.
{
	d <- pnorm(tau0 - theta1) +
		c*exp(-theta0 + theta1*tau0 - theta1^2/2)
	d <- pmax(pmin(d,1),0)
	return(theta1 + qnorm(d))
}


LocalInterpolate <- function(x,d)
# For a vector x with n = length(x) >= 2, this function returns
# a vector xx with n*d - d + 1 components, such that
# xx[(j*d-d+1):(j*d)] == x[j] + (x[j+1] - x[j])*(1:d)/d
# for 1 <= j < n.
{
	n <- length(x)
	xx <- rep(NA,n*d-d+1)
	xx[1] <- x[1]
	for (k in 1:d){
		lambda <- k/d
		xx[1 + k + (0:(n-2))*d] <-
			(1 - lambda)*x[1:(n-1)] + lambda*x[2:n]
	}
	return(xx)
}

### Auxiliary program 0 ----

LocalPrepareData <- function(X,W)
# Replace a data vector X with weights W
# by a vector x with strictly increasing
# components such that {x[.]} = {X[.]}
# and corresponding weights w.
{
	n <- length(X)
	tmp <- order(X)
	X <- X[tmp]
	W <- W[tmp]
	XX <- c(X,Inf)
	WW <- cumsum(W)
	tmp <- which(XX[1:n] < XX[2:(n+1)])
	x <- X[tmp]
	n <- length(x)
	WW <- WW[tmp]
	w <- WW - c(0,WW[1:(n-1)])
	w <- w/WW[n]
	return(list(x=x,w=w,n=n))
}
