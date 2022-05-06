#######################################################
### Active Set Algorithm for the Computation of the ###
### NPMLE of a Log-Convex, Isotonic Density Ratio   ###
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

TailInflation.B0 <- function(X, W=rep(1/length(X),length(X)),
	alpha=1, beta=1,
	xmin=0, xmax=1.1*max(X),
	delta1=10^(-10)/length(X),
	delta2=10^(-3)*sqrt(alpha)/length(X))
# Computation of a log-convex and increasing density ratio with respect to
# the Gamma distribution with shape alpha and rate beta..
# Graphical display of all intermediate steps...
{
	# If X ~ P, where P(dx) = e^{theta(x)} P_{alpha,beta}(dx)
	# with P_{alpha,beta} the Gamma distribution with
	# shape alpha and rate beta, then beta*X ~ Pt, where
	# Pt(dx) = e^thetat(x) P_alpha,1(dx) and thetat(x) = theta(x/beta).
	# Hence we multiply the data by beta, estimate thetat and obtain
	# theta(x) = thetat(x*beta) in the end.
	X <- X*beta
	
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
	mu <- sum(w*x)
	
	if (mu <= alpha){
		tau <- NULL
		theta <- 0
		LL <- 0
		m <- 0
		tmp <- Optimality1.2B(x,w,mu,alpha)
		if (length(tmp$tnew) > 0){
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
		}else{
			ht0 <- 0
		}
	}else{
		tau <- 0
		kappa <- 1 - alpha/mu
		ckappa <- alpha*log(mu/alpha)
		theta <- c(-ckappa, kappa)
		LL <- -ckappa + kappa*mu
		m <- 1
		tmp <- Optimality2.2B(x,w,tau,theta,alpha)
		if (length(tmp$tnew) > 0){
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
		}else{
			ht0 <- 0
		}
	}
	
	# Graphical display:
	message <- readline("First plot:")
	if (m > 0){
		plot(c(xmin,xmax),rep(0,2),type='l',
			lwd=3,col='green',
			xlab=expression(italic(x)),
			ylab=expression(italic(theta(x))),
			xlim=range(x),
			ylim=range(kappa*c(xmin,xmax) - ckappa),
			main=paste('Starting point: ',
				'LL = ',round(LL,digits=9),sep=''))
		lines(c(xmin,xmax),kappa*c(xmin,xmax)-ckappa,lwd=2)
		rug(x)
	}else{
		plot(c(xmin,xmax),rep(0,2),type='l',
			lwd=3,col='green',
			xlab=expression(italic(x)),
			ylab=expression(italic(theta(x))),
			xlim=range(x),
			main=paste('Starting point: ',
				'LL = ',round(LL,digits=9),sep=''))
		lines(c(xmin,xmax),rep(0,2),lwd=2)
		rug(x)
	}
	
	if (max(ht0) <= delta2){
		print(paste('Parametric fit sufficient: Gamma(',
			alpha,',',beta*theta[m+1],
			'). Log-likelihood = ',
			LL,'.',sep=''))
		if (m==1){
			theta[2] <- theta[2]*beta
		}
		return(list(tau=tau,theta=theta,x=x,w=w,m=m,
			LL=LL,
			nr.local.search=0,
			nr.Newton=0))
	}
	
	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0
	
	while (max(ht0) > delta2){
		nr.local.search <- nr.local.search + 1
		# Store current value of LL:
		LL.old <- LL
		# Add new knots
		tmp <- LocalNewKnots.2B(tau,theta,t0,ht0,delta2)
		tau <- tmp$knots
		theta <- tmp$theta
		isnewknot <- tmp$isnewknot
		# Graphical display:
		message <- readline("New knots:")
		xtx <- c(xmin,tau,xmax)
		Dtmp <- LocalLinearSplines1.2B(tau,xtx)
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
		
		wt <- LocalLinearSplines2.2B(tau,x,w)
		proposal <- LocalNewton.2B(wt,tau,theta,alpha,delta1)
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
				main=paste('LL =',round(LL,digits=9)))
			lines(xtx,thetaxtx,lwd=2)
			lines(xtx,thetaxtx.new,lwd=2,col='blue')
			points(tau,thetaxtx.new[2:(length(tau)+1)],col='blue')
			abline(v=tau[isnewknot],col='blue')
			rug(x)
			
			# 2nd step size correction:
			corr <- LocalStepSize2.2B(tau,theta,theta.new,isnewknot)
			theta <- corr$theta
			if (length(corr$tau) < length(tau)){
				LL.ref <- -Inf
				# Update of knots, wt and isnewknot:
				tau <- corr$tau
				isnewknot <- corr$isnewknot
				wt <- LocalLinearSplines2.2B(tau,x,w)
				# Graphical display:
				message <- readline("2nd step size correction:")
				xtx <- c(xmin,tau,xmax)
				Dtmp <- LocalLinearSplines1.2B(tau,xtx)
				thetaxtx <- Dtmp %*% theta
				lines(xtx,thetaxtx,lwd=2,col='magenta')
			}else{
				LL.ref <- LL
			}
			
			# Normalization:
			theta <- LocalNormalize.2B(tau,theta,alpha)
			LL <- sum(wt*theta)
			if (LL > LL.ref){
				proposal <- LocalNewton.2B(wt,tau,theta,alpha,delta1)
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
			main=paste('LL =',round(LL,digits=9)))
		lines(xtx,thetaxtx,lwd=2)
		points(tau,thetaxtx[2:(length(tau)+1)])
		rug(x)
		
		# Check for global optimality:
		if (LL > LL.old){
			tmp <- Optimality2.2B(x,w,tau,theta,alpha)
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
			kk <- (ht0 >= max(delta2, max(ht0)*0.0001))
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
		main=paste('LL =',round(LL,digits=9)))
	lines(xtx,thetaxtx,lwd=2)
	points(tau,thetaxtx[2:(length(tau)+1)])
	rug(x)
	
	# Rescale to account for beta != 1
	x <- x/beta
	tau <- tau/beta
	theta[length(tau)+1] <- beta*theta[length(tau)+1]
	
	return(list(tau=tau,theta=theta,x=x,w=w,m=length(tau),
		LL=LL,
		nr.local.search=nr.local.search,
		nr.Newton=nr.Newton))
}

TailInflation.B <- function(X, W=rep(1/length(X),length(X)),
	alpha=1, beta=1,
	delta1=10^(-10)/length(X),
	delta2=10^(-3)*sqrt(alpha)/length(X))
# Computation of a log-convex and increasing density ratio with respect to
# the Gamma distribution with shape alpha and rate beta..
# Pure computation, no graphical displays...
{
	# If X ~ P, where P(dx) = e^theta(x) P_alpha,beta(dx)
	# with P_alpha,beta the Gamma distribution with
	# shape alpha and rate beta, then beta*X ~ Pt, where
	# Pt(dx) = e^thetat(x) P_alpha,1(dx) and thetat(x) = theta(x/beta).
	# Hence we multiply the data by beta, estimate thetat and obtain
	# theta(x) = thetat(x*beta) in the end.
	X <- X*beta
	
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
	mu <- sum(w*x)
	
	if (mu <= alpha){
		tau <- NULL
		theta <- 0
		LL <- 0
		m <- 0
		tmp <- Optimality1.2B(x,w,mu,alpha)
		if (length(tmp$tnew) > 0){
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
		}else{
			ht0 <- 0
		}
	}else{
		tau <- 0
		kappa <- 1 - alpha/mu
		ckappa <- alpha*log(mu/alpha)
		theta <- c(-ckappa, kappa)
		LL <- -ckappa + kappa*mu
		m <- 1
		tmp <- Optimality2.2B(x,w,tau,theta,alpha)
		if (length(tmp$tnew) > 0){
			ht0 <- tmp$htnew
			t0 <- tmp$tnew
		}else{
			ht0 <- 0
		}
	}
	
	if (max(ht0) <= delta2){
		if (m==1){
			theta[2] <- theta[2]*beta
		}
		return(list(tau=tau,theta=theta,x=x,w=w,m=m,
			LL=LL,
			nr.local.search=0,
			nr.Newton=0))
	}
	
	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0
	
	while (max(ht0) > delta2){
		nr.local.search <- nr.local.search + 1
		# Store current value of LL:
		LL.old <- LL
		# Add new knots
		tmp <- LocalNewKnots.2B(tau,theta,t0,ht0,delta2)
		tau <- tmp$knots
		theta <- tmp$theta
		isnewknot <- tmp$isnewknot
		wt <- LocalLinearSplines2.2B(tau,x,w)
		
		proposal <- LocalNewton.2B(wt,tau,theta,alpha,delta1)
		nr.Newton <- nr.Newton + 1
		while (proposal$dirderiv > delta1){
			theta.new <- proposal$theta.new
			# 2nd step size correction:
			corr <- LocalStepSize2.2B(tau,theta,theta.new,isnewknot)
			theta <- corr$theta
			if (length(corr$tau) < length(tau)){
				LL.ref <- -Inf
				# Update of knots, wt and isnewknot:
				tau <- corr$tau
				isnewknot <- corr$isnewknot
				wt <- LocalLinearSplines2.2B(tau,x,w)
			}else{
				LL.ref <- LL
			}
			# Normalization:
			theta <- LocalNormalize.2B(tau,theta,alpha)
			LL <- sum(wt*theta)
			if (LL > LL.ref){
				proposal <- LocalNewton.2B(wt,tau,theta,alpha,delta1)
				nr.Newton <- nr.Newton + 1
			}else{
				proposal$dirderiv <- 0	
			}
		}
		
		# Check for global optimality:
		if (LL > LL.old){
			tmp <- Optimality2.2B(x,w,tau,theta,alpha)
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
	
	# Rescale to account for beta != 1
	x <- x/beta
	tau <- tau/beta
	theta[length(tau)+1] <- beta*theta[length(tau)+1]
	
	return(list(tau=tau,theta=theta,x=x,w=w,m=length(tau),
		LL=LL,
		nr.local.search=nr.local.search,
		nr.Newton=nr.Newton))
}

### Auxiliary programs 4 ----
### Checking for knew potential knots

LocalNewKnots.2B <- function(tau,theta,t0,ht0,delta=0)
# Deactivation of constraints:
# We look for new knots on the left of the left-most knot, on the right of the
# right-most knot, and in between any consecutive two knots.
{
	if (is.null(tau)){
		k <- which.max(ht0)
		return(list(knots=t0[k],theta=c(0,0),isnewknot=TRUE))
	}
	m <- length(tau)
	if (tau[1] == 0) { # m potential new knots
		knots.new <- rep(NA,2*m)
		knots.new[2*(1:m)-1] <- tau
		theta.new <- rep(NA,2*m+1)
		theta.new[2*(0:m)+1] <- theta
		isnewknot <- rep(FALSE,2*m)
		
		# (i) To the right of the right-most knot
		II <- which(t0 > tau[m])
		if (length(II) > 0 && max(ht0[II]) > delta) {
			k <- which.max(ht0[II])
			i <- II[k]
			knots.new[2*m] <- t0[i]
			isnewknot[2*m] <- TRUE
			theta.new[2*m] <- theta[m] + (t0[i] - tau[m])*theta[m+1]
		}
		# (ii) In between consecutive knots tau[j] and tau[j+1]
		if (m > 1) {
			for (j in 1:(m-1)) {
				II <- which(t0 > tau[j] & t0 < tau[j+1])
				if (length(II) > 0 && max(ht0[II]) > delta){
					k <- which.max(ht0[II])
					i <- II[k]
					knots.new[2*j] <- t0[i]
					isnewknot[2*j] <- TRUE
					theta.new[2*j] <- ((tau[j+1] - t0[i])*theta[j] +
							(t0[i] - tau[j])*theta[j+1])/
						(tau[j+1] - tau[j])
				}
			}
		}
	} else { # m+1 potential new knots
		knots.new <- rep(NA,2*m+1)
		knots.new[2*(1:m)] <- tau
		theta.new <- rep(NA,2*m+2)
		theta.new[2*(1:(m+1))] <- theta
		isnewknot <- rep(FALSE,2*m+1)
		
		# (i) To the left of the left-most knot
		II <- which(t0 < tau[1])
		if (length(II) > 0 && max(ht0[II]) > delta) {
			k <- which.max(ht0[II])
			i <- II[k]
			knots.new[1] <- t0[i]
			isnewknot[1] <- TRUE
			theta.new[1] <- theta[1]
		}
		# (ii) To the right of the right-most knot
		II <- which(t0 > tau[m])
		if (length(II) > 0 && max(ht0[II]) > delta) {
			k <- which.max(ht0[II])
			i <- II[k]
			knots.new[2*m+1] <- t0[i]
			isnewknot[2*m+1] <- TRUE
			theta.new[2*m+1] <- theta[m] + (t0[i] - tau[m])*theta[m+1]
		}
		# (iii) In between consecutive knots tau[j] and tau[j+1]
		if (m > 1) {
			for (j in 1:(m-1)) {
				II <- which(t0 > tau[j] & t0 < tau[j+1])
				if (length(II) > 0 && max(ht0[II]) > delta){
					k <- which.max(ht0[II])
					i <- II[k]
					knots.new[2*j+1] <- t0[i]
					isnewknot[2*j+1] <- TRUE
					theta.new[2*j+1] <- ((tau[j+1] - t0[i])*theta[j] +
							(t0[i] - tau[j])*theta[j+1])/
						(tau[j+1] - tau[j])
				}
			}
		}
	}
	
	tmp <- !is.na(knots.new)
	knots.new <- knots.new[tmp]
	isnewknot <- isnewknot[tmp]
	theta.new <- theta.new[!is.na(theta.new)]
	return(list(knots=knots.new,theta=theta.new,isnewknot=isnewknot))
}


LocalDirDeriv1.2B <- function(t,x,w,alpha=1)
# This procedure computes the directional derivatives
#    h(t[i]) := DL(theta,V_{t[i]})
# for 1 <= i <= length(t) in case of
#    theta(x) = 0 .
# In addition, it returns the right-sided derivatives
#    h'(t[i]+) .
# It is used to illustrate the methods and check the algorithms.
# The estimation procedure uses only the function Optimality1.2B().
{
	z <- c(x,t)
	wz <- c(w,rep(0,length(t)))
	z.is.t <- c(rep(FALSE,length(x)),rep(TRUE,length(t)))
	ord <- order(z)
	z <- z[ord]
	wz <- wz[ord]
	z.is.t <- z.is.t[ord]
	nz <- length(z)
	
	h1t <- cumsum(wz) - pgamma(z,alpha,1)
	h0t <- z*h1t + sum(wz*z) - alpha - cumsum(wz*z) +
		alpha*pgamma(z,alpha+1,1)
	return(cbind("h(t)"=h0t[z.is.t],"h'(t+)"=h1t[z.is.t]))
}

Optimality1.2B <- function(x,w,mu=sum(w*x),alpha=1,
	mindiff=10^(-5))
# Determine knots t0 with locally maximal directional
# derivatives DL(theta,V_{t0}) in case of
# theta(x) = 0
{
	n <- length(x)
	t0 <- rep(0,n)
	ht0 <- rep(NA,n)
	Fhatt0 <- c(0,cumsum(w[1:(n-1)]))
	Yhatt0 <- c(0,cumsum(w[1:(n-1)]*x[1:(n-1)]))
	t0 <- qgamma(Fhatt0,alpha,1)
	Fthetat0 <- pgamma(t0,alpha,1)
	ht0 <- t0*(Fhatt0 - Fthetat0) +
		mu - alpha - Yhatt0 + alpha*pgamma(t0,alpha+1,1)
	dx <- x[2:n] - x[1:(n-1)]
	dt0x <- pmin(t0[2:n] - x[1:(n-1)], x[2:n] - t0[2:n])
	JJ <- (2:n)[dt0x > mindiff*dx]
	return(list(tnew=t0[c(1,JJ)],htnew=ht0[c(1,JJ)]))
}


LocalDirDeriv2.2B <- function(t,x,w,tau,theta,alpha=1)
# This procedure computes the directional derivatives
#    h(t[i]) := DL((tau,theta),V_{t[i],(tau,theta)})
# and the right-sided derivatives
#    h'(t[i]+)
# for 1 <= i <= length(t).
# It is only used to illustrate and check the methods.
# The estimation procedure uses only Optimality2.2B() is used.
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
	
	ib <- 1
	if (tau[1] > 0){
		while (z[ib+1] <= tau[1]){
			ib <- ib+1
		}
		Phatz[1:ib] <- cumsum(wz[1:ib])
		Pthetaz[1:ib] <- exp(theta[1])*pgamma(z[1:ib],alpha,1)
		Yhatz[1:ib] <- cumsum(wz[1:ib]*z[1:ib])
		h1t[1:ib] <- Phatz[1:ib] - Pthetaz[1:ib]
		h00 <- -tau[1]*(Phatz[ib] - Pthetaz[ib]) +
			Yhatz[ib] -
			exp(theta[1])*alpha*pgamma(tau[1],alpha+1,1)
		h0t[1:ib] <- h00 + z[1:ib]*h1t[1:ib] -
			Yhatz[1:ib] +
			exp(theta[1])*alpha*pgamma(z[1:ib],alpha+1,1)
	}
	
	if (m > 1){
		for (j in 1:(m-1)){
			ia <- ib
			while (z[ib+1] <= tau[j+1]){
				ib <- ib+1
			}
				# Now, z[ia] = tau[j] and z[ib] = tau[j+1]
			dj <- tau[j+1] - tau[j]
			Hj <- sum(wz[ia:ib]*(tau[j+1] - z[ia:ib]))/dj -
				J10.2B(theta[j],theta[j+1],tau[j],tau[j+1],alpha)
			thetaz[ia:ib] <- ((tau[j+1] - z[ia:ib])*theta[j] +
					(z[ia:ib] - tau[j])*theta[j+1])/dj
			Pthetaz[ia:ib] <- J00.2B(theta[j],thetaz[ia:ib],
				tau[j],z[ia:ib],alpha)
			Phatz[ia:ib] <- cumsum(wz[ia:ib])
			Yhatz[ia:ib] <- cumsum(wz[ia:ib]*(z[ia:ib] - tau[j]))
			h1t[ia:ib] <- Phatz[ia:ib] - Pthetaz[ia:ib] - Hj
			h0t[ia:ib] <- -Yhatz[ia:ib] +
				(z[ia:ib] - tau[j])*(h1t[ia:ib] +
					J01.2B(theta[j],thetaz[ia:ib],
						tau[j],z[ia:ib],alpha))
      }
	}
	ia <- ib
	ib <- nz
	thetaz[ia:ib] <- theta[m] + (z[ia:ib] - tau[m])*theta[m+1]
	Phatz[(ib-1):ia] <- cumsum(wz[ib:(ia+1)])
	Pthetaz[ia:ib] <- K0.2B(thetaz[ia:ib],theta[m+1],
		z[ia:ib],alpha)
	Yhatz[ia:ib] <- cumsum(wz[ia:ib]*(z[ia:ib] - tau[m]))
	h1t[ia:ib] <- Pthetaz[ia:ib] - Phatz[ia:ib]
	h0t[ia:ib] <- -Yhatz[ia:ib] +
		(z[ia:ib] - tau[m])*
			(h1t[ia:ib] + J01.2B(theta[m],thetaz[ia:ib],
				tau[m],z[ia:ib],alpha))
	
	return(cbind("h(t)"=h0t[z.is.t],"h'(t+)"=h1t[z.is.t]))
}

Optimality2.2B <- function(x,w,tau,theta,alpha=1,
	mindiff=10^(-5))
# Determine knots t0 with locally maximal directional
# derivatives DL(theta,V_{t0,theta}) in case of a
# non-affine current fit given by (tau,theta).
# In case t[1] > 0, t0 = 0 has to be considered.
{
	t <- c(x,tau)
	wt <- c(w,rep(0,length(tau)))
	ord <- order(t)
	t <- t[ord]
	wt <- wt[ord]
	if (t[1] > 0){
		t <- c(0,t)
		wt <- c(0,wt)
	}
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

	# Search for potential new knots to the left of tau[1] > 0: 
	ib <- 1
	if (tau[1] > 0){
		while (t[ib+1] <= tau[1]){
			ib <- ib+1
		}
			# Now, we consider potential new knots at 0 or
			# between t[1] = x[1] and t[ib] = tau[1].
		Phatt[1:ib] <- cumsum(wt[1:ib])
			# Phatt[i] = empirical prob. of [0,t[i]]
		Pthetat[1:ib] <- exp(theta[1])*pgamma(t[1:ib],alpha,1)
			# Pthetat[i] = estimated prob. of [0,t[i]]
		Yhatt[1:ib] <- cumsum(wt[1:ib]*t[1:ib])
			# Yhatt[i] = empirical integral of x
			# over [0,tau[1]]
		tnew[1] <- 0
		htnew[1] <- -tau[1]*(Phatt[ib] - Pthetat[ib]) +
			Yhatt[ib] -
			exp(theta[1])*alpha*pgamma(tau[1],alpha+1,1)
		h1tp[2:(ib-1)] <- Phatt[2:(ib-1)] - Pthetat[2:(ib-1)]
		h1tm[2:(ib-1)] <- Phatt[2:(ib-1)] - Pthetat[3:ib]
		JJ <- (2:(ib-1))[h1tp[2:(ib-1)] > 0 & h1tm[2:(ib-1)] < 0]
			# Which intervals (t[i],t[i+1]) contain
			# a local maximum of h_theta?
		tnew[JJ] <- InverseF0Theta.2B(theta[1],Phatt[JJ],alpha)
		# Eliminate spurious local maxima:
		JJ1 <- JJ + 1
		dtJJ <- t[JJ1] - t[JJ]
		dtnewtJJ <- pmin(tnew[JJ] - t[JJ], t[JJ1] - tnew[JJ])
		failures <- (dtnewtJJ < mindiff*dtJJ)
		if (any(failures)){
			tnew[JJ[failures]] <- NA
			JJ <- JJ[!failures]
		}
		# Evaluate h at tnew[JJ]:
		htnew[JJ] <- htnew[1] - Yhatt[JJ] +
			exp(theta[1])*alpha*pgamma(tnew[JJ],alpha+1,1)
	}
	if (m > 1){
		for (j in 1:(m-1)){
			# Search for potential new knots between
			# tau[j] and tau[j+1]:
			ia <- ib
			while (t[ib+1] <= tau[j+1]){
				ib <- ib+1
			}
				# Now, t[ia] = tau[j] and t[ib] = tau[j+1].
			Phatt[ia:(ib-1)] <- cumsum(wt[ia:(ib-1)])
				# Phatt[i] = empirical prob. of (tau[j],t[i]]
			dj <- tau[j+1] - tau[j]
			Hj <- sum(wt[ia:ib]*(tau[j+1] - t[ia:ib]))/dj -
				J10.2B(theta[j],theta[j+1],
					tau[j],tau[j+1],alpha)
				# Hj = empirical minus estimated integral of
				# j_{10}(x; tau[j],tau[j+1]) over (tau[j],tau[j+1]]
			thetat[ia:ib] <- ((tau[j+1] - t[ia:ib])*theta[j] +
					(t[ia:ib] - tau[j])*theta[j+1])/dj
				# thetat[i] = theta(t[i])
			Pthetat[ia:ib] <- J00.2B(theta[j],thetat[ia:ib],
				tau[j],t[ia:ib],alpha)
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
			theta1j <- (theta[j+1] - theta[j])/dj
			tnew[JJ] <- InverseJ00.2B(theta[j],theta1j,
				tau[j],Phatt[JJ]-Hj,alpha)
			# Eliminate spurious local maxima:
			JJ1 <- JJ + 1
			dtJJ <- t[JJ1] - t[JJ]
			dtnewtJJ <- pmin(tnew[JJ] - t[JJ],
				t[JJ1] - tnew[JJ])
			failures <- (dtnewtJJ < mindiff*dtJJ)
			if (any(failures)){
				tnew[JJ[failures]] <- NA
				JJ <- JJ[!failures]
			}
			# Evaluate h at tnew[JJ]:
			thetastar[JJ] <- ((tau[j+1] - tnew[JJ])*theta[j] +
				(tnew[JJ] - tau[j])*theta[j+1])/dj
			htnew[JJ] <- -Yhatt[JJ] + (tnew[JJ] - tau[j])*
				J01.2B(theta[j],thetastar[JJ],
					tau[j],tnew[JJ],alpha)
		}
	}
	
	# Search for potential new knots to the right of tau[m]: 
	ia <- ib
	ib <- nt
		# Now, t[ia] = tau[m] and t[ib] = x[n].
	Phatt[(ib-1):ia] <- cumsum(wt[ib:(ia+1)])
		# Phatt[i] = empirical prob. of (t[i],Inf)
	thetat[ia:ib] <- theta[m] + (t[ia:ib] - tau[m])*theta[m+1]
		# thetat[i] = theta(t[i])
	Pthetat[ia:ib] <- K0.2B(thetat[ia:ib],theta[m+1],
		t[ia:ib],alpha)
		# Pthetat[i] = estimated prob. of (t[i],Inf)
	Yhatt[ia:ib] <- cumsum(wt[ia:ib]*(t[ia:ib] - tau[m]))
		# Yhatt[i] = empirical integral of (x - tau[m])
		# over (tau[m],t[i]]
	h1tp[ia:(ib-1)] <- Pthetat[ia:(ib-1)] - Phatt[ia:(ib-1)]
	h1tm[ia:(ib-1)] <- Pthetat[(ia+1):ib] - Phatt[ia:(ib-1)]
	JJ <- ia-1 + which(h1tp[ia:(ib-1)] > 0 & h1tm[ia:(ib-1)] < 0)
		# Which intervals (t[i],t[i+1]) contain
		# a local maximum of h_theta?
	tnew[JJ] <- InverseK0.2B(theta[m],theta[m+1],
		tau[m],Phatt[JJ],alpha)
	# Eliminate spurious local maxima:
	JJ1 <- JJ + 1
	dtJJ <- t[JJ1] - t[JJ]
	dtnewtJJ <- pmin(tnew[JJ] - t[JJ], t[JJ1] - tnew[JJ])
	failures <- (dtnewtJJ < mindiff*dtJJ)
	if (any(failures)){
		tnew[JJ[failures]] <- NA
		JJ <- JJ[!failures]
	}
	# Evaluate h at tnew[JJ]:
	thetastar[JJ] <- theta[m] + (tnew[JJ] - tau[m])*theta[m+1]
	htnew[JJ] <- -Yhatt[JJ] + (tnew[JJ] - tau[m])*
			J01.2B(theta[m],thetastar[JJ],tau[m],tnew[JJ],alpha)

	return(list(tnew=tnew[!is.na(tnew)],htnew=htnew[!is.na(tnew)]))
}


### Auxiliary programs 3 ----
### 2nd step size correction

LocalConvexity.2B <- function(tau,theta)
# Changes of slope,
# each slope and change of slope should be > 0 and
# if theta defines a strictly convex and increasing function on {tau[i] : ...}.
{
  m <- length(tau)
  if (m >= 2){
    dtau <- tau[2:m] - tau[1:(m-1)]
    dtheta <- theta[2:m] - theta[1:(m-1)]
    dtdt <- c(0,dtheta/dtau,theta[m+1])
  }
  else{
    dtdt <- c(0,theta[2])
  }
  chofslope <- dtdt[2:(m+1)] - dtdt[1:m]
  return(chofslope)
}

LocalStepSize2.2B <- function(tau,theta,theta.new,
                              isnewknot=rep(FALSE, length(tau)))
# Replace theta with
#    (1 - t0)*theta + t0*theta.new ,
# where t0 is the largest number in [0,1] such that
# the new theta defines still a convex function.
# Return new tuples tau and theta, possibly shorter
# than the original ones.
{
	m <- length(tau)
	chofsl.new <- LocalConvexity.2B(tau,theta.new)
	# 1st case
	if (m == 1 && chofsl.new <= 0) {
		theta <- c(0,0)
			# Constant function equal to 0 with knot at tau[1]
		return(list(tau=tau,theta=theta))
	}
	# 2nd easy (in fact, ideal) case:
	if (min(chofsl.new) > 0) {
		return(list(tau=tau,theta=theta.new))
	}
	# All other cases with min(chofsl.new) <= 0 lead to
	# removal of one knot / activation of one constraint:
	#	
	# (1) If at least one of the new knots is responsible for a
	# negative change of slope in theta.new, this means that the
	# constraint at the worst knot has to be activated, resulting
	# in the deletion of that knot. Then theta.new is 
	# discarded and will have to be recomputed with a Newton step.
	if (any(isnewknot) && min(chofsl.new[isnewknot]) <= 0) {
		tmp  <- which(chofsl.new <= 0 & isnewknot)
		j0 <- which.min(chofsl.new[tmp])
		j0 <- tmp[j0]
		theta <- theta[-j0]
		tau <- tau[-j0]
		isnewknot <- isnewknot[-j0]
		return(list(tau=tau,theta=theta,isnewknot=isnewknot))
	}
	# (2) At this stage, none of the newly added knots is
	# responsible for a negative change of slope in theta.new.
	# Nevertheless, at least one of the old constraints has to
	# be activated and the corresponding knot has to be
	# deleted. We therefore look for the largest t0 such that
	# (1 - t0)*theta + t0*theta.new is feasible.
	chofsl <- LocalConvexity.2B(tau,theta)
	tmp  <- which(chofsl.new <= 0)
	t0 <- chofsl[tmp]/(chofsl[tmp] - chofsl.new[tmp])
	j0 <- which.min(t0)
	t0 <- t0[j0]
	j0 <- tmp[j0]
	theta <- (1 - t0)*theta + t0*theta.new
	tau <- tau[-j0]
	isnewknot <- isnewknot[-j0]
	theta <- theta[-j0]
	m <- m-1
	# Now make sure that phi is even strictly concave at each knot:
	chofsl <- LocalConvexity.2B(tau,theta)
	while (m > 1 && min(chofsl) <= 0){
		j0 <- which.min(chofsl)
		tau <- tau[-j0]
		isnewknot <- isnewknot[-j0]
		theta <- theta[-j0]
		m <- m-1
		chofsl <- LocalConvexity.2B(tau,theta)
	}
	if (m == 1 && chofsl <= 0){
		theta <- c(0,0)
		# Constant function equal to 0 with knot at tau[1]
	}
	return(list(tau=tau,theta=theta,isnewknot=isnewknot))
}

### Auxiliary programs 2 ----
### Normalization, log-likelihood plus derivatives
### and Newton step with 1st step size correction

LocalNormalize.2B <- function(tau,theta,alpha=1)
# Normalize theta such that it defines
# a log-probability density
{
	m <- length(tau)
	integral <- exp(theta[1])*pgamma(tau[1], alpha, 1) +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		integral <- integral +
			sum(J00.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha))
	}
	theta.new <- theta
	theta.new[1:m] <- theta.new[1:m] - log(integral)
	return(theta.new)
}

LocalLL.2B <- function(wt,tau,theta,alpha=1)
# Log-likelihood
{
	m <- length(tau)
	integral <- exp(theta[1])*pgamma(tau[1], alpha, 1) +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		integral <- integral +
			sum(J00.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha))
	}
	LL <- sum(wt*theta) + 1 - integral
	return(LL)
}

LocalLL2.2B <- function(wt,tau,theta,alpha=1)
# Log-likelihood plus its gradient vector
# and minus Hessian matrix
{
	m <- length(tau)
	integral <- exp(theta[1])*pgamma(tau[1], alpha, 1) +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		integral <- integral +
			sum(J00.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha))
	}
	LL <- sum(wt*theta) + 1 - integral
	GLL <- wt
	GLL[1]   <- GLL[1] - exp(theta[1])*pgamma(tau[1], alpha, 1)
	GLL[m]   <- GLL[m] -   K0.2B(theta[m],theta[m+1],tau[m],alpha)
	GLL[m+1] <- GLL[m+1] - K1.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		GLL[1:(m-1)] <- GLL[1:(m-1)] -
			J10.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		GLL[2:m] <- GLL[2:m] -
			J01.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
	}
	MHLL <- matrix(0,m+1,m+1)
	MHLL[1,1] <- exp(theta[1])*pgamma(tau[1], alpha, 1)
	MHLL[m,m] <- MHLL[m,m] +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	MHLL[m,m+1] <- K1.2B(theta[m],theta[m+1],tau[m],alpha)
	MHLL[m+1,m] <- MHLL[m,m+1]
	MHLL[m+1,m+1] <- K2.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		A20 <- cbind(1:(m-1),1:(m-1))
		MHLL[A20] <- MHLL[A20] +
			J20.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		A02 <- cbind(2:m,2:m)
		MHLL[A02] <- MHLL[A02] +
			J02.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		tmp11 <- J11.2B(theta[1:(m-1)],theta[2:m],
			tau[1:(m-1)],tau[2:m],alpha)
		MHLL[cbind(1:(m-1),2:m)] <- tmp11
		MHLL[cbind(2:m,1:(m-1))] <- tmp11
	}
	return(list(LL=LL,GLL=GLL,MHLL=MHLL))
}

LocalNewton.2B <- function(wt,tau,theta,alpha=1,delta0=10^(-11))
# Newton step with (1st) step size correction
{
	m <- length(tau)
	integral <- exp(theta[1])*pgamma(tau[1], alpha, 1) +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		integral <- integral +
			sum(J00.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha))
	}
	LL <- sum(wt*theta) + 1 - integral
	GLL <- wt
	GLL[1]   <- GLL[1] - exp(theta[1])*pgamma(tau[1], alpha, 1)
	GLL[m]   <- GLL[m] -   K0.2B(theta[m],theta[m+1],tau[m],alpha)
	GLL[m+1] <- GLL[m+1] - K1.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		GLL[1:(m-1)] <- GLL[1:(m-1)] -
			J10.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		GLL[2:m] <- GLL[2:m] -
			J01.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
	}
	MHLL <- matrix(0,m+1,m+1)
	MHLL[1,1] <- exp(theta[1])*pgamma(tau[1], alpha, 1)
	MHLL[m,m] <- MHLL[m,m] +
		K0.2B(theta[m],theta[m+1],tau[m],alpha)
	MHLL[m,m+1] <- K1.2B(theta[m],theta[m+1],tau[m],alpha)
	MHLL[m+1,m] <- MHLL[m,m+1]
	MHLL[m+1,m+1] <- K2.2B(theta[m],theta[m+1],tau[m],alpha)
	if (m > 1) {
		A20 <- cbind(1:(m-1),1:(m-1))
		MHLL[A20] <- MHLL[A20] +
			J20.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		A02 <- cbind(2:m,2:m)
		MHLL[A02] <- MHLL[A02] +
			J02.2B(theta[1:(m-1)],theta[2:m],
				tau[1:(m-1)],tau[2:m],alpha)
		tmp11 <- J11.2B(theta[1:(m-1)],theta[2:m],
			tau[1:(m-1)],tau[2:m],alpha)
		MHLL[cbind(1:(m-1),2:m)] <- tmp11
		MHLL[cbind(2:m,1:(m-1))] <- tmp11
	}
	dtheta <- qr.solve(MHLL,GLL)
	dirderiv <- sum(dtheta * GLL)
	theta.new <- theta + dtheta
	LL.new <- LocalLL.2B(wt,tau,theta.new,alpha)
	while (LL.new < LL + dirderiv/3 && dirderiv > delta0) {
		theta.new <- (theta + theta.new)/2
		dirderiv <- dirderiv/2
		LL.new <- LocalLL.2B(wt,tau,theta.new,alpha)
	}
	return(list(theta.new=theta.new,dirderiv=dirderiv,
		LL=LL,LL.new=LL.new))
}

### Auxiliary programs 1 ----
### Linear splines and
### Basic functions J(.,.), K(.,.) plus 1st and 2nd
### partial derivatives

LocalLinearSplines1.2B <- function(tau,x=tau)
# For a vector tau with m >= 1 strictly increasing
# components, this function computes a design matrix Dx
# for the set of piecewise linear functions f on the real
# line with kinks only on {tau[j] : 1 <= j <= m}, evaluated
# at x. The basis corresponds to the parameter vector theta
# with components
#    theta[j]   = f(tau[j]), 1 <= j <= m,
#    theta[m+1] = f'(tau[m] +) .
{
	m <- length(tau)
	n <- length(x)
	Dx <- matrix(0,n,m+1)
	tmp <- (x <= tau[1])
	Dx[tmp,1] <- 1
	if (m > 1){
		dtau <- tau[2:m] - tau[1:(m-1)]
		for (j in 1:(m-1)){
			tmp <- (x > tau[j] & x <= tau[j+1])
			Dx[tmp,j  ] <- (tau[j+1] - x[tmp])/dtau[j]
			Dx[tmp,j+1] <- (x[tmp] - tau[j])/dtau[j]
		}
	}
	tmp <- (x > tau[m])
	Dx[tmp,m] <- 1
	Dx[tmp,m+1] <- x[tmp] - tau[m]
	return(Dx)
}

LocalLinearSplines2.2B <- function(tau,
	x=tau,w=rep(1/length(x),length(x)))
# For linear splines with given knots
#   tau[1] < tau[2] < ... < tau[m]
# and arbitrary vectors x and w of equal length with w >= 0,
# this procedure returns a weight vector wt, such that for
# the basis functions f_1, f_2, ..., f_m, f_{m+1}, f_{m+2}
# described/computed in LocalLinearSplines1(),
#    wt[j] = sum_{i=1}^n w[i]*f_j(x[i]) .
{
	Dx <- LocalLinearSplines1.2B(tau,x)
	wt <- colSums(Dx * w)
	return(wt)
}


G.2B <- function(a,b,s=1) 
# G(a,b;s) = F0(b,s) - F0(a,s).
{
	pgamma(b, s, 1) - pgamma(a, s, 1)
}

J00.2B <- function(x,y,a=0,b=1,alpha=1,epsilon=1e-4) 
# J00.2B(x,y,a,b)
#    =  integral_a^b
#          exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J00xy <- rep(0,length(x))
	II <- (yt < 1)
	J00xy[II] <- exp(xt[II] +
			log(G.2B(at[II], bt[II], alpha))) *
		(1 - yt[II])^(- alpha)
	II <- (yt >= 1)
		# In this case, use an upper bound
		# (exact value in case of yt = 1)
	J00xy[II] <- exp(xt[II] + (yt[II] - 1)*b[II]) *
		(b[II]^alpha - a[II]^alpha) / gamma(alpha + 1)
	return(J00xy)
  # # An alternative option for numerical stability.
  # # If it is used, the last "return(J00xy)" should be deleted.
  # II <- (yt <= 1-epsilon | yt >= 1)
  # if (all(II)) {
  #   return(J00xy)
  # }else{
  #   II <- !II
  #   JJ <- (1:length(xt))[II]
  #   c0 <- (b[II]^ alpha      - a[II]^ alpha)     / alpha
  #   c1 <- (b[II]^(alpha + 1) - a[II]^(alpha + 1))/(alpha + 1)
  #   J00.L <- exp(xt[II] + (yt[II] - 1)*c1/c0) *
  #     c0/gamma(alpha)
  #   J00.U <- exp(xt[II] + (yt[II] - 1)*a[II]) *
  #     c0/gamma(alpha)
  #   KK <- which(J00xy[II] < J00.L | J00.U < J00xy[II])
  #   JJ <- JJ[KK]
  #   if (length(JJ) > 0){
  #     J00xy[JJ] <- (J00.L[KK] + J00.U[KK])/2
  #   }
  #   return(J00xy)
  # }
}

J10.2B <- function(x,y,a=0,b=1,alpha=1) 
# J10.2B(x,y,a,b)
#    =  integral_a^b
#          (1-v(z))*exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J10xy <- exp(xt) * 
		(1 - yt)^(- alpha - 1) * 
			(bt * G.2B(at, bt, alpha) -
				alpha*G.2B(at, bt, alpha + 1)) /
			(b - a + (a == b))
	return(J10xy)
}

J01.2B <- function(x,y,a=0,b=1,alpha=1) 
# J01.2B(x,y,a,b)
#    =  integral_a^b
#          v(z)*exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J01xy <- exp(xt) * 
		(1 - yt)^(- alpha - 1) * 
			(-at * G.2B(at, bt, alpha) +
				alpha*G.2B(at, bt, alpha + 1)) /
			(b - a + (a == b))
	return(J01xy)
}

J20.2B <- function(x,y,a=0,b=1,alpha=1) 
# J20.2B(x,y,a,b,alpha)
#    =  integral_a^b
#          (1-v(z))^2*exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).  
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J20xy <- exp(xt) * 
		(1 - yt)^(-alpha-2) * 
			(             bt^2 * G.2B(at, bt, alpha  ) 
			      - 2*alpha*bt * G.2B(at, bt, alpha+1)
			 + alpha*(alpha+1) * G.2B(at, bt, alpha+2)) / 
			(b - a + (a == b))^2
	return(J20xy)
}

J11.2B <- function(x,y,a=0,b=1,alpha=1) 
# J11.2B(x,y,a,b,alpha)
#    =  integral_a^b
#          (1-v(z))*v(z)*exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).  
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J11xy <- exp(xt) * 
		(1 - yt)^(- alpha - 2) * 
			(          -at*bt  * G.2B(at, bt, alpha  )
			   + alpha*(at+bt) * G.2B(at, bt, alpha+1)
			 - alpha*(alpha+1) * G.2B(at, bt, alpha+2)) /
			(b - a + (a == b))^2
	return(J11xy)
}

J02.2B <- function(x,y,a=0,b=1,alpha=1) 
# J02.2B(x,y,a,b,alpha)
#    =  integral_a^b
#          v(z)^2*exp((1-v(z))x + v(z)y)*pgamma(z, alpha, 1) dz
# with v(z) = (z - a)/(b - a).  
{
	xt <- (b*x - a*y)/(b - a + (a == b))
	yt <-   (y -   x)/(b - a + (a == b))
	bt <- (1 - yt)*b
	at <- (1 - yt)*a
	J02xy <- exp(xt) * 
		(1 - yt)^(- alpha - 2) * 
			(             at^2 * G.2B(at, bt, alpha  ) 
			      - 2*alpha*at * G.2B(at, bt, alpha+1)
			 + alpha*(alpha+1) * G.2B(at, bt, alpha+2)) /
			(b - a + (a == b))^2
	return(J02xy)
}

K0.2B <- function(x,y,a=0,alpha=1)
# K0.2B(x,y,a,alpha)
# = integral_a^Inf exp(x + y*(z-a))*pgamma(z, alpha, 1) dz
{
	if(y >= 1) {
		return(Inf)
	}
	return(exp(x - a*y)*
			(1-pgamma((1-y)*a, alpha, 1))/(1-y)^alpha)
}

K1.2B <- function(x,y,a=0,alpha=1)
# K1.2B(x,y,a,alpha) = dK0.2B(x,y,a,alpha)/dy
{
	at <- (1 - y)*a
	return(exp(x - a*y)*
		(   - at*(1 - pgamma(at, alpha,   1))
		 + alpha*(1 - pgamma(at, alpha+1, 1)))/
		(1 - y)^(alpha + 1))
}

K2.2B <- function(x,y,a=0,alpha=1)
# K2.2B(x,y,a,alpha) = d^2K0.2B(x,y,a)/dy^2
{
	at <- (1 - y)*a
	return(exp(x - a*y)*
		(             at^2*(1 - pgamma(at, alpha,   1))
		      - 2*alpha*at*(1 - pgamma(at, alpha+1, 1))
		 + alpha*(alpha+1)*(1 - pgamma(at, alpha+2, 1)))/
		(1 - y)^(alpha + 2))
}

InverseF0Theta.2B <- function(theta0, c, alpha=1)
# Determines the unique solution tau of
#    exp(theta0)*pgamma(tau,alpha,1) = c.
{
	d <- pmin(pmax(c*exp(-theta0),0),1)
	return(qgamma(d, alpha, 1))
}

InverseJ00.2B <- function(theta0, theta1, tau0, c, alpha=1)
# Determines the unique solution tau of
#    J00(theta0, theta0 + theta1*(tau-tau0),
#        tau0, tau, alpha) = c.
{ 
	if (theta1 < 0 || theta1 >= 1){
		return(NULL)
	}
	d <- c*(1-theta1)^alpha*exp(theta1*tau0 - theta0) +
		pgamma((1 - theta1)*tau0, alpha, 1)
	d <- pmin(pmax(d,0),1)
	return(qgamma(d, alpha, 1)/(1 - theta1))
}

InverseK0.2B <- function(theta0, theta1, tau0, c, alpha=1)
# Determines the unique solution tau of
#    K0.2B(theta0 + theta1*(tau-tau0),theta1,tau,alpha) = c.
{
	if (theta1 < 0 || theta1 >= 1){
		return(NULL)
	}
	d <- c*(1-theta1)^alpha*exp(theta1*tau0 - theta0)
	d <- pmin(pmax(d,0),1)
	return(qgamma(1-d, alpha, 1)/(1 - theta1))
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

### Auxiliary programs 0 ----
### Data preparation,
### data simulation and
### plotting functions

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

Evaluate.2B <- function(t, tau, theta) 
# For plotting purposes only.
# For a vector t of nt = length(t) non-negative elements, 
# this function returns the value of theta(t). An active 
# knot at 0 is artificially added in case tau[1] > 0.
{
	if (tau[1] > 0) {
		tau <- c(0,tau)
		theta <- c(theta[1], theta)
	}
	m <- length(tau)
	dtau <- tau[2:m] - tau[1:(m-1)]
	dtheta <- theta[2:m] - theta[1:(m-1)]
	dtdt <- c(dtheta/dtau,theta[m+1])
	tau <- c(tau, Inf)
	nt <- length(t)
	res <- rep(0, nt)
	for(i in 1:nt) {
		j <- which(tau[1:m] <= t[i] & t[i] < tau[2:(m+1)])
		res[i] <- theta[j] + (t[i]-tau[j])*dtdt[j]
	}
	return(res)
}

Simulate.2B <- function(n,tau,theta,alpha=1,beta=1,sim.num=FALSE)
# The function returns simulated data 
# (X_1, X_2, ..., X_N) in R^N with density
#      f_theta(x) = f_o(x) * exp(theta(x)) / M(theta)
# where f_o(x) is the density of Gamma(alpha,beta) and
# theta is a piecewise linear convex and increasing
# function, with values theta at knots tau. The last
# element of theta is the last slope of theta(x) and it has
# to be smaller than beta.
# If sim.num is TRUE, returns the number of needed
# simulations and the vectors tau and theta.
{
	# Add a knot at 0 if needed
	if (tau[1] > 0) {
		tau <- c(0,tau)
		theta <- c(theta[1], theta)
	}
	m <- length(tau)
	if (is.infinite(tau[m])) {
		tau <- tau[1:(m-1)]
		m <- m - 1
	}
	# Knots, slopes, values of theta at the knots:
	if (m == 1) {
		theta_prime <- theta[m+1]
	} else {
		theta_prime <- rep(0, m)
		theta_prime[1:(m-1)] <- (theta[2:m]-theta[1:(m-1)])/
			(tau[2:m]-tau[1:(m-1)])
		theta_prime[m] <- theta[m+1]
	}
	tau <- c(tau, Inf)
	theta <- theta[1:m]
	beta_b <- beta - theta_prime[m]
	
	# g function, g(x) in [0,1] for all x
	g <- function(x){
		k_x <- (1:m)[(tau[1:m] <= x & x < tau[2:(m+1)])]
		res <- (theta_prime[k_x] - theta_prime[m])*(x - tau[k_x]) -
			theta_prime[m]*tau[k_x] + theta[k_x] - theta[1]
		return(exp(res))
	}
	
	# Acceptance rejection:
	X <- rep(0,n)
	i <- 1
	nn <- 0
	while(i <= n){
		Y <- rgamma(1, shape = alpha, rate = beta_b)
		U <- runif(1)
		if (U <= g(Y)) {
			X[i] <- Y
			i <- i+1
		}
		nn <- nn + 1
	}
	if (sim.num == TRUE) {
		return(list(X = X,sim.num = nn,
			theta = c(theta, theta_prime[m]),tau = tau))
	} else {
		return(X)
	}
}

