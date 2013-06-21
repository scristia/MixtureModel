gibbs <-
function(r, nn,
		  alpha=rep(1,3), ## dirichlet prior
		  mu0=c(-4, -0.5, 0), ## theta ~ normal(mu0, tau20)
		  tau20=0.1,  ## variance of theta
		  nu0=1, ## number of prior observations for precision
		  kappa0=1, ## number of prior observations for mu0
		  sigma20=0.1, ## prec ~ gamma(nu0/2, nu0/2*sigma20)
		  K=3, ## assume 3 components (for now)
		  S=100
		  ){

	if(missing(nn)) stop("nn missing")
	if(sum(nn) != length(r)) stop("nn must sum to length(r)")
	z <- rep(NA, length(r))
	s2n <- tau2n <- mun <- numeric(K)
	p <- matrix(NA, sum(nn), K)
	nun <- nu0+nn
	
	## for sampling from constrained full conditionals
	a0 <- min(r)
	b0 <- max(r)

	if(missing(mu0)){
		## assume the mode is copy number 2
		mu0 <- numeric(K)
		rr <- r[!is.na(r)]
		d <- density(rr)
		mu0[3] <- d$x[which.max(d$y)]
		mu0[2] <- mu0[1] - 0.5
		mu0[1] <- mu0[3] - 2
	}

	s2 <- var(r, na.rm=TRUE)
	precs <- means <- matrix(NA, S, K)
	## simulate from prior
	rbar <- means[1, ] <- rnorm(K, mu0, tau20)
	Z <- matrix(NA, length(r), S-1)
	## just use marginal variance as guess of variance -- very diffuse
	precs[1,] <- 1/rep(s2, K)
	PI <- matrix(NA,S, K)
	theta <- c()
	for(s in 2:S){
		##
		## simulate pi from its multinomial posterior
		##
		pi <- rdirichlet(1, alpha+nn)
		##
		## update mu_n and tau_n^2
		##
		tau2n <- 1/(1/tau20 + nn*precs[s-1, ])
		nun <- nu0+nn
		s2n <- 1/nun * (nu0*sigma20 + (nn-1)*s2 + kappa0*nn/(kappa0+nn)*(rbar - mu0)^2)
		mun <- (1/tau20)/(1/tau20 + nn*precs[s-1, ])*mu0 + nn*precs[s-1, ]/(1/tau20 + nn*precs[s-1, ])*rbar
		##
		## simulate from full conditional for theta
		## samples from the constrained distributions for each theta
		## 
		for(k in 1:K){
			if(k==1){
				a <- a0
				b <- means[s-1, k+1]
				theta[k] <- constr.draw(mun[k], tau2n[k],a, b)
			}
			else if(k==K){
				a <- theta[k-1]
				b <- b0
				theta[k] <- constr.draw(mun[k], tau2n[k], a, b)
			}
			else{
				a <- theta[k-1]
				b <- means[s-1, k+1]
				theta[k] <- constr.draw(mun[k], tau2n[k], a, b)
			}
		}	
		## check: make sure contstrained sampling is working properly
		inc <- diff(theta) > 0
		cnt <- 0
		while(!all(inc)){
			theta <- rnorm(K, mun, sqrt(tau2n))
			inc <- diff(theta) > 0
			cnt <- cnt+1
			if(cnt > 20) stop("theta's not ordered")
		}

		##
		## simulate precision from full conditional
		##
		prec <- rgamma(K, nun/2, nun/2 * s2n)
		means[s, ] <- theta
		precs[s, ] <- prec

		## simulate latent variable
		d <- matrix(NA, nrow = length(r), ncol = K)
		for(i in 1:K) d[,i] <- pi[i]*dnorm(r, theta[i], sqrt(1/prec[i]))
		p <- d/apply(d, 1, sum)
	
		#z <- rMultinom(p,1) - 1 ## Requires Hmisc package
		
		u <- runif(length(r))
		tmp <- p[, 1]
		z[u < tmp] <- 0
		for(i in 2:K){
			z[tmp < u & u < tmp + p[, i]] <- i-1
			tmp <- tmp + p[, i]
		}	
		
		##
		## update [pi|data]
		##
		for(i in 1:K) nn[i] <- sum(z==(i-1)) 

		if(any(nn < 10 | is.na(nn))){
		browser()
		}
		##
		rbar <- sapply(split(r, z), mean, na.rm=TRUE)
		s2 <- sapply(split(r, z), var, na.rm=TRUE)
		## for identifiability
		if(any(diff(rbar) < 0)){
			i <- which(diff(rbar) < 0)
			rbar[i+1] <- rbar[i]+0.01
		}
		Z[, s-1] <- z
		PI[s, ] <- pi
	}
	list(P=PI, means=means, precs=precs, Z=Z)
}
