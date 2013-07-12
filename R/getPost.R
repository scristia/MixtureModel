getPost <-
function(r, kmin=2, kmax=5, delta=0.15, S=200, plot=F, burnin=1:100, main=""){
	r <- r[!is.na(r)] ## remove missing values
	n <- length(r)
	loglik <- c()
	bic <- c()
	kpost <- list()
	for(k in kmin:kmax){
#		nn <- rep(length(r)/k,k
#		mus <- kmeans(r, centers=k, nstart=15)$centers
#		mus <- sort(mus)
		inits <- inits(r,k)
		mu0 <- inits$mu0
		sigma20 <- inits$sigma20
		nn <- inits$nn
		print(nn)
		alpha <- rep(1,k)
		print(mu0)

		post <- gibbs(r=r,nn=nn, tau20=0.1, K=k, alpha=alpha, sigma20=sigma20, mu0=mu0, S=S, delta=delta)
		loglik[k-kmin+1] <- post$loglik
		bic[k-kmin+1] <- post$bic
#		loglik[k-kmin+1] <- loglik.normmix(r, post, K=k, burnin=burnin)
#		bic[k-kmin+1] <- -2*loglik[k-kmin+1] + (3*k-1)*log(n)
		#save posterior
		kpost[[k-kmin+1]] <- post
	}
	## print table showing log-likelihoods and BIC
	results <- data.frame("K"=seq(kmin,kmax), "log-likelihood"=loglik, "BIC"=bic)
	print(results)
	## Return posterior of model with largest BIC
	if(plot){
		plotPosts(r=r, posts=kpost, burnin=burnin, main=main)
#		dens.comp <- function(y, comp, p, mix.mean, mix.prec) {
#				p[comp]*dnorm(y, mean=mix.mean[comp], sd=1/sqrt(mix.prec[comp]))
#			}
#		y <- seq(min(r), max(r), len=1000)
#		## Sort by BIC
#		best.bic <- results[order(results[ ,3])[1:2], ]
#		subs <- best.bic$K - kmin + 1
#		## Get posterior estimates of two best fitting models
#		if(ncol(kpost[[subs[1]]][["P"]]) == 1){
#			p.1 <- mean(kpost[[subs[1]]][["P"]][-burnin, ])
#			mean.1<- mean(kpost[[subs[1]]][["means"]][-burnin, ])
#			prec.1 <- mean(kpost[[subs[1]]][["precs"]][-burnin, ])
#		}
#		else{
#			p.1 <- apply(kpost[[subs[1]]][["P"]][-burnin, ], 2, mean)
#			mean.1<- apply(kpost[[subs[1]]][["means"]][-burnin, ], 2, mean)
#			prec.1 <- apply(kpost[[subs[1]]][["precs"]][-burnin, ], 2, mean)
#		}
#			## second
#		if(ncol(kpost[[subs[2]]][["P"]]) == 1){
#			p.2 <- mean(kpost[[subs[2]]][["P"]][-burnin, ])
#			mean.2<- mean(kpost[[subs[2]]][["means"]][-burnin, ])
#			prec.2 <- mean(kpost[[subs[2]]][["precs"]][-burnin, ])
#		}
#		else{
#			p.2 <- apply(kpost[[subs[2]]][["P"]][-burnin, ], 2, mean)
#			mean.2<- apply(kpost[[subs[2]]][["means"]][-burnin, ], 2, mean)
#			prec.2 <- apply(kpost[[subs[2]]][["precs"]][-burnin, ], 2, mean)
#		}
#		d.1 <- rowSums(sapply(1:length(p.1), dens.comp, y=y, 
#				      p=p.1, mix.mean=mean.1, mix.prec=prec.1), na.rm=TRUE)
#		d.2 <- rowSums(sapply(1:length(p.2), dens.comp, y=y, 
#				      p=p.2, mix.mean=mean.2, mix.prec=prec.2), na.rm=TRUE)
#		## Draw histogram of data
#		hist(r, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main=main)
#		## Overlay densities in different colors
#		lines(y, d.1, col="darkgreen", lwd=2)
#		lines(y, d.2, col=rgb(24,167,181, max=255), lty=2, lwd=2)
#		## Write key with K and BIC
#		legend("topleft", legend=paste(paste("K = ", best.bic[,1], " BIC = ", round(best.bic[,3],2))), col=c("darkgreen", rgb(24,167,181,max=255)), lty=1:2, bty="n")
	}
#	if(!plot) return(kpost)
	return(kpost)
}
