getPost <-
function(r, kmin=2, kmax=5, S=200){
	r <- r[!is.na(r)] ## remove missing values
	n <- length(r)
	loglik <- c()
	bic <- c()
	kpost <- list()
	for(k in kmin:kmax){
		nn <- rep(length(r)/k,k)
		alpha <- rep(1,k)
		mus <- kmeans(r, centers=k, nstart=15)$centers
		mus <- sort(mus)
		print(k); print(mus)

		post <- gibbs(r=r,nn=nn, K=k, alpha=alpha, mu0=mus, S=S)
		loglik[k-kmin+1] <- loglik.normmix(r, post)
		bic[k-kmin+1] <- -2*loglik[k-kmin+1] + (3*k-1)*log(n)
		#save posterior
		kpost[[k-kmin+1]] <- post
	}
	## print table showing log-likelihoods and BIC
	results <- data.frame("K"=seq(kmin,kmax), "log-likelihood"=loglik, "BIC"=bic)
	print(results)
	## Return posterior of model with largest BIC
	return(kpost)
}
