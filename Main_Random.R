

# ----------------------------------------------------------------------
# R = 1
# ----------------------------------------------------------------------
lam.vec <- 10^seq(2, -2, length=100)
delta1.init <- list()
beta1.init <- list()
ours1.val <- rep(-Inf,  length(lam.vec))
ours1.KL <- rep(Inf, length(lam.vec))

for(kk in 1:length(lam.vec)){

  	if (kk == 1) {
  		temp <- EM_SparseGroup(X, Y, R = 1, lambda = lam.vec[kk], 
  			delta.init = NULL, 
  			beta.init = NULL, 
  			max.iter.EM = 1e3, gamma0 = 5,
  			obj.tol = 1e-7)
  	} else {
  		temp <- EM_SparseGroup(X, Y, R = 1, lambda = lam.vec[kk], 
  			delta.init = delta1.init[[kk-1]], 
  			beta.init = beta1.init[[kk-1]], 
  			max.iter.EM = 1e3, gamma0 = 5,
  			obj.tol = 1e-7)
  	}
		ours1.val[kk] <- obs.data.loglik(Yval, Xval, temp$beta.list, temp$deltas, R = 1, lambda = 0, weights = rep(0, p))
		delta1.init[[kk]] <- temp$deltas
		beta1.init[[kk]] <- temp$beta.list


		cat("# ------------------ ", "\n")
		cat(kk, ":", ours1.val[kk], "\n")
		cat("# ------------------ ", "\n")
		# ---------------------------------------------------------
		# End if validation error has increased last five lambda
		# -------------------------------------------------------
		if (kk > 20) {
			if(sum(ours1.val[kk-(10:1)] > ours1.val[kk-(9:0)]) == 10){
				break
			}
		}
}

# ----------------------------------------------
# Get the best fitted model
# ----------------------------------------------
best.ind <- which(ours1.val == max(ours1.val))[1]
delta1 <- delta1.init[[best.ind]]
beta1 <- beta1.init[[best.ind]]
temp <- list("beta.list" = beta1, "deltas" = delta1)
ours1 <- predict.EM(Xnew, temp)
ours1.KL <- sum(ours1$probs*(log(ours1$probs)-log(pop$probs)))
ours1.HELL <- mean(apply(sqrt(ours1$probs) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)
Results <- append(Results, namedList(ours1.KL, ours1.HELL, delta1, beta1, ours1.val, delta1.init))

# ------------------------------
# Clear for memory 
# ------------------------------
beta1.init <- NULL
delta1.init <- NULL




# ----------------------------------------------------------------------
# R = 2
# ----------------------------------------------------------------------
lam.vec <- 10^seq(2, -2, length=100)
delta2.init <- list()
ours2.val <- rep(-Inf,  length(lam.vec))
ours2.KL <- rep(Inf, length(lam.vec))
beta2.init <- list()


for(kk in 1:length(lam.vec)){
	
	temp <- EM_SparseGroup(X, Y, R = 2, lambda = lam.vec[kk], 
		delta.init = NULL, 
		beta.init = NULL, 
		max.iter.EM = 1e3, gamma0 = 5, 
		obj.tol = 1e-7, epochs=1)

	ours2.val[kk] <- obs.data.loglik(Yval, Xval, temp$beta.list, temp$deltas, R=2, lambda=0, weights= rep(0, p))
	delta2.init[[kk]] <- temp$deltas
	beta2.init[[kk]] <- temp$beta.list

	cat("# ------------------ ", "\n")
	cat(kk, ":", ours2.val[kk], "; deltas = ", temp$deltas, "\n")
	cat("# ------------------ ", "\n")
	# ---------------------------------------------------------
	# End if validation error has increased last five lambda
	# -------------------------------------------------------
	if (kk > 10) {
	  if(sum(ours2.val[(kk-9):kk] < max(ours2.val)) == 10){
	    break
	  }
	}
}

# ----------------------------------------------
# Get the best fitted model
# ----------------------------------------------
best.ind <- which(ours2.val == max(ours2.val))[1]
delta2 <- delta2.init[[best.ind]]
beta2 <- beta2.init[[best.ind]]
temp <- list("beta.list" = beta2, "deltas" = delta2)
ours2 <- predict.EM(Xnew, temp)
ours2.KL <- sum(ours2$probs*(log(ours2$probs)-log(pop$probs)))
ours2.HELL <- mean(apply(sqrt(ours2$probs) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)
Results <- append(Results, namedList(ours2.KL, ours2.HELL, delta2, beta2, ours2.val, delta2.init))

# ------------------------------
# Clear for memory 
# ------------------------------
beta2.init <- NULL
delta2.init <- NULL

 




# ----------------------------------------------------------------------
# R = 3
# ----------------------------------------------------------------------
lam.vec <- 10^seq(2, -2, length=100)
delta3.init <- list()
beta3.init <- list()
ours3.val <- rep(-Inf,  length(lam.vec))
ours3.KL <- rep(Inf, length(lam.vec))


for(kk in 1:length(lam.vec)){

	temp <- EM_SparseGroup(X, Y, R = 3, lambda = lam.vec[kk], 
		delta.init = NULL, 
		beta.init = NULL, 
		max.iter.EM = 1e3, gamma0 = 5,
		obj.tol = 1e-7, epochs=1)

	ours3.val[kk] <- obs.data.loglik(Yval, Xval, temp$beta.list, temp$deltas, R = 3, lambda = 0, weights = rep(0, p))
	delta3.init[[kk]] <- temp$deltas
	beta3.init[[kk]] <- temp$beta.list


	cat("# ------------------ ", "\n")
	cat(kk, ":", ours3.val[kk], "\n")
	cat("# ------------------ ", "\n")
	# ---------------------------------------------------------
	# End if validation error has increased last five lambda
	# -------------------------------------------------------
	if (kk > 20) {
	  if(sum(ours3.val[(kk-9):kk] < max(ours3.val)) == 10){
	    break
	  }
	}
}

# ----------------------------------------------
# Get the best fitted model
# ----------------------------------------------
best.ind <- which(ours3.val == max(ours3.val))[1]
delta3 <- delta3.init[[best.ind]]
beta3 <- beta3.init[[best.ind]]
temp <- list("beta.list" = beta3, "deltas" = delta3)
ours3 <- predict.EM(Xnew, temp)
ours3.KL <- sum(ours3$probs*(log(ours3$probs)-log(pop$probs)))
ours3.HELL <- mean(apply(sqrt(ours3$probs) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)
Results <- append(Results, namedList(ours3.KL, ours3.HELL, delta3, beta3, ours3.val, delta3.init))

# ------------------------------
# Clear for memory 
# ------------------------------
beta3.init <- NULL
delta3.init <- NULL




# ----------------------------------------------------------------------
# R = 4
# ----------------------------------------------------------------------
lam.vec <- 10^seq(2, -2, length=100)
delta4.init <- list()
beta4.init <- list()
ours4.val <- rep(-Inf,  length(lam.vec))
ours4.KL <- rep(Inf, length(lam.vec))

for(kk in 1:length(lam.vec)){

	temp <- EM_SparseGroup(X, Y, R = 4, lambda = lam.vec[kk], 
		delta.init = NULL, 
		beta.init = NULL, 
		max.iter.EM = 1e3, gamma0 = 5,
		obj.tol = 1e-7, epochs=1)

	ours4.val[kk] <- obs.data.loglik(Yval, Xval, temp$beta.list, temp$deltas, R = 4, lambda = 0, weights = rep(0, p))
	delta4.init[[kk]] <- temp$deltas
	beta4.init[[kk]] <- temp$beta.list


	cat("# ------------------ ", "\n")
	cat(kk, ":", ours4.val[kk], "\n")
	cat("# ------------------ ", "\n")
	# ---------------------------------------------------------
	# End if validation error has increased last five lambda
	# -------------------------------------------------------
	if (kk > 20) {
	  if(sum(ours4.val[(kk-9):kk] < max(ours4.val)) == 10){
	    break
	  }
	}
}

# ----------------------------------------------
# Get the best fitted model
# ----------------------------------------------
best.ind <- which(ours4.val == max(ours4.val))[1]
delta4 <- delta4.init[[best.ind]]
beta4 <- beta4.init[[best.ind]]
temp <- list("beta.list" = beta4, "deltas" = delta4)
ours4 <- predict.EM(Xnew, temp)
ours4.KL <- sum(ours4$probs*(log(ours4$probs)-log(pop$probs)))
ours4.HELL <- mean(apply(sqrt(ours4$probs) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)
Results <- append(Results, namedList(ours4.KL, ours4.HELL, delta4, beta4, ours4.val, delta4.init))

# ------------------------------
# Clear for memory 
# ------------------------------
beta4.init <- NULL
delta4.init <- NULL





# ------------------------------------------------
# Separate group lasso
# ------------------------------------------------
test.preds <- list()
betaGroup.glmnet <- list()
for(kk in 1:M){

	t0 <- glmnet(y = Y[[kk]], x = X, family="multinomial", intercept=FALSE, type.multinomial = "grouped")
	preds <- predict(t0, Xval, type="response")
	errs <- rep(0, length(t0$lambda))
	for(ll in 1:length(t0$lambda)){
		errs[ll] <-  -sum(log(rowSums(Yval[[kk]]*preds[,,ll])))
	}
	best.ind <- which.min(errs)
	test.preds[[kk]] <- predict(t0, Xnew, type="response")[,,best.ind]
	betaGroup.glmnet[[kk]] <- Reduce(cbind, coef(t0, s = t0$lambda[best.ind]))

}


predArray <- array(0, dim = c(dim(Xnew)[1], Ydims))
for(hh in 1:dim(Xnew)[1]){
	predArray[hh,,,,] <- test.preds[[1]][hh,]%o%test.preds[[2]][hh,]%o%test.preds[[3]][hh,]%o%test.preds[[4]][hh,]
}

sepGroup.KL <- sum(predArray*(log(predArray)-log(pop$probs)))
sepGroup.HELL <- mean(apply(sqrt(predArray) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)

Results <- append(Results, namedList(betaGroup.glmnet, sepGroup.KL, sepGroup.HELL))


# ------------------------------------------------
# Separate lasso
# ------------------------------------------------
test.preds <- list()
betaUngroup.glmnet <- list()
for(kk in 1:M){

	t0 <- glmnet(y = Y[[kk]], x = X, family="multinomial", intercept=FALSE, type.multinomial = "ungrouped")
	preds <- predict(t0, Xval, type="response")
	errs <- rep(0, length(t0$lambda))
	for(ll in 1:length(t0$lambda)){
		errs[ll] <-  -sum(log(rowSums(Yval[[kk]]*preds[,,ll])))
	}
	best.ind <- which.min(errs)
	test.preds[[kk]] <- predict(t0, Xnew, type="response")[,,best.ind]
	betaUngroup.glmnet[[kk]] <- Reduce(cbind, coef(t0, s = t0$lambda[best.ind]))

}


predArray <- array(0, dim = c(dim(Xnew)[1], Ydims))
for(hh in 1:dim(Xnew)[1]){
	predArray[hh,,,,] <- test.preds[[1]][hh,]%o%test.preds[[2]][hh,]%o%test.preds[[3]][hh,]%o%test.preds[[4]][hh,]
}

sepUngroup.KL <- sum(predArray*(log(predArray)-log(pop$probs)))
sepUngroup.HELL <- mean(apply(sqrt(predArray) - sqrt(pop$probs), 1, function(x){sqrt(sum(x^2))}))/sqrt(2)

Results <- append(Results, namedList(betaUngroup.glmnet, sepUngroup.KL, sepUngroup.HELL))

saveRDS(Results, file = savename)
