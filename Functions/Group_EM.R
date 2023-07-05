# -------------------------------------
# EM (mixture of latent independents)
# -------------------------------------
# X: n x p matrix of predictors (without intercept)
# Y: list of length M -- each component is n x cm matrix with 1 in the column corresponding to observed categorical response
# delta.init: prob of mixture components
# beta.init: initial values for beta's
# ---------------------------------------

EM_SparseGroup <- function(X, 
                           Y, R, lambda, weights = NULL,
                           delta.init = NULL, 
                           beta.init = NULL, 
                           max.iter.EM = 1000, 
                           obj.tol = 1e-8, gamma0 = 1, 
                           epochs = 1) {
  
  # -----------------------
  # Preliminaries
  # -----------------------
  n <- dim(X)[1]
  p <- dim(X)[2]
  M <- length(Y)
  
  # ------------------------------
  # Set initial values for beta
  # ------------------------------
  if (is.null(weights)){
    weights <- rep(1, p)
  }

  if (is.null(beta.init)) {
    # temp <- EM_Init(X, Y, R, lambda = lambda, weights = weights,
    #                        beta.init = NULL, 
    #                        max.iter.EM = max.iter.EM, 
    #                        obj.tol = obj.tol, gamma0 = gamma0, 
    #                        epochs = 1) 
    # beta.list <- temp$beta.list
    beta.list <- list()
    for (j in 1:M) {
      beta.list[[j]] <- list()
      for (k in 1:R) {
        beta.list[[j]][[k]] <- matrix(runif(p*ncol(Y[[j]]), -0.1, 0.1), nrow=p, ncol=ncol(Y[[j]]))
      }
    }
  } else {
    beta.list <- beta.init
    beta.init <- NULL
  }
  
  
  # ------------------------------
  # Set initial values for delta
  # ------------------------------
  if(is.null(delta.init)){
    deltas <- seq(1, 1, length=R); 
    deltas <- deltas/sum(deltas)
  } else {
    deltas <- delta.init
  }
  
  # -------------------------------------
  # Could standardize X here, if needed
  # -------------------------------------
  Xin <- X
  
  # -------------------------
  # EM initialization
  # -------------------------
  k.iter <- 1
  iterating <- TRUE
  obj.prev <- obs.data.loglik(Y, Xin, beta.list, deltas, R, lambda, weights)
  S <- NULL
  
  while (iterating) {
    
    # -------------------------------
    # E-step 
    # -------------------------------
    if (k.iter == 1) {
      pis <- pis.update(Y, Xin, beta.list, deltas, R)
    } else {
      pis <- pis.update.S(Y, S, deltas, R, indices.set)
    }
    
    # --------------------------------
    # M-step
    # --------------------------------
    temp <- beta.update(Y, Xin, pis, beta.list, R, S, lambda, weights, gamma0, epochs=epochs)
    beta.list.new <- temp$beta.list
    indices.set <- temp$indices.set
    S <- temp$S # beta.list is initializer
    temp <- NULL
    deltas.new <- delta.update(pis)
    
    # ---------------------------------
    # Check convergence
    # ----------------------------------
    obj.current <- obs.data.loglik.S(Y, S, beta.list.new, deltas.new, R, lambda, weights, indices.set)
    if (k.iter %% 100 == 0) {
      cat(obj.current, "\n")
    }
    residual <- obj.current - obj.prev
    if (residual < obj.tol) {
      iterating <- FALSE
    } else {
      k.iter <- k.iter + 1
      if (k.iter > max.iter.EM) {
        iterating <- FALSE
        warning("EM algorithm did not converge!")
      }
    }
    
    # ----------------------------------
    # Update iterates
    # ----------------------------------
    deltas <- deltas.new
    beta.list <- beta.list.new
    obj.prev <- obj.current
    
  } 
  
  return(list(
    "pis" = pis,
    "deltas" = deltas,
    "beta.list" = beta.list))
  
}




pis.update <- function(Y, Xin, beta.list, deltas, R) {
  
  eval.logf.univariate <- function(Ymat, Xin, beta){
    lleval <- exp(crossprod(t(Xin), beta))
    Pmat <- lleval/(tcrossprod(colSums(t(lleval)), rep(1, dim(Ymat)[2])))
    inner <- log(colSums(t(Ymat*Pmat)))
    return(inner)
  }
  
  eval.f.overall <- function(Y, Xin, beta.list, R) {
    
    fr <- matrix(0, nrow = dim(Xin)[1], ncol = R)
    M <- length(Y)
    
    for (k in 1:R) {
      t1 <- rep(0, dim(Xin)[1])
      for (j in 1:M) {
        t1 <- t1 + eval.logf.univariate(Y[[j]], Xin, beta.list[[j]][[k]])
      }
      fr[,k] <- exp(t1)
    }
    
    return(fr)
    
  }
  
  temp <- eval.f.overall(Y, Xin, beta.list, R)
  t0 <- tcrossprod(rep(1,dim(Xin)[1]),deltas)*temp
  pis <- t0/colSums(t(t0))
  return(pis)
  
}




pis.update.S <- function(Y, S, deltas, R, indices.set) {
  
  eval.f.univariate <- function(Ymat, S) {
    lleval <- exp(S)
    Pmat <- lleval/(tcrossprod(rowSums(lleval), rep(1, dim(Ymat)[2])))
    return(rowSums(Ymat*Pmat))
  }
  
  eval.f.overall <- function(Y, S, R, indices.set) {
    
    fr <- matrix(1, nrow = dim(Y[[1]])[1], ncol = R)
    M <- length(Y)
    uu <- 1
    for (k in 1:M) {
      for (j in 1:R) {
        fr[,j] <- fr[,j]*eval.f.univariate(Y[[k]], S[,indices.set[[uu]]])
        uu <- uu + 1
      }
    }
    return(fr)
  }
  
  temp <- eval.f.overall(Y, S, R, indices.set)
  t0 <- tcrossprod(rep(1, dim(Y[[1]])[1]),deltas)*temp
  pis <- t0/colSums(t(t0))
  return(pis)
  
}




delta.update <- function(pis) {
  colSums(pis)/(dim(pis)[1])
}


beta.update <- function(Y, Xin, pis, beta.list,  R, S, lambda, weights, gamma0 = 1, epochs = epochs) {
  
  p <- dim(Xin)[2]
  M <- length(Y)
  cats <- unlist(lapply(Y, ncol))
  indices.set <- list()
  uu <- 1
  for (j in 1:M) {
    if (j == 1) {
      indices.set[[uu]] <- 1:cats[j]
    }
    for (k in 1:R) {
      if (k != R) {
        indices.set[[uu+1]] <- indices.set[[uu]][length(indices.set[[uu]])] + 1:cats[j]
      } else {
        if (j != M) {
          indices.set[[uu+1]] <- indices.set[[uu]][length(indices.set[[uu]])] + 1:cats[j + 1]
        }
      }
      uu <- uu + 1
    }
  }
  
  # -------------------------------------
  # Construct S_{rm} = XB_{rm} matrices
  # -------------------------------------
  beta.list.out <- list()
  if(is.null(S)){
    S <- matrix(0, nrow=dim(Xin)[1], ncol = R*sum(cats))
    uu <- 1
    for (j in 1:M) {
      for (k in 1:R) {
        S[, indices.set[[uu]]] <- crossprod(t(Xin), beta.list[[j]][[k]])
        uu <- uu + 1
      }
    }
  }
  
  eval.logf.S <- function(Ymat, S){
    lleval <- exp(S)
    Pmat <- lleval/(tcrossprod(rowSums(lleval), rep(1, dim(Ymat)[2])))
    return(log(rowSums(Ymat*Pmat)))
  }
  
  eval.logf.overall.S <- function(Y, S, R, pis, indices.set) {
    M <- length(Y)
    t1 <- 0
    uu <- 1
    for (k in 1:M) {
      for (j in 1:R) {
        t1 <- t1 - sum(pis[,j]*eval.logf.S(Y[[k]], S[,indices.set[[uu]]]))
        uu <- uu + 1
      }
    }
    return(t1)
  }
  
  grad.j <- function(j, Y, Xin, pis, R, M, S, indices.set){
    cats <- unlist(lapply(Y, ncol))
    n <- dim(Xin)[1]
    P.minusY <- matrix(0, nrow=dim(Xin)[1], ncol = R*sum(cats))
    uu <- 1
    for (kk in 1:M) {
      for (jj in 1:R) {
        t0 <- exp(S[,indices.set[[uu]]])
        P.minusY[,indices.set[[uu]]] <- pis[,jj]*((t0/tcrossprod(rowSums(t0), rep(1, cats[kk]))) - Y[[kk]])
        uu <- uu + 1
      }
    }
    return(crossprod(Xin[,j], P.minusY))
  }
  
  obj.prev <- eval.logf.overall.S(Y, S, R, pis, indices.set)
  obj.up <- obj.prev
  t0 <- sample(1:p)

  for (lll in 1:epochs) {
    for (jj in t0) {
      
      # --------------------------------------------------------
      # Obtain one proximal gradient descent update for jjth predictor
      # --------------------------------------------------------
      beta.prev <- rep(0, R*sum(cats))
      uu <- 1
      for(j in 1:M){
        for(k in 1:R){
          beta.prev[indices.set[[uu]]] <- beta.list[[j]][[k]][jj,]
          uu <- uu + 1
        }
      }
      
      # -------------------------------------
      # compute objective at beta.temp
      # -------------------------------------
      # if (jj == t0[1]) {
      #   obj.prev <- eval.logf.overall.S(Y, S, R, pis, indices.set)
      # } else {
      #   obj.prev <- obj.up # from last update
      # }
      obj.prev <- obj.up
      gamma <- gamma0
      linesearch <- TRUE
      temp <- c(grad.j(jj, Y, Xin, pis, R, M, S, indices.set))	
      
      while(linesearch){
        
        # -------------------------------------
        # get proximal gradient update
        # -------------------------------------
        u <- beta.prev - gamma*temp
        beta.up <- max(1 - (weights[jj]*lambda*gamma)/sqrt(sum(u^2)), 0)*u
        
        # --------------------------------------
        # check line search condition
        # --------------------------------------
        S.up <- S
        uu <- 1
        for (j in 1:M) {
          for (k in 1:R) {
            S.up[,indices.set[[uu]]] <- S.up[,indices.set[[uu]]] - tcrossprod(Xin[,jj], beta.list[[j]][[k]][jj,] - beta.up[indices.set[[uu]]]) 
            uu <- uu + 1
          }
        }
        
        obj.up <- eval.logf.overall.S(Y, S.up, R, pis, indices.set)
        if(is.na(obj.up)){
          obj.up <- Inf
        }
        if (obj.up <= obj.prev + sum(temp*(beta.up - beta.prev)) + (1/(2*gamma))*sum((beta.up - beta.prev)^2)) {
          linesearch <- FALSE
        } else {
          gamma <- gamma/2
        }
      }
      
      # ------------------------------------
      # Update S and beta with new iterates
      # ------------------------------------
      S <- S.up
      uu <- 1
      for (j in 1:M) {
        for (k in 1:R) {
          beta.list[[j]][[k]][jj,] <-  beta.up[indices.set[[uu]]]
          uu <- uu + 1
        }
      }
    }
  }
  
  return(list("beta.list" = beta.list, "S" = S, "indices.set" = indices.set))
  
}


obs.data.loglik <- function(Y, Xin, beta.list, deltas, R, lambda, weights){
  
  eval.f.univariate <- function(Ymat, X, beta) {
    lleval <- exp(crossprod(t(X), beta))
    Pmat <- lleval/(tcrossprod(rowSums(lleval), rep(1, dim(Ymat)[2])))
    return(rowSums(Ymat*Pmat))
  }
  
  eval.f.overall <- function(Y, X, beta.list, R) {
    
    fr <- matrix(0, nrow = dim(X)[1], ncol = R)
    M <- length(Y)
    
    for (k in 1:R) {
      t1 <- rep(1, dim(X)[1])
      for (j in 1:M) {
        t1 <- t1*eval.f.univariate(Y[[j]], X, beta.list[[j]][[k]])
      }
      fr[,k] <- t1
    }
    return(fr)
  }
  
  temp <- eval.f.overall(Y, Xin, beta.list, R)
  t0 <- log(rowSums((rep(1,dim(Xin)[1])%*%t(deltas))*temp))
  pen <- 0
  M <- length(Y)
  pen <- lambda*sum(weights*apply(Reduce(cbind, (lapply(beta.list, function(x){Reduce(cbind, x)}))), 1, function(x){sqrt(sum(x^2))}))
  return(sum(t0) - pen)
}

obs.data.loglik.S <- function(Y, S, beta.list, deltas, R, lambda, weights, indices.set){
  
  eval.f.univariate <- function(Ymat, S) {
    lleval <- exp(S)
    Pmat <- lleval/(tcrossprod(rowSums(lleval), rep(1, dim(Ymat)[2])))
    return(rowSums(Ymat*Pmat))
  }
  
  eval.f.overall <- function(Y, S, R, indices.set) {
    
    fr <- matrix(1, nrow = dim(Y[[1]])[1], ncol = R)
    M <- length(Y)
    uu <- 1
    for (k in 1:M) {
      for (j in 1:R) {
        fr[,j] <- fr[,j]*eval.f.univariate(Y[[k]], S[,indices.set[[uu]]])
        uu <- uu + 1
      }
    }
    return(fr)
  }
  
  temp <- eval.f.overall(Y, S, R, indices.set)
  t0 <- log(rowSums((rep(1,dim(Y[[1]])[1])%*%t(deltas))*temp))
  pen <- 0
  M <- length(Y)
  pen <- lambda*sum(weights*apply(Reduce(cbind, (lapply(beta.list, function(x){Reduce(cbind, x)}))), 1, function(x){sqrt(sum(x^2))}))
  return(sum(t0) - pen)
}







predict.EM <- function(Xnew, EMfit) {
  
  prob.univariate <- function(X, beta){
    lleval <- exp(crossprod(t(X), beta))
    Pmat <- lleval/(tcrossprod(rowSums(lleval), rep(1, dim(beta)[2])))
    return(Pmat)
  }
  
  prob.overall <- function(X, EMfit, R) {
    
    M <- length(EMfit$beta)
    R <- length(EMfit$beta[[1]])
    Ydims <- unlist(lapply(lapply(EMfit$beta, `[[`, 1), 
                           function(x){dim(x)[2]}))
    prob.tensor <- array(0, dim=c(dim(X)[1],Ydims))
    
    for (k in 1:R) {
      
      probs <- list()
      for (j in 1:M) {
        probs[[j]] <- prob.univariate(X, EMfit$beta[[j]][[k]])
      }
      t0.array <- NULL
      
      for(i in 1:dim(X)[1]){
        t0 <- NULL
        for(j in 1:M){
          if(j == 1){
            t0 <- probs[[j]][i,]
          } else {
            t0 <- t0 %o% probs[[j]][i,]
          }
        }
        if(i > 1){
          t0.array <- abind(t0.array, t0, along=1)
        } else {
          t0.array <- abind(t0.array, t0, along=0)
        }
        
      }
      
      prob.tensor <- prob.tensor + EMfit$deltas[k]*t0.array
      
    }
    
    return(prob.tensor)
  }
  
  temp <- prob.overall(Xnew, EMfit)
  
  return(list(
    "probs" = temp,
    "preds" = t(apply(temp, 1, function(x){which(x == max(x), arr.ind=TRUE)}))))	
}





















# EM_Init <- function(X, Y, R, lambda, weights = NULL,
#                            beta.init = NULL, 
#                            max.iter.EM = 1000, 
#                            obj.tol = 1e-8, gamma0 = 1, 
#                            epochs = 1) {
  
#   # -----------------------
#   # Preliminaries
#   # -----------------------
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   M <- length(Y)
  
#   # ------------------------------
#   # Set initial values for beta
#   # ------------------------------
#   if (is.null(beta.init)) {
#     beta.list <- list()
#     for (j in 1:M) {
#       beta.list[[j]] <- list()
#       for (k in 1:R) {
#         beta.list[[j]][[k]] <- matrix(runif(p*ncol(Y[[j]]), -0.1, 0.1), nrow=p, ncol=ncol(Y[[j]]))
#       }
#     }
#   } else {
#     beta.list <- beta.init
#     beta.init <- NULL
#   }
  
#   if (is.null(weights)){
#     weights <- rep(1, p)
#   }
  
#   # ------------------------------
#   # Set initial values for delta
#   # ------------------------------
#   deltas <- seq(1, 1, length=R); 
#   deltas <- deltas/sum(deltas)
  
#   # -------------------------------------
#   # Could standardize X here, if needed
#   # -------------------------------------
#   Xin <- X
  
#   # -------------------------
#   # EM initialization
#   # -------------------------
#   k.iter <- 1
#   iterating <- TRUE
#   obj.prev <- obs.data.loglik(Y, Xin, beta.list, deltas, R, lambda, weights)
#   S <- NULL
  
#   while (iterating) {
    
#     # -------------------------------
#     # E-step 
#     # -------------------------------
#     if (k.iter == 1) {
#       pis <- pis.update(Y, Xin, beta.list, deltas, R)
#     } else {
#       pis <- pis.update.S(Y, S, deltas, R, indices.set)
#     }
    
#     # --------------------------------
#     # M-step
#     # --------------------------------
#     temp <- beta.update(Y, Xin, pis, beta.list, R, S, lambda, weights, gamma0, epochs=epochs)
#     beta.list.new <- temp$beta.list
#     indices.set <- temp$indices.set
#     S <- temp$S # beta.list is initializer
#     temp <- NULL
#     #deltas.new <- delta.update(pis)
#     deltas.new <- deltas
    
#     # ---------------------------------
#     # Check convergence
#     # ----------------------------------
#     obj.current <- obs.data.loglik.S(Y, S, beta.list.new, deltas.new, R, lambda, weights, indices.set)
#     if (k.iter %% 100 == 0) {
#       cat(obj.current, "\n")
#     }
#     residual <- obj.current - obj.prev
#     if (residual < obj.tol) {
#       iterating <- FALSE
#     } else {
#       k.iter <- k.iter + 1
#       if (k.iter > max.iter.EM) {
#         iterating <- FALSE
#         warning("EM algorithm did not converge!")
#       }
#     }
    
#     # ----------------------------------
#     # Update iterates
#     # ----------------------------------
#     deltas <- deltas.new
#     beta.list <- beta.list.new
#     obj.prev <- obj.current
    
#   } 
  
#   return(list(
#     "pis" = pis,
#     "deltas" = deltas,
#     "beta.list" = beta.list))
  
# }




