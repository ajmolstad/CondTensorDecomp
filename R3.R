# -----------------------------------------------------
# Set slurm array ID
# -----------------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


# -----------------------------------------------------
# queue results list
# -----------------------------------------------------
namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}


library(abind)
library(glmnet)
source("~/blue/LatentIndep/Functions/Group_EM.R")

R <- 3
p <- 100
nreps <- 100
t1 <- expand.grid("n" = rep(c(150, 300, 450), each=nreps),
	"delta1" = seq(1/3, 2/3, length=5)[1:4], 
	"sigma" = c(1,2))
n <- t1[uu,1]
delta <- c(t1[uu,2], 1 - (1/3 + t1[uu,2]), 1/3)
nval <- 200
savename <- paste("~/blue/LatentIndep/Simulations/R3/Results/Results_", uu, ".RDS", sep="")


set.seed(uu)
SigmaX <- matrix(0, p, p)
for(j in 1:p){
  for(k in 1:p){
    SigmaX[j,k] <- .5^abs(j-k)
  }
}
eo <- eigen(SigmaX)
SigmaXsqrt <- eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%SigmaXsqrt
Xval <- matrix(rnorm(nval*p), nrow=nval, ncol=p)%*%SigmaXsqrt
Ydims <- c(4, 4, 4, 4)
M <- length(Ydims)
preds <- sample(1:p, 5)
beta <- list()
for(j in 1:M){
	beta[[j]] <- list()
	for(k in 1:R){
		beta[[j]][[k]] <- matrix(0, p, Ydims[j])
		beta[[j]][[k]][preds,] <- rnorm(Ydims[j], sd = t1[uu,3])
	}
}

ntest <- 1000
Xnew <- matrix(rnorm(ntest*p), nrow=ntest, ncol=p)%*%SigmaXsqrt

makeDat <- list()
makeDat$beta.list <- beta
makeDat$deltas <- delta

# -------------------------------
# Generate training data
# -------------------------------
t0 <- predict.EM(X, makeDat)
Y <- list()
for(k in 1:length(Ydims)){
	Y[[k]] <- matrix(0, nrow=n, ncol=Ydims[k])
}

for(j in 1:n){
	pr <- t0$probs[j,,,,]
	temp <- sample(1:prod(Ydims), 1, prob =c(pr))
	observed <- which(t0$probs[j,,,,] == pr[temp], arr.ind=TRUE)
	for(k in 1:length(Ydims)){
		Y[[k]][j,observed[k]] <- 1
	}
}

# -------------------------------
# Generate validation data
# -------------------------------
t0 <- predict.EM(Xval, makeDat)
Yval <- list()
for(k in 1:length(Ydims)){
	Yval[[k]] <- matrix(0, nrow=nval, ncol=Ydims[k])
}

for(j in 1:nval){
	pr <- t0$probs[j,,,,]
	temp <- sample(1:prod(Ydims), 1, prob =c(pr))
	observed <- which(t0$probs[j,,,,] == pr[temp], arr.ind=TRUE)
	for(k in 1:length(Ydims)){
		Yval[[k]][j,observed[k]] <- 1
	}
}

# -----------------------------------------------
# Get probabilities from testing data for KL 
# -----------------------------------------------
pop <- predict.EM(Xnew, makeDat)
Results <- namedList(beta, delta, n)


source("~/blue/LatentIndep/Simulations/R3/Main_Random.R")
