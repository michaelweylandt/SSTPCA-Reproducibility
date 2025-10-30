Rcpp::sourceCpp("../../../tensor_factorizations.cpp")
library(dplyr)
library(reshape2)
library(ggplot2)


## THINGS TO SHOW IN SIMULATION
# - Consistency of MF-SST-PCA
# - Rank Selection Effectiveness for (r, K) -- this simulation
# - Misspecification

# For this plot, since we're messing with the rank, our error target will
# be |X_star - X_hat|_F^2
#
# We consider three approaches
# 1) Oracle K, r_K
# 2) Oracle K, BIC - r_k
# 3) MSE K, BIC - r_k

TGRID <- c(50, 100, 150, 200, 250, 300)#, 400, 500)
p <- 50
SDGRID <- c(0.1, 0.2, 0.5, 1, 2)
NREPS <- 30

ALL_ERROR <- data.frame()

DMIN_GRID <- seq(1, 12, by=1)
DFACTORS <- c(3, 2, 1)
K <- length(DFACTORS)
r <- c(3, 2, 5)
rr <- cumsum(r)
R <- sum(r)

TSCALE <- TRUE

# Oracle initialization for the oracle case - best possible here
for(dmin in DMIN_GRID){
    for(sd in SDGRID){
        for(T in TGRID){
            cat("T=", T, "scenario = DOUBLE ORACLE sd=", sd, "dmin=", dmin, "\n")
            ## Generate X from a rank-3 SST-PCA model

            d <- dmin * DFACTORS

            if(TSCALE)  d <- d * sqrt(T) / sqrt(min(TGRID))

            ERROR_U <- matrix(NA, ncol=K, nrow=NREPS)
            ERROR_V <- matrix(NA, ncol=K, nrow=NREPS)

            vnorm <- function(x) sqrt(sum(x^2))
            vec_angle <- function(u, v) acos(abs(sum(u * v)) / (vnorm(u) * vnorm(v))) / pi * 180

            mat_angle <- function(U, V) acos(min(svd(crossprod(U, V), nu=0, nv=0)$d)) / pi * 180
            rQ <- function(p, k) qr.Q(qr(matrix(rnorm(p * k), nrow=p)))

            ss_noise <- function(T, p){
                replicate(T, simplify = "array", {
                    E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                    E <- (E + t(E)) / 2
                })
            }

            for(rep in seq(NREPS)){
                U <- rQ(T, K)
                V <- rQ(p, R)

                FACTORS <- lapply(1:K, function(i)
                    list(u=U[,i,drop=FALSE],
                         V = if(i == 1) V[,1:rr[1],drop=FALSE] else V[,(rr[i-1]+1):rr[i],drop=FALSE],
                         d = d[i]))
                X_TERMS <- lapply(FACTORS, function(f)f$d * (tcrossprod(f$V) %o% as.vector(f$u)))

                X_STAR <- Reduce(`+`, X_TERMS)
                E <- sd * ss_noise(T, p)
                X_OBS  <- X_STAR + E

                X_WORKING <- X_OBS
                FITS <- vector("list", K)
                for(k in 1:K){
                    FITS[[k]] <- ss_tpm(X=X_WORKING, rank=r[k], u_init=U[,k], u_init_strategy=2)
                    X_WORKING <- X_WORKING - FITS[[k]]$X_hat
                }

                r_hats <- r

                X_HAT <- Reduce(`+`, lapply(FITS, function(f) f$X_hat))
                X_ERR <- mean((X_HAT - X_STAR)^2)
                X_RERR <- X_ERR / mean(X_STAR^2)

                ALL_ERROR <- rbind(ALL_ERROR,
                                   data.frame(X_ERR = X_ERR,
                                              X_RERR=X_RERR,
                                              T = T,
                                              p = p,
                                              sd=sd,
                                              rep=rep,
                                              dmin=dmin,
                                              strategy="DOUBLE_ORACLE",
                                              r_hats=r_hats,
                                              r=r,
                                              r_ix=1:3,
                                              K=3,
                                              K_hat=3))
            }
        }

        saveRDS(ALL_ERROR, "select_rank_BIC_error.rds")
    }
}

R_MAX <- 10


for(dmin in DMIN_GRID){
    for(sd in SDGRID){
        for(T in TGRID){
            cat("T=", T, "scenario = ORACLE INIT sd=", sd, "dmin=", dmin, "\n")
            ## Generate X from a rank-3 SST-PCA model

            d <- dmin * DFACTORS
            if(TSCALE)  d <- d * sqrt(T) / sqrt(min(TGRID))

            ERROR_U <- matrix(NA, ncol=K, nrow=NREPS)
            ERROR_V <- matrix(NA, ncol=K, nrow=NREPS)

            vnorm <- function(x) sqrt(sum(x^2))
            vec_angle <- function(u, v) acos(abs(sum(u * v)) / (vnorm(u) * vnorm(v))) / pi * 180

            mat_angle <- function(U, V) acos(min(svd(crossprod(U, V), nu=0, nv=0)$d)) / pi * 180
            rQ <- function(p, k) qr.Q(qr(matrix(rnorm(p * k), nrow=p)))

            ss_noise <- function(T, p){
                replicate(T, simplify = "array", {
                    E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                    E <- (E + t(E)) / 2
                })
            }

            for(rep in seq(NREPS)){
                U <- rQ(T, K)
                V <- rQ(p, R)

                FACTORS <- lapply(1:K, function(i)
                    list(u=U[,i,drop=FALSE],
                         V = if(i == 1) V[,1:rr[1],drop=FALSE] else V[,(rr[i-1]+1):rr[i],drop=FALSE],
                         d = d[i]))
                X_TERMS <- lapply(FACTORS, function(f)f$d * (tcrossprod(f$V) %o% as.vector(f$u)))

                X_STAR <- Reduce(`+`, X_TERMS)
                E <- sd * ss_noise(T, p)
                X_OBS  <- X_STAR + E

                X_WORKING <- X_OBS
                FITS <- vector("list", K)
                r_hats <- numeric(K)
                for(k in 1:K){
                    BIC <- numeric(R_MAX)
                    for(r_hat in seq(R_MAX)){
                        fit_hat <-  ss_tpm(X=X_WORKING, rank=r_hat, u_init_strategy=2, u_init=U[,k])
                        sse <- sum((fit_hat$X_hat - X_WORKING)^2)
                        BIC[r_hat] <- p^2 * T * log(sse) + p * r_hat * log(p^2 * T)
                    }
                    r_hat_BIC <- which.min(BIC)
                    r_hats[k] <- r_hat_BIC
                    FITS[[k]] <- ss_tpm(X=X_WORKING, rank=r_hat_BIC, u_init=U[,k], u_init_strategy=2)
                    X_WORKING <- X_WORKING - FITS[[k]]$X_hat
                }

                X_HAT <- Reduce(`+`, lapply(FITS, function(f) f$X_hat))
                X_ERR <- mean((X_HAT - X_STAR)^2)
                X_RERR <- X_ERR / mean(X_STAR^2)

                ALL_ERROR <- rbind(ALL_ERROR,
                                   data.frame(X_ERR = X_ERR,
                                              X_RERR=X_RERR,
                                              T = T,
                                              p = p,
                                              rep=rep,
                                              sd=sd,
                                              dmin=dmin,
                                              strategy="INIT_ORACLE",
                                              r_hats=r_hats,
                                              r=r,
                                              r_ix=1:3,
                                              K=3,
                                              K_hat=3))
            }
        }

        saveRDS(ALL_ERROR, "select_rank_BIC_error.rds")
    }
}

# No oracle usage here, but we know K
for(dmin in DMIN_GRID){
    for(sd in SDGRID){
        for(T in TGRID){
            cat("T=", T, "STRATEGY=NO_ORACLE_KNOWN_K sd=", sd, "dmin=", dmin, "\n")
            ## Generate X from a rank-3 SST-PCA model

            d <- dmin * DFACTORS
            if(TSCALE)  d <- d * sqrt(T) / sqrt(min(TGRID))

            ERROR_U <- matrix(NA, ncol=K, nrow=NREPS)
            ERROR_V <- matrix(NA, ncol=K, nrow=NREPS)

            vnorm <- function(x) sqrt(sum(x^2))
            vec_angle <- function(u, v) acos(abs(sum(u * v)) / (vnorm(u) * vnorm(v))) / pi * 180

            mat_angle <- function(U, V) acos(min(svd(crossprod(U, V), nu=0, nv=0)$d)) / pi * 180
            rQ <- function(p, k) qr.Q(qr(matrix(rnorm(p * k), nrow=p)))

            ss_noise <- function(T, p){
                replicate(T, simplify = "array", {
                    E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                    E <- (E + t(E)) / 2
                })
            }

            for(rep in seq(NREPS)){
                U <- rQ(T, K)
                V <- rQ(p, R)

                FACTORS <- lapply(1:K, function(i)
                    list(u=U[,i,drop=FALSE],
                         V = if(i == 1) V[,1:rr[1],drop=FALSE] else V[,(rr[i-1]+1):rr[i],drop=FALSE],
                         d = d[i]))
                X_TERMS <- lapply(FACTORS, function(f)f$d * (tcrossprod(f$V) %o% as.vector(f$u)))

                X_STAR <- Reduce(`+`, X_TERMS)
                E <- sd * ss_noise(T, p)
                X_OBS  <- X_STAR + E

                X_WORKING <- X_OBS
                FITS <- vector("list", K)
                r_hats <- numeric(K)
                for(k in 1:K){
                    BIC <- numeric(R_MAX)
                    for(r_hat in seq(R_MAX)){
                        fit_hat <-  ss_tpm(X=X_WORKING, rank=r_hat, u_init_strategy=0)
                        sse <- sum((fit_hat$X_hat - X_WORKING)^2)
                        BIC[r_hat] <- p^2 * T * log(sse) + p * r_hat * log(p^2 * T)
                    }
                    r_hat_BIC <- which.min(BIC)
                    r_hats[k] <- r_hat_BIC
                    FITS[[k]] <- ss_tpm(X=X_WORKING, rank=r_hat_BIC)
                    X_WORKING <- X_WORKING - FITS[[k]]$X_hat
                }

                X_HAT <- Reduce(`+`, lapply(FITS, function(f) f$X_hat))
                X_ERR <- mean((X_HAT - X_STAR)^2)
                X_RERR <- X_ERR / mean(X_STAR^2)

                ALL_ERROR <- rbind(ALL_ERROR,
                                   data.frame(X_ERR = X_ERR,
                                              X_RERR=X_RERR,
                                              T = T,
                                              p = p,
                                              rep=rep,
                                              sd=sd,
                                              dmin=dmin,
                                              strategy="R_BIC",
                                              r_hats=r_hats,
                                              r=r,
                                              r_ix=1:3,
                                              K=3,
                                              K_hat=3))
            }
        }

        saveRDS(ALL_ERROR, "select_rank_BIC_error.rds")
    }
}




# Hardest Case - don't even know K here
for(dmin in DMIN_GRID){
    for(sd in SDGRID){
        for(T in TGRID){
            cat("T=", T, "STRATEGY=NO_ORACLE_UNKNOWN_K sd=", sd, "dmin=", dmin, "\n")
            ## Generate X from a rank-3 SST-PCA model

            d <- dmin * DFACTORS
            if(TSCALE)  d <- d * sqrt(T) / sqrt(min(TGRID))

            ERROR_U <- matrix(NA, ncol=K, nrow=NREPS)
            ERROR_V <- matrix(NA, ncol=K, nrow=NREPS)

            vnorm <- function(x) sqrt(sum(x^2))
            vec_angle <- function(u, v) acos(abs(sum(u * v)) / (vnorm(u) * vnorm(v))) / pi * 180

            mat_angle <- function(U, V) acos(min(svd(crossprod(U, V), nu=0, nv=0)$d)) / pi * 180
            rQ <- function(p, k) qr.Q(qr(matrix(rnorm(p * k), nrow=p)))

            ss_noise <- function(T, p){
                replicate(T, simplify = "array", {
                    E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                    E <- (E + t(E)) / 2
                })
            }

            for(rep in seq(NREPS)){
                U <- rQ(T, K)
                V <- rQ(p, R)

                FACTORS <- lapply(1:K, function(i)
                    list(u=U[,i,drop=FALSE],
                         V = if(i == 1) V[,1:rr[1],drop=FALSE] else V[,(rr[i-1]+1):rr[i],drop=FALSE],
                         d = d[i]))
                X_TERMS <- lapply(FACTORS, function(f)f$d * (tcrossprod(f$V) %o% as.vector(f$u)))

                X_STAR <- Reduce(`+`, X_TERMS)
                E <- sd * ss_noise(T, p)
                X_OBS  <- X_STAR + E

                X_MAT <- array(X_OBS, dim=c(p^2, T));
                X_MAT_d <- svd(X_MAT, nu=0, nv=0)$d;
                d_thresh <- 2.858 * median(X_MAT_d)
                K_hat <- sum(X_MAT_d >= d_thresh)

                if(K_hat >= 1){

                    X_WORKING <- X_OBS
                    FITS <- vector("list", K_hat)
                    r_hats <- numeric(K_hat)
                    for(k in 1:K_hat){
                        BIC <- numeric(R_MAX)
                        for(r_hat in seq(R_MAX)){
                            fit_hat <-  ss_tpm(X=X_WORKING, rank=r_hat, u_init_strategy=0)
                            sse <- sum((fit_hat$X_hat - X_WORKING)^2)
                            BIC[r_hat] <- p^2 * T * log(sse) + p * r_hat * log(p^2 * T)
                        }
                        r_hat_BIC <- which.min(BIC)
                        r_hats[k] <- r_hat_BIC
                        FITS[[k]] <- ss_tpm(X=X_WORKING, rank=r_hat_BIC)
                        X_WORKING <- X_WORKING - FITS[[k]]$X_hat
                    }

                    X_HAT <- Reduce(`+`, lapply(FITS, function(f) f$X_hat))
                    X_ERR <- mean((X_HAT - X_STAR)^2)
                    X_RERR <- X_ERR / mean(X_STAR^2)

                    ALL_ERROR <- rbind(ALL_ERROR,
                                       data.frame(X_ERR = X_ERR,
                                                  X_RERR=X_RERR,
                                                  T = T,
                                                  p = p,
                                                  rep=rep,
                                                  sd=sd,
                                                  dmin=dmin,
                                                  strategy="K_MSE_R_BIC",
                                                  r_hats=r_hats,
                                                  r=r[seq(K_hat)],
                                                  r_ix=seq(K_hat),
                                                  K=3,
                                                  K_hat=K_hat))
                } else {
                    ALL_ERROR <- rbind(ALL_ERROR,
                                       data.frame(X_ERR=mean(X_STAR^2),
                                                  X_RERR=1,
                                                  T = T,
                                                  p = p,
                                                  rep=rep,
                                                  sd=sd,
                                                  dmin=dmin,
                                                  strategy="K_MSE_R_BIC",
                                                  r_hats=0,
                                                  r=0,
                                                  r_ix=0,
                                                  K=3,
                                                  K_hat=0))
                }
            }
        }

        saveRDS(ALL_ERROR, "select_rank_BIC_error.rds")
    }
}



saveRDS(ALL_ERROR, "select_rank_BIC_error.rds")
