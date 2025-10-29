Rcpp::sourceCpp("../../../tensor_factorizations.cpp")
Rcpp::sourceCpp("../../../tensor_factorizations_simul.cpp")
library(dplyr)

TGRID <- c(50, 100, 150, 200, 250, 300)#, 400, 500)

ALL_ERROR <- data.frame()

NREPS <- 2

min_svd <- function(X) min(svd(X, nu=0, nv=0)$d)
mat_angle <- function(X, Y) min_svd(crossprod(X, Y)) / (min_svd(X) * min_svd(Y))

DMIN_GRID <- seq(1, 12, by=1)
DFACTORS <- c(3, 2, 1)
K <- length(DFACTORS)
r <- c(3, 2, 1)
rr <- cumsum(r)
R <- sum(r)
p <- 50
SDGRID <- c(0.1, 0.2, 0.5, 1)

INITIALIZATION <- c(0, 1, 2, -1, -2, -3)
SCENARIOS <- c(1, 2, 3)
ORTHO_Q <- c(TRUE, FALSE)

for(dmin in DMIN_GRID){
    for(sd in SDGRID){
        for(init in INITIALIZATION){
            for(T in TGRID){
                for(oq in ORTHO_Q){
                    for(scenario in SCENARIOS){
                        cat("SCENARIO=", scenario, "ORTHOQ=", oq, " T=", T, "init=", init, "sd=", sd, "dmin=", dmin, "\n")
                        ## Generate X from a rank-3 SST-PCA model

                        d <- dmin * DFACTORS

                        ERROR_U <- numeric(NREPS)
                        ERROR_V <- numeric(NREPS)
                        X_ERR <- numeric(NREPS)
                        X_RERR <- numeric(NREPS)

                        rQ <- function(p, k) qr.Q(qr(matrix(rnorm(p * k), nrow=p)))

                        ss_noise <- function(T, p){
                            replicate(T, simplify = "array", {
                                E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                                E <- (E + t(E)) / 2
                            })
                        }

                        for(rep in seq(NREPS)){
                            if(scenario == 1){
                                U <- rQ(T, K)
                                V <- rQ(p, R)
                            } else if (scenario == 2){
                                U <- replicate(K, rQ(T, 1)) %>% drop
                                V <- rQ(p, R)
                            } else if (scenario == 3){
                                U <- replicate(K, rQ(T, 1)) %>% drop
                                V <- do.call(cbind, lapply(r, function(.) rQ(p, .)))
                            }


                            FACTORS <- lapply(1:K, function(i)
                                list(u=U[,i,drop=FALSE],
                                     V = if(i == 1) V[,1:rr[1],drop=FALSE] else V[,(rr[i-1]+1):rr[i],drop=FALSE],
                                     d = d[i]))
                            X_TERMS <- lapply(FACTORS, function(f)f$d * (tcrossprod(f$V) %o% as.vector(f$u)))

                            X_STAR <- Reduce(`+`, X_TERMS)
                            E <- sd * ss_noise(T, p)
                            X_OBS  <- X_STAR + E

                            #orthog <- function(X) qr.Q(qr(X))
                            if(init >= 0){
                                FIT <- ss_tpm_simul(X=X_OBS, n_factors=3, ranks=r, u_init_strategy=init, u_init=U, ortho_q = oq)

                                ERROR_U[rep] <- mat_angle(FIT$u_hat, U)
                                ERROR_V[rep] <- mat_angle(do.call(cbind, FIT$v_hat), V)
                                X_ERR[rep] <- mean((Reduce(`+`, FIT$X_hat) - X_STAR)^2)
                                X_RERR[rep] <- mean(X_STAR^2)
                            } else {
                                # Need to handle the deflation case here...
                                X_WORKING <- X_OBS
                                FITS <- vector("list", K)
                                for(k in 1:K){
                                    FITS[[k]] <- ss_tpm(X=X_WORKING, rank=r[k], u_init=U[,k], u_init_strategy=abs(init)-1)
                                    X_WORKING <- X_WORKING - FITS[[k]]$X_hat
                                }

                                X_HAT_ACCUM <- X_OBS - X_WORKING

                                U_hat <- do.call(cbind, lapply(FITS, function(f) f$u_hat))
                                V_hat <- do.call(cbind, lapply(FITS, function(f) f$v_hat))

                                ERROR_U[rep] <- mat_angle(U_hat, U)
                                ERROR_V[rep] <- mat_angle(V_hat, V)

                                X_ERR[rep] <- mean((X_HAT_ACCUM - X_STAR)^2)
                                X_RERR[rep] <- mean((X_HAT_ACCUM - X_STAR)^2) / mean(X_STAR^2)
                            }
                        }

                        ERROR_DF <- tibble(cosU = ERROR_U,
                                           cosV = ERROR_V,
                                           X_ERR = X_ERR,
                                           X_RERR= X_RERR,
                                           OQ = oq,
                                           T=T,
                                           p = p,
                                           init=init,
                                           sd=sd,
                                           dmin=dmin,
                                           K=K,
                                           d=list(d),
                                           rep=1:NREPS,
                                           ranks=list(r),
                                           scenario=scenario)

                        ALL_ERROR <- rbind(ALL_ERROR, ERROR_DF)

                        saveRDS(ALL_ERROR, "multi_rank_error_simul.rds")
                    }
                }
            }
        }
    }
}
