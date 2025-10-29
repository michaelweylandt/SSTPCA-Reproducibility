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
# 3) Deflation K, BIC - r_k

TGRID <- c(50, 100, 150, 200, 250, 300, 400, 500)
p <- 50
SDGRID <- c(0.1, 0.2, 0.5, 1)
NREPS <- 25

ALL_ERROR <- data.frame()

DMIN_GRID <- seq(1, 12, by=1)

PROBLEMS <- list(
    P1 = list(r = c(1, 2, 3),
              dfactors = c(1, 2, 3)),
    P2 = list(r = c(5),
              dfactors = c(1)),
    P3 = list(r = c(3, 2, 1),
              dfactors = c(1, 2, 3)),
    P4 = list(r = c(5, 5, 5),
              dfactors = c(1, 2, 3)),
    P5 = list(r = c(10, 1, 1),
              dfactors = c(1, 2, 3)),
    P6 = list(r = c(10, 1, 1, 10),
              dfactors = c(1, 2, 3, 4)),
    P7 = list(r = c(10, 1, 1, 10),
              dfactors = c(1, 2, 2, 1)),
    P8 = list(r = c(10, 1, 5, 1, 10),
              dfactors = c(1, 2, 3, 2, 1)),
    P9 = list(r = c(5, 10),
              dfactors = c(1, 3)),
    P10 = list(r = c(5, 10),
               dfactors = c(2, 2))
)

for(sd in SDGRID){
    for(T in TGRID){
        for(dmin in DMIN_GRID){
            for(p_ix in seq_along(PROBLEMS)){
                r <- PROBLEMS[[p_ix]]$r
                dfactors <- PROBLEMS[[p_ix]]$dfactors

                K <- length(r)
                rr <- cumsum(r)
                R <- sum(r)

                cat("T=", T, "sd=", sd, "dmin=", dmin, "problem=", p_ix, "\n")
                d <- dmin * dfactors

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

                    normE <- sqrt(mean(E^2))
                    normX <- sapply(X_TERMS, function(x) sqrt(mean(x^2)))
                    min_normX <- min(normX)
                    normX_gt_normE <- sum(normX > normE)

                    X_MAT <- array(X_OBS, dim=c(p^2, T));
                    X_MAT_d <- svd(X_MAT, nu=0, nv=0)$d;
                    d_thresh <- 2.858 * median(X_MAT_d)
                    K_hat <- sum(X_MAT_d >= d_thresh)

                    ALL_ERROR <- rbind(ALL_ERROR,
                                       data.frame(K = K,
                                                  K_hat = K_hat,
                                                  K_err = abs(K - K_hat),
                                                  normE=normE,
                                                  min_normX=min_normX,
                                                  normX_gt_normE=normX_gt_normE,
                                                  T = T,
                                                  p = p,
                                                  problem=p_ix,
                                                  sd=sd,
                                                  dmin=dmin,
                                                  strategy="DonohoSVT_M3"))

                    saveRDS(ALL_ERROR, "select_K_error.rds")
                }
            }
        }
    }
}
