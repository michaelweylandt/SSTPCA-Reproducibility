library(Rcpp)
library(tidyverse)

Rcpp::sourceCpp("sstpca_tracking.cpp")

vnorm <- function(x) sqrt(sum(x^2))
rmse <- function(x) vnorm(x) / length(x)
SVD <- function(x, ...){if(anyNA(x)) list(d = NA) else svd(x, ...)}

T <- 20;
p <- 200;

max_iter <- 60

u_star <- rep(1, times = T) / sqrt(T)

RESULTS <- tibble(T = integer(),
                  p = integer(),
                  rank = integer(),
                  d = numeric(),
                  iterations = integer(),
                  max_iter = integer(),
                  U_U_STAR = list(),
                  U_U_HAT  = list(),
                  V_V_STAR = list(),
                  V_V_HAT  = list(),
                  n = integer())

for(r in 1:5){
    d <- 15 * r^(-1/4)
    for(n in 1:100){
        generate_V_star <- function(p, rank){
            A_svd <- svd(matrix(rnorm(p^2), p, p), nu = rank, nv = 0)$u
        }

        ss_noise <- function(T, p){
            replicate(T, simplify = "array", {
                E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                E <- (E + t(E)) / 2
            })
        }

        V_star <- generate_V_star(p, r)

        X_star <- d * (tcrossprod(V_star) %o% u_star)
        X      <- X_star + ss_noise(T, p)

        generate_u_init <- function(T){
            u <- abs(rnorm(T))
            u / vnorm(u)
        }

        FIT <- ss_tpm(X, u_init = generate_u_init(T), rank = r, max_iter = max_iter)

        ## Compute accuracy measures
        u_star_conv <- zoo::na.locf(sapply(seq(max_iter), function(i) sum(FIT$U_history[i,] * u_star)))
        u_hat_conv  <- zoo::na.locf(sapply(seq(max_iter), function(i) sum(FIT$U_history[i,] * FIT$u_hat)))

        # Two robustness things here:
        # -- SVD = base::svd but aborting on NAs rather than erroring
        # -- min(d, 1) to prevent numerically unstable weirdness
        V_star_conv <- zoo::na.locf(sapply(seq(max_iter), function(i) min(c(SVD(crossprod(FIT$V_history[,,i], V_star), nu = 0, nv = 0)$d, 1))))
        V_hat_conv <- zoo::na.locf(sapply(seq(max_iter), function(i) min(c(SVD(crossprod(FIT$V_history[,,i], FIT$v_hat), nu = 0, nv = 0)$d, 1))))

        u_star_F <- zoo::na.locf(sapply(seq(max_iter), function(i) min(vnorm(FIT$U_history[i,] - u_star),
                                                                       vnorm(FIT$U_history[i,] + u_star))))
        u_hat_F  <- zoo::na.locf(sapply(seq(max_iter), function(i) min(vnorm(FIT$U_history[i,] - FIT$u_hat),
                                                                       vnorm(FIT$U_history[i,] + FIT$u_hat))))

        v_star_F <- zoo::na.locf(sapply(seq(max_iter), function(i) vnorm(tcrossprod(FIT$V_history[,,i]) - tcrossprod(V_star))))
        v_hat_F  <- zoo::na.locf(sapply(seq(max_iter), function(i) vnorm(tcrossprod(FIT$V_history[,,i]) - tcrossprod(V_star))))

        RESULT <- tibble(T = T,
                         p = p,
                         rank = r,
                         d = d,
                         actual_iterations = FIT$iterations,
                         max_iterations = max_iter,
                         iteration = list(seq(max_iter)),
                         U_U_STAR = list(u_star_conv),
                         U_U_HAT  = list(u_hat_conv),
                         V_V_STAR = list(V_star_conv),
                         V_V_HAT  = list(V_hat_conv),
                         U_U_STAR_F = list(u_star_F),
                         U_U_HAT_F  = list(u_hat_F),
                         V_V_STAR_F = list(v_star_F),
                         V_V_HAT_F  = list(v_hat_F),
                         n = n)

        RESULTS <- rbind(RESULTS, RESULT)
        saveRDS(RESULTS, "simulation_computational_convergence_results.rds")
    }
}
