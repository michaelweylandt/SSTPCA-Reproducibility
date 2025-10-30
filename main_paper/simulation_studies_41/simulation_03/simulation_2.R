library(tidyverse)

Rcpp::sourceCpp("../../../tensor_factorizations.cpp")

vnorm <- function(x) sqrt(sum(x^2))
mse <- function(x) vnorm(x) / length(x)

n <- 50

RESULTS <- data.frame()
for(T in seq(10, 110, by = 10)){
    for(p in seq(10, 110, by = 10)){
        for(d in seq(1, 20)){
            results <- do.call(rbind, replicate(n, {
                u <- rnorm(T); u <- u / vnorm(u)
                v <- rnorm(p); v <- v / vnorm(v)

                ss_noise <- function(T, p){
                    replicate(T, simplify = "array", {
                        E <- matrix(rnorm(p * p, sd = 1/sqrt(2)), ncol = p, nrow = p)
                        E <- (E + t(E)) / 2
                    })
                }

                X     <- d * (v %o% v %o% u)
                X_hat <- X + ss_noise(T, p)

                fit_stable_init <- ss_tpm(X_hat, u_init_strategy = 0L)
                fit_randg_init  <- ss_tpm(X_hat, u_init_strategy = 1L)
                fit_oracle_init <- ss_tpm(X_hat, u_init_strategy = 2L, u_init = u)

                make_results <- function(fit, label){
                    data.frame(T = T,
                               p = p,
                               d = d,
                               u_svd = norm(fit$u_hat %o% fit$u_hat - u %o% u, type = "2"),
                               u_acos = acos(abs(sum(u * fit$u_hat))) * 180 / pi,
                               v_svd = norm(fit$V_hat - v %o% v, type = "2"),
                               v_acos = acos(abs(sum(v * fit$v_hat))) * 180 / pi,
                               u_mse = mse(fit$u_hat - u),
                               v_mse = mse(fit$V_hat - v %o% v),
                               X_mse = mse(fit$X_hat - X),
                               X_err_scaled = vnorm(fit$X_hat - X) / vnorm(X),
                               label = label)
                }

                rbind(make_results(fit_stable_init, label = "stable"),
                      make_results(fit_randg_init, label = "random"),
                      make_results(fit_oracle_init, label = "oracle"))
            }, simplify = FALSE))

            RESULTS <- rbind(RESULTS, results |> group_by(label) |> summarize(across(everything(), mean)))
            saveRDS(RESULTS, "simulation_2.rds")
        }
    }
}
