library(tidyverse)
library(rTensor)
library(igraph)
library(Rcpp)

Rcpp::sourceCpp("sstpca.cpp")

vnorm <- function(x) sqrt(sum(x^2))
rmse <- function(x) vnorm(x) / length(x)
SVD <- function(x, ...){if(anyNA(x)) list(d = NA) else svd(x, ...)}

truncate_matrix <- function(X, rank = Matrix::rankMatrix(X)){
    s <- svd(X, nu = rank, nv = rank)
    s$u %*% diag(s$d, rank, rank) %*% t(s$v)
}

RESULTS <- tibble(T = integer(),
                  p = integer(),
                  rank = integer(),
                  VV_L2 = double(),
                  method = double(),
                  n = integer(),
                  model = character())

### DEFAULT GENERATIVE MODEL PARAMETERS FOR NOW
T <- 20; p <- 200; k <- 5; REPS <- 50
sbm_p <- 0.8; sbm_q <- 0.2

for(T in c(20, 40, 60, 80, 100)){
    for(k in c(3, 5, 7)){
        for(p in 105 * c(1, 2, 3, 4, 5)){
            for(n in seq(REPS)){
                make_sbm <- function() {
                    g <- sample_sbm(p, toeplitz(c(sbm_p, rep(sbm_q, k - 1))), block.sizes = rep(p / k, k))
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }
                G_star <- Reduce(`+`, replicate(500, make_sbm(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_sbm(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="SBM") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }

            for(n in seq(REPS)){
                positions <- sample_dirichlet(p, alpha = rep(0.3, k))
                make_rdpg <- function(){
                    g <- sample_dot_product(positions)
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }

                G_star <- Reduce(`+`, replicate(500, make_rdpg(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_rdpg(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="Dirichlet-RDPG") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }

            for(n in seq(REPS)){
                positions <- sample_sphere_surface(dim = 5, n = p, positive = TRUE)
                make_rdpg <- function(){
                    g <- sample_dot_product(positions)
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }

                G_star <- Reduce(`+`, replicate(500, make_rdpg(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_rdpg(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="Positive-Sphere-RDPG") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }

            for(n in seq(REPS)){
                positions <- sample_sphere_surface(dim = 5, n = p)
                make_rdpg <- function(){
                    g <- sample_dot_product(positions)
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }

                G_star <- Reduce(`+`, replicate(500, make_rdpg(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_rdpg(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="Sphere-RDPG") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }


            for(n in seq(REPS)){
                positions <- sample_sphere_volume(dim = 5, n = p, positive = TRUE)
                make_rdpg <- function(){
                    g <- sample_dot_product(positions)
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }

                G_star <- Reduce(`+`, replicate(500, make_rdpg(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_rdpg(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="Positive-Ball-RDPG") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }

            for(n in seq(REPS)){
                positions <- sample_sphere_volume(dim = 5, n = p)
                make_rdpg <- function(){
                    g <- sample_dot_product(positions)
                    g <- as_edgelist(g)
                    M <- matrix(0, p, p); M[as.matrix(g)] <- 1;
                    M <- M + t(M)
                    M
                }

                G_star <- Reduce(`+`, replicate(500, make_rdpg(), simplify = FALSE))
                V_star <- eigen(truncate_matrix(G_star, k))$vectors[,1:k]

                X <- replicate(T, make_rdpg(), simplify = FALSE); X <- simplify2array(X)

                X_vec <- matrix(X, nrow = T, ncol = p^2, byrow = TRUE)
                X_rTensor <- as.tensor(X) # tensor representation used by rTensor package

                SS_TPCA_FIT <- ss_tpm(X, u_init = rep(1, T), rank = k)
                PCA_FIT     <- svd(X_vec, nu = 1, nv = 1)
                CP_FIT_1    <- cp(X_rTensor, 1, max_iter = 50)
                CP_FIT_k    <- cp(X_rTensor, k, max_iter = 50)
                HOSVD_FIT   <- hosvd(X_rTensor, c(k, k, 1))      # HOSVD
                HOOI_FIT    <- tucker(X_rTensor, c(k, k, 1), 50) # HOOI

                VV_STAR      <- tcrossprod(V_star)
                PCA_VV       <- matrix(PCA_FIT$v, nrow = p, ncol = p)
                TRUNC_PCA_VV <- truncate_matrix(matrix(PCA_FIT$v, nrow = p, ncol = p), k)
                SSTPCA_VV    <- tcrossprod(SS_TPCA_FIT$v_hat)

                CP1_VV       <- CP_FIT_1$lambda[1] * tcrossprod(CP_FIT_1$U[[1]], CP_FIT_1$U[[2]])
                CP1_VV       <- CP1_VV + t(CP1_VV)
                CP1_VV       <- CP1_VV / norm(CP1_VV, "2")

                CPK_VV       <- matrix(0, p, p)
                for(i in seq(k)){
                    temp   <- CP_FIT_k$lambda[i] * tcrossprod(CP_FIT_k$U[[1]][,k], CP_FIT_k$U[[2]][,k])
                    temp   <- temp + t(temp)
                    temp   <- temp / norm(temp, "2")
                    CPK_VV <- CPK_VV + temp
                }

                HOSVD_VV     <- tcrossprod(HOSVD_FIT$U[[1]], HOSVD_FIT$U[[2]]) # These two u components are the same, so not strictly necessary
                HOOI_VV      <- tcrossprod(HOOI_FIT$U[[1]], HOOI_FIT$U[[1]]) # These two u components should be the same, but not enforced as such by estimator
                # so arbitrarily selecting #1

                tibble(method = c("PCA", "Truncated-PCA", "SS-TPCA", "CP-1", "CP-K", "HOSVD", "HOOI"),
                       VV_L2  = c(vnorm(VV_STAR - PCA_VV),
                                  vnorm(VV_STAR - TRUNC_PCA_VV),
                                  vnorm(VV_STAR - SSTPCA_VV),
                                  vnorm(VV_STAR - CP1_VV),
                                  vnorm(VV_STAR - CPK_VV),
                                  vnorm(VV_STAR - HOSVD_VV),
                                  vnorm(VV_STAR - HOOI_VV)),
                       T = T,
                       p = p,
                       rank = k,
                       n = n,
                       model="Ball-RDPG") -> RESULT

                RESULTS <- rbind(RESULTS, RESULT)
            }

            saveRDS(RESULTS, "simulation_methods_comparison_results")
        }
    }
}
