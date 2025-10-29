#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ttv3(const arma::cube& X,
               const arma::vec& u){
    // NOTE: We're not specializing to the semi-symmetric case here
    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword q = X.n_slices;

    if(u.n_elem != q){
        Rcpp::stop("dimension mismatch in ttv3.");
    }

    arma::mat R(n, p, arma::fill::zeros);
    for(arma::uword i = 0; i < n; i++){
        for(arma::uword j = 0; j < p; j++){
            arma::vec X_fiber = X(arma::span(i), arma::span(j), arma::span::all);
            R(i, j) = arma::dot(X_fiber, u);
        }
    }

    return R;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec ttv12(const arma::cube& X,
                const arma::vec& v){
    // NOTE: We're not specializing to the semi-symmetric case here
    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword q = X.n_slices;

    if(v.n_elem != n){
        Rcpp::stop("dimension mismatch in ttv12.");
    }

    if(p != n){
        Rcpp::stop("dimension mismatch in ttv12.");
    }

    arma::vec res(q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; i++){
        res(i) = arma::as_scalar(v.t() * X.slice(i) * v);
    }

    return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec ttm12(const arma::cube& X,
                const arma::mat& v){
    // NOTE: We're not specializing to the semi-symmetric case here
    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword q = X.n_slices;

    if(v.n_rows != n){
        Rcpp::stop("dimension mismatch in ttm12.");
    }

    if(p != n){
        Rcpp::stop("dimension mismatch in ttm12.");
    }

    arma::vec res(q, arma::fill::zeros);

    for(arma::uword i = 0; i < q; i++){
        res(i) = arma::trace(v.t() * X.slice(i) * v);
    }

    return res;
}

#define TOPK_THRESHOLD 1
#define HARD_THRESHOLD 2
#define SOFT_THRESHOLD 3

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec truncated_power_method(const arma::mat& X_raw,
                                 double param,
                                 arma::uword method,
                                 double eps = 1e-6,
                                 bool shift = true) {
    // Implement the Truncated Power Method of Yuan and Zhang (JMLR 2013)
    // Essentially just the classical power method for (leading) eigenvectors
    // with an added truncation step.g
    //
    // Compared with their approach, we allow for possibly indefinite
    // X by actually applying the method to X + omega I for omega sufficiently
    // large to make (X + omega I) positive definite
    //
    // This is kosher because the eigenvectors of X and (X + omega I) are the same
    // because X and (X + omega I) are simultaneously diagonalizable (trivially)

    arma::mat X = X_raw;
    arma::uword p = X.n_rows;

    if(!X.is_sympd()){
        double omega = arma::min(arma::eig_sym(X));
        X += 1.2 * std::abs(omega) * arma::eye(p, p);
    }

    arma::vec v(p, arma::fill::randn); v = arma::normalise(v);
    arma::vec v_old = v;

    arma::uword iter = 0;

    arma::uword k;
    double lambda;

    switch(method){
    case(TOPK_THRESHOLD):{
        k = param;
        if(k == 0){
            Rcpp::stop("k must be at least 1.");
        }

        if(k > p){
            Rcpp::stop("k must be at most the number of rows / columns of X.");
        }
        break;
    }
    default:{
        lambda = param;
        if(lambda < 0){
            Rcpp::stop("lambda must be non-negative.");
        }
        break;
    }
    }

    do{
        v_old = v;
        v = X * v;
        iter++;

        switch(method){
        case(TOPK_THRESHOLD):{
            // TOP-K (by absolute value) TRUNCATION
            arma::uvec sorted_ix = arma::sort_index(arma::abs(v), "descend");
            for(arma::uword i = k; i < p; i++){
                v(sorted_ix(i)) = 0;
            }
            break;
        }
        case(HARD_THRESHOLD):{
            for(arma::uword i = 0; i < p; i++){
            v(i) = std::abs(v(i)) > lambda ? v(i) : 0;
        }
            break;
        }
        case(SOFT_THRESHOLD):{
            for(arma::uword i = 0; i < p; i++){
            v(i) = std::abs(v(i)) > lambda ? (arma::sign(v(i)) * (std::abs(v(i)) - lambda)) : 0;
        }
            break;
        }
        default:{
            Rcpp::stop("Unknown penalty method!");
        }
        }

        if(arma::norm(v) == 0){
            Rcpp::stop("v_hat is zero - that's bad!");
        }
        v = arma::normalise(v);

        if(iter > 1000){
            break;
        }

    } while (arma::norm(v - v_old) > eps);
    return v;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ss_tensor_power_method_large(const arma::cube& X,
                                        const arma::uword u_init_strategy = 0,
                                        const arma::uword rank = 1,
                                        const double eps = 1e-6,
                                        const arma::vec& u_init = 0){
    // Given a p x p x T semi-symmetric tensor X
    // (i.e., X(:, :, k) is a symmetric matrix for all k)
    // calculate a sparse CP-type decomposition which respects the
    // structure of X while inducing sparsity in the symmetric components
    //
    // That is, solve
    //   arg_max{v in B_p, u in B_T} <X, v o v o u> - lambda * P(v)
    // where
    //   - B_p, B_T are the (Euclidean) unit spheres of dimension p, T respectively
    //   - < > is the Frobenius type inner product
    //   - o denotes the tensor outer product
    //   - P can be the L0 or L1 penalty. [We also include a top K thresholding option, but that's not a normal penalty]
    //
    // Even though this is CP like, the non-convexity of the problem means that
    // we need to explicitly state the algorithm used: normal CP (i.e., CP-ALS)
    // won't preserve semi-symmetry so we use a modified version of the
    // Tensor Power Method proposed in "Regularized Tensor Factorizations and Higher-Order
    // Principal Components Analysis" (GI Allen, ArXiv 1202.2476; unpublished)
    // which consists of the following iterates, repeated to convergence:
    //
    // u_hat = C x1 v_hat x2 v_hat / |C x1 v_hat x2 v_hat|
    // v_hat = arg_max{v in B_p} v^T (C x3 u_hat) v - lambda P(v)
    //
    // where x1, x2, x3 are the (contractive) tensor-vector products along
    // the first, second, third mode respectively
    //
    // To solve the second problem, we use the truncated power method
    // of Yuan and Zhang implemented above
    //
    // FIXME: The penalized version is not presently implemented!
    //
    // The implementation is pretty easy, but armadillo's tensor (cube) type
    // is a bit limited, so we have to do book-keeping manually
    // using the ttv12(C, v_hat) -> C x1 v_hat x2 v_hat and
    //           ttv3(C, u_hat) -> C x3 u_hat
    // helper functions

    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword q = X.n_slices;

    if(rank != 1){
        Rcpp::stop("This function only supports rank-1 v.");
    }

    if(n != p){
        Rcpp::stop("X needs to be semi-symmetric.");
        // NOTE: We're not actually checking semi-symmetry: only that the dimensions work
    }

    arma::vec u(q);

    if(u_init_strategy == 0){
        u.fill(1.0);
    } else if(u_init_strategy == 1){
        u.randn();
    } else if(u_init_strategy == 2){
        u = u_init;
    } else {
        Rcpp::stop("Unknown u-initialization scheme!");
    }

    u = arma::normalise(u);

    arma::vec v(n, arma::fill::zeros);

    arma::vec u_old;
    arma::vec v_old;

    arma::uword iter = 0;

    do{
        u_old = u;
        v_old = v;
        iter++;

        const arma::mat Xu = ttv3(X, u);

        v = truncated_power_method(Xu, 0, 3); // Lasso w/0 thresholding = normal power method

        u = ttv12(X, v);
        u = arma::normalise(u);

        if(iter > 1000){
            break;
        }

    } while (arma::norm(v - v_old)  + arma::norm(u - u_old) > eps);

    return Rcpp::List::create(Rcpp::Named("u_hat") = u,
                              Rcpp::Named("v_hat") = v,
                              Rcpp::Named("iterations") = iter);

}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ss_tensor_power_method(const arma::cube& X,
                                  const arma::uword u_init_strategy = 0,
                                  const arma::uword rank = 1,
                                  const double eps = 1e-6,
                                  const arma::vec& u_init = 0,
                                  const arma::uword max_iter=1000){
    // Given a p x p x T semi-symmetric tensor X
    // (i.e., X(:, :, k) is a symmetric matrix for all k)
    // calculate a sparse CP-type decomposition which respects the
    // structure of X while inducing sparsity in the symmetric components
    //
    // That is, solve
    //   arg_max{v in B_p, u in B_T} <X, v o v o u> - lambda * P(v)
    // where
    //   - B_p, B_T are the (Euclidean) unit spheres of dimension p, T respectively
    //   - < > is the Frobenius type inner product
    //   - o denotes the tensor outer product
    //   - P can be the L0 or L1 penalty. [We also include a top K thresholding option, but that's not a normal penalty]
    //
    // Even though this is CP like, the non-convexity of the problem means that
    // we need to explicitly state the algorithm used: normal CP (i.e., CP-ALS)
    // won't preserve semi-symmetry so we use a modified version of the
    // Tensor Power Method proposed in "Regularized Tensor Factorizations and Higher-Order
    // Principal Components Analysis" (GI Allen, ArXiv 1202.2476; unpublished)
    // which consists of the following iterates, repeated to convergence:
    //
    // u_hat = C x1 v_hat x2 v_hat / |C x1 v_hat x2 v_hat|
    // v_hat = arg_max{v in B_p} v^T (C x3 u_hat) v - lambda P(v)
    //
    // where x1, x2, x3 are the (contractive) tensor-vector products along
    // the first, second, third mode respectively
    //
    // To solve the second problem, we use the truncated power method
    // of Yuan and Zhang implemented above
    //
    // FIXME: The penalized version is not presently implemented!
    //
    // The implementation is pretty easy, but armadillo's tensor (cube) type
    // is a bit limited, so we have to do book-keeping manually
    // using the ttv12(C, v_hat) -> C x1 v_hat x2 v_hat and
    //           ttv3(C, u_hat) -> C x3 u_hat
    // helper functions

    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword q = X.n_slices;

    if(n != p){
        Rcpp::stop("X needs to be semi-symmetric.");
        // NOTE: We're not actually checking semi-symmetry: only that the dimensions work
    }

    arma::vec u(q);

    if(u_init_strategy == 0){
        u.fill(1.0);
    } else if(u_init_strategy == 1){
        u.randn();
    } else if(u_init_strategy == 2){
        u = u_init;
    } else {
        Rcpp::stop("Unknown u-initialization scheme!");
    }

    u = arma::normalise(u);

    arma::mat v(n, rank, arma::fill::zeros);

    arma::vec u_old;
    arma::mat v_old;

    arma::uword iter = 0;

    do{
        u_old = u;
        v_old = v;
        iter++;

        const arma::mat Xu = ttv3(X, u);
        
        // Deal with small asymmetries due to floating point
        const arma::mat Xu_sym = (Xu + Xu.t()) / 2; 

        arma::mat eigenvectors_Xu;
        arma::vec eigenvalues_Xu;

        eig_sym(eigenvalues_Xu, eigenvectors_Xu, Xu_sym);

        v = eigenvectors_Xu.tail_cols(rank);

        u = ttm12(X, v);
        u = arma::normalise(u);
        
        if(iter % 20 == 0){
          Rcpp::checkUserInterrupt();
        }

        if(iter > max_iter){
            break;
        }

    } while (arma::norm(v - v_old)  + arma::norm(u - u_old) > eps);

    return Rcpp::List::create(Rcpp::Named("u_hat") = u,
                              Rcpp::Named("v_hat") = v,
                              Rcpp::Named("iterations") = iter);

}

/*** R
ss_tpm <- function(X, ...){
    fit <- ss_tensor_power_method(X, ...)
    fit$X_hat <- tcrossprod(fit$v_hat) %o% as.vector(fit$u_hat)
    fit$u_hat <- as.vector(fit$u_hat)
    fit$V_hat <- tcrossprod(fit$v_hat)
    fit$v_hat <- fit$v_hat
    fit$d <- sum(fit$X_hat * X) / sum(fit$X_hat^2)
    fit$X_hat <- fit$X_hat * fit$d
    rownames(fit$V_hat) <- colnames(fit$V_hat) <- rownames(X[,,1])

    fit
}

ss_tpm_large <- function(X, ...){
    fit <- ss_tensor_power_method_large(X, ...)
    fit$X_hat <- tcrossprod(fit$v_hat) %o% as.vector(fit$u_hat)
    fit$u_hat <- as.vector(fit$u_hat)
    fit$V_hat <- tcrossprod(fit$v_hat)
    fit$v_hat <- fit$v_hat
    fit$d <- sum(fit$X_hat * X) / sum(fit$X_hat^2)
    fit$X_hat <- fit$X_hat * fit$d
    rownames(fit$V_hat) <- colnames(fit$V_hat) <- rownames(X[,,1])

    fit
}
*/


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube tensor_cusum_cpp(const arma::cube& X){
    // Given a p x p x T semi-symmetric tensor X
    // (i.e., X(:, :, k) is a symmetric matrix for all k)
    // calculate a CUSUM chart elementwise (across the final mode T)
    //
    // There seem to be (many!) alternate definitions in the literature, so here
    // we follow the convention given by Wang and Samworth (2018) whereby
    //
    // C_{ijt} = sqrt(T / (t * (T - t))) * (t / T * sum_{r = 1}^T X_{ijr} - sum_{r = 1}^t X_{ijr})
    //
    // This is pretty easy to implement: we use a bit of armadillo bradcoasting,
    // cumulative sums and cube (tensor) subassignment to avoid redundant work,
    // but it's not necessary

    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    const arma::uword T = X.n_slices;

    if(n != p){
        Rcpp::stop("X needs to be semi-symmetric.");
        // NOTE: We're not actually checking semi-symmetry: only that the dimensions work
    }

    arma::cube C(n, p, T - 1, arma::fill::zeros);
    arma::vec t = arma::regspace(1, T - 1);

    for(arma::uword i = 0; i < n; i++){
        for(arma::uword j = 0; j < p; j++){
            arma::vec fiber = X(arma::span(i), arma::span(j), arma::span::all);
            double sum = arma::accu(fiber);
            arma::vec cumsum = arma::cumsum(fiber);
            cumsum = cumsum.head(T - 1);

            arma::vec cusum_stat = arma::sqrt(T / (t % (T - t))) % (t / T * sum - cumsum);

            C(arma::span(i), arma::span(j), arma::span::all) = cusum_stat;
        }
    }

    return C;
}

/*** R
tensor_cusum <- function(X){
    C <- tensor_cusum_cpp(X)
    dimnames(C)[[1]] <- dimnames(C)[[2]] <- colnames(X[,,1])
    C
}
*/
