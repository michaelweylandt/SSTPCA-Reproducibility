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
Rcpp::List ss_tensor_power_method_simul(const arma::cube& X,
                                        const arma::uword n_factors,
                                        const arma::vec& ranks,
                                        const arma::uword u_init_strategy = 1,
                                        const double eps = 1e-6,
                                        const arma::mat& u_init = 0,
                                        const arma::uword max_iter=1000,
                                        const bool ortho_q = false
                                          ){
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

    arma::mat U(q, n_factors);

    if(u_init_strategy == 0){
        U.fill(1.0);
    } else if(u_init_strategy == 1){
        U.randn();
    } else if(u_init_strategy == 2){
        U = u_init;
    } else {
        Rcpp::stop("Unknown u-initialization scheme!");
    }

    U = arma::normalise(U);

    arma::field<arma::mat> Vs(n_factors);
      
    for(arma::uword f_ix=0; f_ix < n_factors; f_ix++){
      arma::mat V_shell(n, ranks(f_ix), arma::fill::zeros);
      Vs(f_ix) = V_shell;
    }
      
    arma::mat U_old;
    arma::uword iter = 0;

    do{
        U_old = U;
        iter++;
        
        for(arma::uword f_ix = 0; f_ix < n_factors; f_ix++){
          arma::vec u = U.col(f_ix); 
          arma::uword rank_f = ranks(f_ix); 
          const arma::mat Xu = ttv3(X, u);
        
          // Deal with small asymmetries due to floating point
          const arma::mat Xu_sym = (Xu + Xu.t()) / 2; 

          arma::mat eigenvectors_Xu;
          arma::vec eigenvalues_Xu;
  
          eig_sym(eigenvalues_Xu, eigenvectors_Xu, Xu_sym);

          arma::mat v = eigenvectors_Xu.tail_cols(rank_f);
          Vs(f_ix) = v; 
          
          u = ttm12(X, v);
          U.col(f_ix) = arma::normalise(u); 
        }
        
        if(ortho_q){
          arma::mat Q; 
          arma::mat R; 
          arma::qr_econ(Q, R, U);
          U = Q;
        }
        
        if(iter % 20 == 0){
          Rcpp::checkUserInterrupt();
        }

        if(iter > max_iter){
            break;
        }

    } while (arma::norm(U - U_old) > eps);
    
    if(arma::rank(U) < n_factors){
      Rcpp::Rcout << "Possible redundant factors found - check fit carefully" << std::endl;
    }

    return Rcpp::List::create(Rcpp::Named("u_hat") = U,
                              Rcpp::Named("v_hat") = Vs,
                              Rcpp::Named("iterations") = iter,
                              Rcpp::Named("n_factors")=n_factors,
                              Rcpp::Named("ranks")=ranks);

}

/*** R
ss_tpm_simul <- function(X, ...){
    fit <- ss_tensor_power_method_simul(X, ...)
    # fit$u_hat is fine
    # Convert fit$v_hat from arma::field (essentially a list with dimensions) 
    # to a regular list
    fit$v_hat <- apply(fit$v_hat, 2, identity)[[1]]
    
    ## Compute baseline X PaRMs
    n_factors <- fit$n_factors
    fit$X_hat <- vector("list", n_factors)
    fit$d <- numeric(n_factors)
    for(f in seq(n_factors)){
      fit$X_hat[[f]] <- tcrossprod(fit$v_hat[[f]]) %o% as.vector(fit$u_hat[,f])
    }
    # Least squares fit given our available vectors
    m <- lm(as.vector(X) ~ 0 + do.call(cbind, lapply(fit$X_hat, as.vector)))
    fit$resid <- sum(resid(m)^2)
    fit$d <- as.vector(coef(m))
    for(f in seq(n_factors)){
      fit$X_hat[[f]] <- fit$d[f] * fit$X_hat[[f]]
    }
    fit
}
*/

