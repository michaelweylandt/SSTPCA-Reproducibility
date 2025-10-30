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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ss_tensor_power_method(const arma::cube& X,
                                  const arma::vec& u_init = 0,
                                  const arma::uword rank = 1,
                                  const double eps = 1e-6,
                                  const arma::uword max_iter = 50
                                  ){
    // Given a p x p x T semi-symmetric tensor X
    // (i.e., X(:, :, k) is a symmetric matrix for all k)
    // calculate a SS-TPCA decomposition which respects the structure of X
    // structure of X
    //
    // This file implements a "tracked" SS-TPCA which returns the u and v iterates
    //
    // The implementation is pretty easy, but armadillo's tensor (cube) type
    // is a bit limited, so we have to do book-keeping manually
    // using the ttv12(C, v_hat) -> C x1 v_hat x2 v_hat and
    //           ttv3(C, u_hat) -> C x3 u_hat
    // helper functions

    const arma::uword p = X.n_rows;
    const arma::uword p1 = X.n_cols;
    const arma::uword T = X.n_slices;

    if(p != p1){
        Rcpp::stop("X needs to be semi-symmetric.");
        // NOTE: We're not actually checking semi-symmetry: only that the dimensions work
    }

    if(T != u_init.n_elem){
        Rcpp::stop("Initialization vector U wrong length");
    }

    arma::vec u = u_init;
    u = arma::normalise(u);

    arma::mat V(p, rank, arma::fill::zeros);

    arma::vec u_old;
    arma::mat V_old;

    arma::uword iter = 0;

    arma::mat  u_storage(max_iter, T, arma::fill::value(arma::datum::nan));
    arma::cube V_storage(p, rank, max_iter, arma::fill::value(arma::datum::nan));

    do{
        u_old = u;
        V_old = V;

        iter++;


        const arma::mat Xu = ttv3(X, u);

        arma::mat eigenvectors_Xu;
        arma::vec eigenvalues_Xu;

        eig_sym(eigenvalues_Xu, eigenvectors_Xu, Xu);

        V = eigenvectors_Xu.tail_cols(rank);

        u = ttm12(X, V);
        u = arma::normalise(u);

        u_storage.row(iter - 1) = u.t();
        V_storage.slice(iter - 1) = V;

        if(iter >= max_iter){
            break;
        }

    } while (arma::norm(V - V_old)  + arma::norm(u - u_old) > eps);

    return Rcpp::List::create(Rcpp::Named("u_hat") = u,
                              Rcpp::Named("v_hat") = V,
                              Rcpp::Named("iterations") = iter,
                              Rcpp::Named("max_iter") = max_iter,
                              Rcpp::Named("rank") = rank,
                              // In order to facilitate comparisons among runs, let's just keep the full "max_iter" history
                              // even though the tails will be dull for most iterations
                              //
                              // TODO: Fill forward on the R side - see commented out code below for truncation steps
                              //
                              Rcpp::Named("U_history") = u_storage,  //u_storage.rows(0, iter),
                              Rcpp::Named("V_history") = V_storage); //V_storage.slices(0, iter));

}

/*** R
ss_tpm <- function(X, ...){
    fit <- ss_tensor_power_method(X, ...)
    fit$X_hat <- tcrossprod(fit$v_hat) %o% as.vector(fit$u_hat)
    fit$u_hat <- as.vector(fit$u_hat)
    fit$V_hat <- tcrossprod(fit$v_hat)
    fit$v_hat <- fit$v_hat
    fit$U_history <- fit$U_history
    fit$V_history <- array(fit$V_history, c(NROW(fit$v_hat), fit$rank, fit$max_iter))
    fit$d <- sum(fit$X_hat * X) / sum(fit$X_hat^2)
    fit$X_hat <- fit$X_hat * fit$d
    rownames(fit$V_hat) <- colnames(fit$V_hat) <- rownames(X[,,1])

    fit
}

*/
