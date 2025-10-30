suppressPackageStartupMessages({
    library(tidyverse)
    library(glue)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(ggnewscale)
    library(patchwork)

    ## Used for some methods comparisons
    library(rTensor)
})
set.seed(100)

Rcpp::sourceCpp("../../tensor_factorizations.cpp")

trunc_mat <- function(X, k){
    # Perform a rank-k hard thresholding of a matrix
    sX <- svd(X, nu=k, nv=k)
    sX$u %*% diag(sX$d, k, k) %*% t(sX$v)
}

## BEGIN DATA PREP ##

if(!file.exists("High-School_data_2013.csv.gz")){
    download.file("http://www.sociopatterns.org/wp-content/uploads/2015/07/High-School_data_2013.csv.gz",
                  "High-School_data_2013.csv.gz"
    )
}

# Despite the file name, this is space, not comma, separated data
#
# Read this into R and do a bit of basic clean-up:
#  1) Create E1, E2 as the vertices for each edge (to+from but we symmetrize later)
#  2) Pull out the hour since we're creating one graph per hour of the study
#  3) Create a variable 'HD' which measures the time of the study (day + hour)
TBL <- read.table("High-School_data_2013.csv.gz") |>
    as_tibble() |>
    rename(E1=V2, E2=V3) |>
    mutate(Time=as.POSIXct(V1),
           Hour=lubridate::hour(Time),
           Hour=sprintf("%02d", Hour),
           Day=lubridate::day(Time),
           n_nodes = n_distinct(E1),
           HD = glue('{Day}H{Hour}'))

# Add up the number of edges in each 'time slice' to confirm that
# hourly aggregation gives reasonable density.
TBL_SUM <- TBL |>
    group_by(HD) |>
    summarize(n=n(),
              n_nodes=n_distinct(E1),
              density=n / choose(n_nodes, 2))

# Now we're going to transform this data into a suitable semi-symmetric tensor
# First, we need to find the total number of students in the data set (since
# not all students appear in each time slice) and the total number of time slices
# to pre-allocate our tensor
TIMES <- TBL |> pull(HD) |> unique() |> sort()
N_TIMES <- length(TIMES)

ALL_USERS <- unique(base::union(TBL$E1, TBL$E2)) |> sort()

TBL <- TBL |> mutate(E1ix = match(E1, ALL_USERS),
                     E2ix = match(E2, ALL_USERS),
                     timeix = match(HD, TIMES))
N_USERS <- length(ALL_USERS)

TENSOR <- array(0, dim=c(N_USERS, N_USERS, N_TIMES))


# A bit of fancy indexing to transform this long time series into
# (some of) the edges of our semi-symmetric tensor. Each row of "EDGE_IDS"
# is a tuple (v1, v2, time) which identifies where we want a 1 in our final tensor
EDGE_IDS <- TBL |> select(E1ix, E2ix, timeix) |> as.matrix()
TENSOR[EDGE_IDS] <- 1

for(t in seq(N_TIMES)){
    # Symmetrize each slice separately using an "OR" rule
    TENSOR[,,t] <- TENSOR[,,t] | t(TENSOR[,,t])
}

## END DATA PREP ##


## Estimate a rank based on Gavish & Donoho Optimal Singular Value Thresholding
# See comments in Supplement E.2
TMAT <- array(TENSOR, dim=c(N_USERS^2, N_TIMES))
TMAT_d <- svd(TMAT, nu=0, nv=0)$d
TMAT_KHAT <- sum(TMAT_d > 2.858 * median(TMAT_d))

# Now fit a range of possible SST-PCA ranks to determine the optimal rank
R_MAX <- 20
FITS <- vector("list", R_MAX)
for(r in seq(R_MAX)){
    FITS[[r]] <- ss_tpm(TENSOR, rank=r, u_init_strategy=0L)
    FITS[[r]]$V_hat_w <- trunc_mat(ttv3(TENSOR, FITS[[r]]$u_hat), r)
    FITS[[r]]$v_hat_w <- eigen(FITS[[r]]$V_hat_w)$vectors[,seq(r)] * matrix(eigen(FITS[[r]]$V_hat_w)$values[seq(r)], ncol=r, nrow=dim(TENSOR)[1], byrow=TRUE)
    FITS[[r]]$X_hat_w <- FITS[[r]]$V_hat_w %o% as.vector(FITS[[r]]$u_hat)
}

# Scree plot suggests rank 9
if(interactive()){
    plot(sapply(FITS, function(f) sum(f$X_hat^2) / sum(TENSOR^2)))
}
R_HAT <- 9

# To verify our results, we need to pull in associated metadata about students
# The ID numbers here match the original contact data file
META <- read.table("http://www.sociopatterns.org/wp-content/uploads/2015/09/metadata_2013.txt") |>
    rename(ID=V1, Class=V2, Gender=V3)

ID_IX_MAP <- TBL |> select(ID=E1, IX=E1ix) |> distinct()
ID_IX_MAP <- data.frame(ID=ALL_USERS, IX=seq_along(ALL_USERS))

# Perform K-means on the estimated hat(v) to divide into communities
# Use a large number of starts for stability
K_CLUST <- R_HAT
KMEANS <- kmeans(FITS[[R_HAT]]$v_hat, K_CLUST, nstart=250)

if(interactive()){
    # Optional: Inspect k-means output
    heatmap(FITS[[R_HAT]]$V_hat_w)
    print(KMEANS$size)
    print(KMEANS$cluster)
}


## To assess accuracy of clustering, we

RESULTS <- left_join(TBL,
                     cbind(ID_IX_MAP, ClustID = KMEANS$cluster),
                     by=join_by(E1ix==IX)) |>
    select(-timeix, -ID) |>
    rename(E1_ClustID=ClustID) |>
    left_join(cbind(ID_IX_MAP, ClustID = KMEANS$cluster),
              by=join_by(E2ix==IX)) |>
    select(-ID) |>
    rename(E2_ClustID=ClustID) |>
    select(-V1)

## Examine Cluster Confusions
RESULTS |>
    group_by(E2) |>
    slice_head(n=1) |>
    ungroup() |>
    select(E2_ClustID, V5) -> SST_PCA_Results

with(SST_PCA_Results, mclust::adjustedRandIndex(E2_ClustID, V5))

SST_PCA_Results |> table() |> t()


`%not.in%` <- Negate(`%in%`)
confused_ix <- which(KMEANS$cluster == which.max(KMEANS$size))
confused_id <- ID_IX_MAP |> filter(IX %not.in% confused_ix) |> pull(ID)
# Visualize U-vector: mainly loaded at the same time each day
if(interactive()){
    plot(FITS[[R_HAT]]$u_hat, col= TBL |> select(HD, Day, Hour) |> distinct() |> arrange(HD) |> pull(Hour), pch=16)
}

#### Alternative Approach
SSTPN_A <- abs(FITS[[R_HAT]]$V_hat) > 0.01
diag(SSTPN_A) <- 0
SSTPN_IG <- graph_from_adjacency_matrix(SSTPN_A, mode="undirected")

PN_SCD_membership <- igraph::cluster_leading_eigen(SSTPN_IG)$membership

PN_SCD_DF <- data.frame(IX = 1:N_USERS, PN_SCD_ClusterId = PN_SCD_membership)

left_join(RESULTS, PN_SCD_DF, join_by(E1ix == IX)) |>
    select(E1, V4, E1ix, PN_SCD_ClusterId, E1_ClustID) |>
    group_by(E1) |>
    slice_head(n=1) |>
    ungroup() |>
    select(PN_SCD_ClusterId, V4) -> PN_SCD_RESULTS
PN_SCD_RESULTS |> table()

if(interactive()){
    print(with(PN_SCD_RESULTS, mclust::adjustedRandIndex(PN_SCD_ClusterId, V4)))
}

## Temporal Aggregation
ANYTENSOR <- apply(TENSOR, 1:2, any)
SUMTENSOR <- apply(TENSOR, 1:2, sum)

# spectral comm detect
for(k in 1:40){
    scd_membership <- igraph::cluster_leading_eigen(
        igraph::graph_from_adjacency_matrix(SUMTENSOR >= k, mode="undirected")
    )$membership

    SCD_DF <- data.frame(IX = 1:N_USERS, SCD_ClusterId = scd_membership)

    left_join(RESULTS, SCD_DF, join_by(E1ix == IX)) |>
        select(E1, V4, E1ix, SCD_ClusterId, E1_ClustID) |>
        group_by(E1) |>
        slice_head(n=1) |>
        ungroup() |>
        select(SCD_ClusterId, V4) -> SCD_RESULTS
    SCD_RESULTS |> table()

    if(interactive()){
        print(k)
        print(with(SCD_RESULTS, mclust::adjustedRandIndex(SCD_ClusterId, V4)))
    }
}

## Sub-tensor - re analyze confused_ix
SUBTENSOR <- TENSOR[confused_ix, confused_ix, ]
STMAT <- array(SUBTENSOR, dim=c(length(confused_ix)^2, N_TIMES))
STMAT_d <- svd(STMAT, nu=0, nv=0)$d

STMAT_KHAT <- sum(STMAT_d > 2.858 * median(STMAT_d))
R_MAX <- 20

SFITS <- vector("list", R_MAX)
for(r in seq(R_MAX)){
    SFITS[[r]] <- ss_tpm(SUBTENSOR, rank=r, u_init_strategy=0L)
    SFITS[[r]]$V_hat_w <- trunc_mat(ttv3(SUBTENSOR, FITS[[r]]$u_hat), r)
    SFITS[[r]]$X_hat_w <- SFITS[[r]]$V_hat_w %o% as.vector(SFITS[[r]]$u_hat)
}

if(interactive()){
    plot(sapply(SFITS, function(f) sum(f$X_hat^2) / sum(SUBTENSOR^2)))
}


TENSOR_SCALED <- TENSOR / sqrt(array(rep(apply(TENSOR, 3, sum), each=dim(TENSOR)[1]^2), dim=dim(TENSOR)))

SSTPCA_SCALED <- ss_tpm(TENSOR_SCALED, rank=R_HAT, u_init_strategy=0L)

SSTPCA_SCALED$V_hat_w <- trunc_mat(ttv3(TENSOR, SSTPCA_SCALED$u_hat), r)
SSTPCA_SCALED$v_hat_w <- eigen(SSTPCA_SCALED$V_hat_w)$vectors[,seq(r)] * matrix(eigen(SSTPCA_SCALED$V_hat_w)$values[seq(r)], ncol=r, nrow=dim(TENSOR)[1], byrow=TRUE)
SSTPCA_SCALED$X_hat_w <- SSTPCA_SCALED$V_hat_w %o% as.vector(SSTPCA_SCALED$u_hat)


SCALED_SSTPN_A <- abs(SSTPCA_SCALED$V_hat) > 0.01
diag(SCALED_SSTPN_A) <- 0
SCALED_SSTPN_IG <- graph_from_adjacency_matrix(SCALED_SSTPN_A, mode="undirected")

SCALED_PN_SCD_membership <- igraph::cluster_leading_eigen(SCALED_SSTPN_IG)$membership

SCALED_PN_SCD_DF <- data.frame(IX = 1:N_USERS, SCALED_PN_SCD_ClusterId = SCALED_PN_SCD_membership)

left_join(RESULTS, SCALED_PN_SCD_DF, join_by(E1ix == IX)) |>
    select(E1, V4, E1ix, SCALED_PN_SCD_ClusterId, E1_ClustID) |>
    group_by(E1) |>
    slice_head(n=1) |>
    ungroup() |>
    select(SCALED_PN_SCD_ClusterId, V4) -> SCALED_PN_SCD_RESULTS
SCALED_PN_SCD_RESULTS |> table()

if(interactive()){
    print(with(SCALED_PN_SCD_RESULTS, mclust::adjustedRandIndex(SCALED_PN_SCD_ClusterId, V4)))
}


TENSOR_NONSCALED_SUM <- apply(TENSOR, 1:2, sum)
NONSCALED_SUM_GRAPH <- graph_from_adjacency_matrix(TENSOR_NONSCALED_SUM, mode="undirected")

NONSCALED_SUM_MEMBERSHIP <- igraph::cluster_leading_eigen(NONSCALED_SUM_GRAPH)$membership
NONSCALED_SUM_SCD_DF <- data.frame(IX = 1:N_USERS, NONSCALED_SUM_SCD_ClusterId = NONSCALED_SUM_MEMBERSHIP)

left_join(RESULTS, NONSCALED_SUM_SCD_DF, join_by(E1ix == IX)) |>
    select(E1, V4, E1ix, NONSCALED_SUM_SCD_ClusterId, E1_ClustID) |>
    group_by(E1) |>
    slice_head(n=1) |>
    ungroup() |>
    select(NONSCALED_SUM_SCD_ClusterId, V4) -> NONSCALED_SUM_SCD_RESULTS
NONSCALED_SUM_SCD_RESULTS |> table()

ARI_NONSCALED_SUM <- with(NONSCALED_SUM_SCD_RESULTS, mclust::adjustedRandIndex(NONSCALED_SUM_SCD_ClusterId, V4))



TENSOR_SCALED_SUM <- apply(TENSOR_SCALED, 1:2, sum)
SCALED_SUM_GRAPH <- graph_from_adjacency_matrix(TENSOR_SCALED_SUM, mode="undirected")

SCALED_SUM_MEMBERSHIP <- igraph::cluster_leading_eigen(SCALED_SUM_GRAPH)$membership
SCALED_SUM_SCD_DF <- data.frame(IX = 1:N_USERS, SCALED_SUM_SCD_ClusterId = SCALED_SUM_MEMBERSHIP)

left_join(RESULTS, SCALED_SUM_SCD_DF, join_by(E1ix == IX)) |>
    select(E1, V4, E1ix, SCALED_SUM_SCD_ClusterId, E1_ClustID) |>
    group_by(E1) |>
    slice_head(n=1) |>
    ungroup() |>
    select(SCALED_SUM_SCD_ClusterId, V4) -> SCALED_SUM_SCD_RESULTS
SCALED_SUM_SCD_RESULTS |> table()

ARI_SCALED_SUM <- with(SCALED_SUM_SCD_RESULTS, mclust::adjustedRandIndex(SCALED_SUM_SCD_ClusterId, V4))



## SST-PCA

SSTPN_A <- abs(FITS[[R_HAT]]$V_hat) > 0.01
diag(SSTPN_A) <- 0
## Residual (second moment) Analysis
##
## This will basically match the figure for the first moment analysis we put
## in the paper
##
## The trick to getting this to work however is to 0/1-ify the residual
## graphs before performing the analysis
RESID_TENSOR <- TENSOR - FITS[[R_HAT]]$X_hat
RESID_TENSOR <- RESID_TENSOR > 0.5

RESID_FITS <- vector("list", R_MAX)
for(r in seq(R_MAX)){
    RESID_FITS[[r]] <- ss_tpm(RESID_TENSOR, rank=r, u_init_strategy=0L)
    RESID_FITS[[r]]$V_hat_w <- trunc_mat(ttv3(RESID_TENSOR, RESID_FITS[[r]]$u_hat), r)
    RESID_FITS[[r]]$v_hat_w <- eigen(RESID_FITS[[r]]$V_hat_w)$vectors[,seq(r)] * matrix(eigen(RESID_FITS[[r]]$V_hat_w)$values[seq(r)], ncol=r, nrow=dim(RESID_TENSOR)[1], byrow=TRUE)
    RESID_FITS[[r]]$X_hat_w <- RESID_FITS[[r]]$V_hat_w %o% as.vector(RESID_FITS[[r]]$u_hat)
}
if(interactive()){
    plot(sapply(RESID_FITS, function(f) sum(f$X_hat^2) / sum(RESID_TENSOR^2)))
}
R_HAT <- 10

RESID_GRAPH <- graph_from_adjacency_matrix(abs(RESID_FITS[[R_HAT]]$V_hat_w) > 2e-1, mode="undirected")
RESID_MEMBERSHIP <- igraph::cluster_leading_eigen(RESID_GRAPH)$membership
RESID_SCD_DF <- data.frame(IX = 1:N_USERS, RESID_SCD_ClusterId = RESID_MEMBERSHIP)


## Methods comparison
##
## Here, we compare _many_ different aggregation strategies against SST-PCA
## on the community detection problem associated with this data. Before getting
## into it, we create function that takes a matrix roughly corresponding to
## the SST-PCA principal network and computes the associated ARI.

matrix_to_ARI <- function(A){
    diag(A) <- 0
    A_GRAPH <- graph_from_adjacency_matrix(A, mode="undirected")

    A_MEMBERSHIP <- igraph::cluster_leading_eigen(A_GRAPH)$membership
    A_SCD_DF <- data.frame(IX = 1:N_USERS, A_SCD_ClusterId = A_MEMBERSHIP)

    left_join(RESULTS, A_SCD_DF, join_by(E1ix == IX)) |>
        select(E1, V4, E1ix, A_SCD_ClusterId, E1_ClustID) |>
        group_by(E1) |>
        slice_head(n=1) |>
        ungroup() |>
        select(A_SCD_ClusterId, V4) -> A_SCD_RESULTS

    with(A_SCD_RESULTS, mclust::adjustedRandIndex(A_SCD_ClusterId, V4))
}

## We collect our results in the follwing data frame:
ARI_RESULTS <- data.frame(method=character(),
                          ARI=numeric())


SSTPN_A <- abs(FITS[[R_HAT]]$V_hat) > 0.01
ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="SST-PCA - 1st Factor",
                                ARI = matrix_to_ARI(SSTPN_A)))

RESID_SSTPN_A <- abs(RESID_FITS[[R_HAT]]$V_hat_w) > 0.2
ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="SST-PCA - 2nd Factor",
                                ARI = matrix_to_ARI(RESID_SSTPN_A)))


SUMTENSOR <- apply(TENSOR, 1:2, sum)

ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="ANY",
                                ARI = matrix_to_ARI(SUMTENSOR > 0)))

ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="ALL",
                                ARI = matrix_to_ARI(SUMTENSOR == N_TIMES)))

ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="MULTI",
                                ARI = matrix_to_ARI(SUMTENSOR)))

for(k in seq(1, N_TIMES)){
    ARI_RESULTS <- rbind(ARI_RESULTS,
                         data.frame(method=sprintf("COUNT-%02d", k),
                                    ARI = matrix_to_ARI(SUMTENSOR >= k)))
}

for(k in seq(1, N_TIMES)){
    ARI_RESULTS <- rbind(ARI_RESULTS,
                         data.frame(method=sprintf("SLICE-%02d", k),
                                    ARI = matrix_to_ARI(TENSOR[,,k])))
}

## TWIST
TWIST_P <- function(M, r, delta){
    M_norm <- matrix(sqrt(rowSums(M^2)), nrow=NROW(M), ncol=NCOL(M))
    M_adj <- M * pmin(M_norm, delta) / M_norm
    trunc_mat(M_adj, r)
}


m <- 5
r <- 9

rTENSOR <- as.tensor(TENSOR)

to_mat <- function(x) x@data

# This gets rather nasty with the transposes since K&B notation has an implicit
# transpose and there's no convention on how to organize the output of SVD
# But this seems to work
U <- U0 <- SUMTENSOR |> svd(nu=r) |> _$u |> t()
W <- W0 <- rTENSOR |> ttm(U, 1) |> ttm(U, 2) |> unfold(row_idx=3, col_idx=c(1, 2)) |> to_mat() |> svd(nu=r) |> _$u

for(k in 1:20){
    U_til <- TWIST_P(U, r, delta=0.001)
    W_til <- TWIST_P(W, r, delta=0.001)

    U <- rTENSOR |> ttm(U_til, 2) |> ttm(t(W_til), 3) |> unfold(row_idx=1, col_idx=c(2, 3)) |> to_mat() |> svd(nu=r) |> _$u |> t()
    W <- rTENSOR |> ttm(U_til, 1) |> ttm(U_til, 2) |> unfold(row_idx=3, col_idx=c(1, 2)) |> to_mat() |> svd(nu=r) |> _$u
}

U_TWIST <- t(U)

TWIST_ARI <- data.frame(IX = 1:N_USERS, TWIST_ClusterID = kmeans(U_TWIST, 9)$cluster) |>
    right_join(RESULTS, join_by(IX == E1ix)) |>
    select(E1, V4, IX, TWIST_ClusterID, E1_ClustID) |>
    group_by(E1) |>
    slice_head(n=1) |>
    ungroup() |>
    select(TWIST_ClusterID, V4) |>
    with(mclust::adjustedRandIndex(TWIST_ClusterID, V4))

ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="TWIST" , ARI = TWIST_ARI))


# This works well, but it isn't TWIST per se
ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="TWIST*",
                                # Hand optimized
                                ARI = matrix_to_ARI(crossprod(t(U_TWIST)) > 0.01)))


TWIST_UU <- graph_from_adjacency_matrix(tcrossprod(U_TWIST) > 0.01) |>
    as_tbl_graph() |>
    mutate(ix = row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    geom_edge_link(alpha=0.05) +
    ggtitle("TWIST Principal Network")

## COSIE / MASE
d <- 9
ASE <- function(A){
    eA <- eigen(A, symmetric = TRUE)
    D <- diag(eA$values[1:d])
    V <- eA$vectors[,1:d]
    V %*% sqrt(D)
}
ASEs <- lapply(seq(N_TIMES), function(t) ASE(TENSOR[,,t]))

U <- do.call(cbind, ASEs)
V_hat <- svd(U, nu=d, nv=0)$u

VV = tcrossprod(V_hat)

VV <- VV > 0.004 # Hand-optimized


ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="COSIE",
                                ARI = matrix_to_ARI(VV)))

COSIE_VV <- graph_from_adjacency_matrix(VV) |>
    as_tbl_graph() |>
    mutate(ix = row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    geom_edge_link(alpha=0.05) +
    ggtitle("COSIE Principal Network")


## JEG

## I'm having trouble following the dimensions in Equations (5) and (6)
## of the JEG paper, so I'm just going to derive a descent algorithm directly
## from Equations (1) and (2)

## From Equation (1), to update each lambda_i, the objective has only
## one relevant term |A_i - sum_j lambda[i, j] * h[,j]h[,j]^T|^2
## Take the derivative with respect to lambda[i,j] (a scalar) to get
##
## deriv = -2 h[,j]^T (A_i - sum_{j'} lambda[i, j'] * h[,j']h[,j']^T) h[,j]
##
## Setting deriv = 0, we get
##
## 0 = h[,j]^TA_ih[,j] - sum_{j'} lambda[i, j'] (h[,j']^T h[,j'])^2
## => lambda[i, j] = (h[,j]^TA_ih[,j]- sum_{j' != j} lambda[i, j'] (h[,j]^T h[,j'])^2) / |(h[,j])|^4
##
## since h[,j] are normalized, this is just
## lambda[i, j] = h[,j]^TA_ih[,j]- sum_{j' != j} lambda[i, j'] (h[,j]^T h[,j])^2)
##
## Updating h[,j] is a bit harder, so we'll use an iterative method
##
## The gradient with respect to h[,j] is
## grad = -4 sum_i lambda[i, j] * (A_i - sum_{j'} lambda[i, j'] * h[,j'] * h[,j']^T) * h[,j]

## Per the paper, we'll solve this greedily on the columns of H, starting with LAMBDA = 0

N_CLUST <- 9
H      <- matrix(rnorm(N_USERS * N_CLUST), ncol=N_CLUST)
H      <- H / matrix(sqrt(colSums(H^2)), ncol=N_CLUST, nrow=N_USERS, byrow=TRUE)
LAMBDA <- matrix(0, ncol=N_CLUST, nrow=N_TIMES)

for(j in seq(N_CLUST)){
    for(k in 1:500){
        # Update LAMBDA in closed form
        for(t in seq(N_TIMES)){
            LAMBDA[t, j] = crossprod(H[,j], TENSOR[,,t] %*% H[,j])
            for(j2 in seq(N_CLUST)){
                if(j != j2){
                    LAMBDA[t, j] <- LAMBDA[t, j] - (LAMBDA[t, j2] * crossprod(H[,j], H[,j2]))
                }
            }
        }
        # Take a gradient step on H[,j]

        grad <- matrix(0, ncol=1, nrow=NROW(H))
        # TERM: (A_i - sum_{j'} lambda[i, j'] * h[,j'] * h[,j']^T)
        term <- TENSOR[,,j]
        for(j2 in seq(N_CLUST)){
            term <- term - LAMBDA[t, j2] * tcrossprod(H[, j2])
        }
        for(t in seq(N_TIMES)){
            grad <- grad + LAMBDA[t, j] * term %*% H[,j]
        }

        grad <- -4 * grad

        H[,j] <- H[,j] - 0.01 * grad

        # Renormalize
        H[,j] <- H[,j] / sqrt(sum(H[,j]^2))
    }
}


ARI_RESULTS <- rbind(ARI_RESULTS,
                     data.frame(method="JEG",
                                # Hand optimized threshold
                                ARI = matrix_to_ARI(tcrossprod(H) > 0.012)))

JEG_HH <- graph_from_adjacency_matrix(tcrossprod(H) > 0.012) |>
    as_tbl_graph() |>
    mutate(ix = row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    geom_edge_link(alpha=0.05) +
    ggtitle("JEG Principal Network")

ARI_SSTPCA <- ARI_RESULTS |> filter(method=="SST-PCA - 2nd Factor") |> pull(ARI)

ARI_RESULTS |> mutate(
    category = case_when(
        str_detect(method, "SST-PCA") ~ "SST-PCA",
        str_detect(method, "COUNT") ~ "COUNT",
        str_detect(method, "SLICE") ~ "SLICE",
        str_detect(method, "TWIST") ~ "TWIST",
        TRUE ~ method
    )
) |> mutate(
    method = factor(method, levels=rev(method), ordered=TRUE)
) |>
    ggplot(aes(x = method, fill=category, y=ARI)) +
    geom_col()  +
    theme_bw() +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_brewer(name="Method Type", type="qual", palette=3) +
    geom_hline(yintercept=ARI_SSTPCA, lwd=1, linetype=2) +
    theme(legend.position="bottom") +
    theme(axis.text.y = element_text(angle = 00, vjust = 0.5, hjust=1, size=16),
          axis.text.x = element_text(angle = 00, vjust = 0.5, hjust=1, size=16),
          axis.title = element_text(size=20)) +
    xlab(NULL) +
    ylab("Community Detection Accuracy - Adjusted Rand Index") +
    coord_flip() -> G_ARI

ggsave("figure_A3.pdf", G_ARI, width=12, height=18)


if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A3.pdf")
} else {
    system("xdg-open figure_A3.pdf")
}


as_tbl_graph(SSTPN_IG) |> activate(nodes) |>
    mutate(ix=row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    mutate(CD = KMEANS$cluster) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    scale_color_discrete(name="Class") +
    new_scale_color() +
    geom_edge_link(alpha=0.05)  +
    ggtitle("SST-PCA Principal Network - First Factor") +
    guides(color = guide_legend(override.aes=list(size=2))) -> G2_NO_CD



as_tbl_graph(RESID_GRAPH) |> activate(nodes) |>
    mutate(ix=row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    mutate(CD = KMEANS$cluster) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    scale_color_discrete(name="Class") +
    new_scale_color() +
    geom_edge_link(alpha=0.05)  +
    theme(legend.position="bottom") +
    ggtitle("SST-PCA Principal Network - Second Factor") +
    guides(color = guide_legend(override.aes=list(size=2))) -> M2_G2_NO_CD




G_EMBEDDINGS <- G2_NO_CD + M2_G2_NO_CD +  TWIST_UU + COSIE_VV +JEG_HH +
    plot_layout(guides="collect", ncol = 2) & theme(legend.position="bottom")

ggsave("figure_A4.pdf", G_EMBEDDINGS, width=14, height=14)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A4.pdf")
} else {
    system("xdg-open figure_A4.pdf")
}

