suppressPackageStartupMessages({
    library(tidyverse)
    library(glue)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(ggnewscale)
    library(patchwork)
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
plot(sapply(FITS, function(f) sum(f$X_hat^2) / sum(TENSOR^2)))
R_HAT <- 9

# Visualize u component
M1_G1 <- TBL |>
    select(Hour, Day, HD) |>
    distinct() |>
    mutate(U = FITS[[R_HAT]]$u_hat) |>
    ggplot(aes(x=HD, y=U)) +
    geom_line(group=1) +
    geom_point(aes(color=factor(Hour)), size=3) +
    scale_color_brewer(name="Hour of Study (U Factor)", type="qual", palette=6,
                       labels=1:10)  +
    theme_bw() +
    ylab(expression(hat(u))) +
    theme(legend.position="bottom") +
    xlab("Day of Study") +
    scale_x_discrete(breaks=c("2H08", "3H06", "4H06", "5H06", "6H06"),
                     labels=c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")) +
    geom_vline(alpha=0.5, xintercept=5.5, linetype=2) +
    geom_vline(alpha=0.5, xintercept=14.5, linetype=2) +
    geom_vline(alpha=0.5, xintercept=23.5, linetype=2) +
    geom_vline(alpha=0.5, xintercept=32.5, linetype=2)

### REVISION FIGURES

# Signal identified by u-vector
N_EDGES <- apply(TENSOR, 3, sum) / 2
U_HAT   <- FITS[[R_HAT]]$u_hat

TBL |>
    select(Hour, Day, HD) |>
    distinct() |>
    mutate(u = FITS[[R_HAT]]$u_hat) |>
    mutate(n_edges = N_EDGES) |>
    mutate(edge_density = n_edges / choose(dim(TENSOR)[1], 2)) |>
    ggplot(aes(x = u, y = n_edges, color=factor(Hour))) +
    geom_point(size=3) +
    theme_bw() +
    xlab(expression("SST-PCA Estimated Time Factor" ~ abs(hat(u)))) +
    scale_y_continuous("Number of Student\nInteractions (Edges)",
                       sec.axis = sec_axis(~ . / choose(dim(TENSOR)[1], 2),
                                           name = "Edge Density",
                                           labels=scales::percent)) +
    scale_color_brewer(name="Hour of Study (U Factor)", type="qual", palette=6,
                       labels=1:10) +
    theme(legend.position="bottom") -> SuppB


TENSOR_SCALED <- TENSOR / sqrt(array(rep(apply(TENSOR, 3, sum), each=dim(TENSOR)[1]^2), dim=dim(TENSOR)))

SSTPCA_SCALED <- ss_tpm(TENSOR_SCALED, rank=9, u_init_strategy=0L)

SSTPCA_SCALED$V_hat_w <- trunc_mat(ttv3(TENSOR, SSTPCA_SCALED$u_hat), r)
SSTPCA_SCALED$v_hat_w <- eigen(SSTPCA_SCALED$V_hat_w)$vectors[,seq(r)] * matrix(eigen(SSTPCA_SCALED$V_hat_w)$values[seq(r)], ncol=r, nrow=dim(TENSOR)[1], byrow=TRUE)
SSTPCA_SCALED$X_hat_w <- SSTPCA_SCALED$V_hat_w %o% as.vector(SSTPCA_SCALED$u_hat)

U_HAT_SCALED <- SSTPCA_SCALED$u_hat

TBL |>
    select(Hour, Day, HD) |>
    distinct() |>
    mutate(u = FITS[[R_HAT]]$u_hat,
           u_scaled = U_HAT_SCALED) |>
    ggplot(aes(x = u, y = u_scaled, color=factor(Hour))) +
    geom_point(size=3) +
    theme_bw() +
    scale_color_brewer(name="Hour of Study (U Factor)", type="qual", palette=6,
                       labels=1:10) +
    ggtitle(expression("Estimated Time Factor (" ~ hat(u) ~ ")"))+
    xlab("Original Data") +
    ylab("Normalized Data") -> SuppC


SuppFig_U <- M1_G1 + SuppB + SuppC + plot_annotation(tag_levels="A") +
    plot_layout(guides="collect", nrow=2) & theme(legend.position="bottom", text=element_text(size=24))

ggsave("figure_A2.pdf", SuppFig_U, width=18, height=12)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A2.pdf")
} else {
    system("xdg-open figure_A2.pdf")
}
