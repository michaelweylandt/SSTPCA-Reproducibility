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

if(FALSE){
    # Optional: Inspect k-means output
    heatmap(FITS[[R_HAT]]$V_hat_w)
    print(KMEANS$size)
    print(KMEANS$cluster)
}

# Figure A1 -- this is basically the same code as Figure 6, but
# here applied directly to `TENSOR` instead of `RESID_TENSOR` so
# that we're looking at structure in the first moment ("low-rank approximation")

# Truncate out small edges and self-loops for legibility
M1_SSTPN_A <- abs(FITS[[R_HAT]]$V_hat) > 0.01
diag(M1_SSTPN_A) <- 0
M1_SSTPN_IG <- graph_from_adjacency_matrix(M1_SSTPN_A, mode="undirected")

# Visualize SST-PCA estimated principal network
as_tbl_graph(M1_SSTPN_IG) |> activate(nodes) |>
    mutate(ix=row_number()) |>
    left_join(ID_IX_MAP, join_by(ix == IX)) |>
    mutate(CD = KMEANS$cluster) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    scale_color_discrete(name="Student Class (Outer)") +
    new_scale_color() +
    geom_node_point(aes(color=factor(CD)), size=1) +
    scale_color_brewer(name="Estimated Community (Inner)", type="qual", palette=3) +
    geom_edge_link(alpha=0.05)  +
    theme(legend.position="bottom") +
    #ggtitle("SST-PCA Principal Network") +
    guides(color = guide_legend(override.aes=list(size=2))) -> M1_G2

# Compare to temporal "ANY" aggregation, followed by spectral community detection
ANYTENSOR <- apply(TENSOR, 1:2, any)
ANYTENSOR_IG <- igraph::graph_from_adjacency_matrix(ANYTENSOR, mode="undirected")

as_tbl_graph(ANYTENSOR_IG) |> activate(nodes) |>
    mutate(ix = row_number()) |>
    left_join(ID_IX_MAP, join_by(ix==IX)) |>
    mutate(CD=cluster_leading_eigen(ANYTENSOR_IG)$membership) |>
    left_join(META, join_by(ID == ID)) |>
    ggraph(layout="kk") +
    geom_node_point(aes(color=Class), size=2.5) +
    scale_color_discrete(name="Student Class (Outer)") +
    new_scale_color() +
    geom_node_point(aes(color=factor(CD)), size=1) +
    scale_color_brewer(name="Estimated Community (Inner)", type="qual", palette=3) +
    geom_edge_link(alpha=0.05)  +
    theme(legend.position="bottom") +
    #ggtitle("Temporally Aggregated Graph") +
    guides(color = guide_legend(override.aes=list(size=2)))-> M1_G3

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

# Scree plot for ranks
data.frame(rerr=sapply(FITS, function(f) sum(f$X_hat^2) / sum(TENSOR^2)),
           rank=seq(R_MAX)) |>
    ggplot(aes(x=rank, y=rerr)) + geom_point(size=2) + geom_line() + xlab("Rank of V Factor") +
    ylab("%VE") +
    theme_bw() +
    geom_vline(xintercept=R_HAT, linetype=2, linewidth=2, color="red4") +
    scale_y_continuous(labels=scales::percent)-> M1_G4

library(patchwork)
M1_G_COMP <- (M1_G2 + M1_G3 + plot_layout(guides="collect") & theme(legend.position="bottom") &
                  theme(legend.text = element_text(size = 14),
                        legend.title=element_text(size=16)))+
    (M1_G1 / M1_G4 + plot_layout(guides="keep") & theme(axis.text = element_text(size=12),
    )) +
    plot_annotation(tag_levels="A")

ggsave("figure_A1.pdf", M1_G_COMP, width=18, height=5)
browseURL("figure_A1.pdf")
