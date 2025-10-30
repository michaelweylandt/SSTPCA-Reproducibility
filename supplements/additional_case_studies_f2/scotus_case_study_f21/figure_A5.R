library(Rcpp)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(patchwork)
library(reshape2)

Rcpp::sourceCpp("../../../tensor_factorizations.cpp")


# Helper functions for visualization
visualize_court <- function(M, threshold = 0.8, true_names = NULL){
    base_vector <- c("CJ", paste0("AJ", 1:8))
    if(!is.null(true_names)){
        rownames(M) <- colnames(M) <- paste(colnames(M), true_names, sep = " - ")
        base_vector <- paste(base_vector, true_names, sep = " - ")
    }


    as.data.frame(M) |> rownames_to_column("y") |>
        reshape2::melt(variable.name = "x") |>
        mutate(pct_value = paste0(round(100 * value, 0), "%"),
               x = factor(x, levels = base_vector),
               y = factor(y, levels = base_vector)) |>
        filter(x != y) |>
        ggplot(aes(x = x, y = y, fill = value)) + geom_tile() + theme_bw() +
        xlab(NULL) + ylab(NULL) +
        (if(min(M) < 0){
            scale_fill_gradient2(low = "red", mid = "white", high = "blue")
        } else {
            scale_fill_gradient(low = "white", high = "red")
        }) +
        guides(fill = "none") + geom_text(aes(label = pct_value)) -> G1


    if(!is.null(true_names)){
        G1 <- G1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }


    G2 <- M |>
        as_tbl_graph() |>
        activate(edges) |>
        filter(abs(weight) < 1,  abs(weight) > threshold) |>
        ggraph() +
            geom_node_point(color = "blue", size = 2) +
            geom_node_text(aes(label = name), repel = TRUE) +
            geom_edge_link(aes(width = weight), alpha = 0.5) +
            scale_edge_width(range = c(0, 1), guide = "none") +
            theme_graph()


    list(G1, G2)


}


visualize_court_no_abs <- function(M, threshold = 0.8, true_names = NULL){
    base_vector <- c("CJ", paste0("AJ", 1:8))
    if(!is.null(true_names)){
        rownames(M) <- colnames(M) <- paste(colnames(M), true_names, sep = " - ")
        base_vector <- paste(base_vector, true_names, sep = " - ")
    }


    as.data.frame(M) |> rownames_to_column("y") |>
        reshape2::melt(variable.name = "x") |>
        mutate(pct_value = paste0(round(100 * value, 0), "%"),
               x = factor(x, levels = base_vector),
               y = factor(y, levels = base_vector)) |>
        filter(x != y) |>
        ggplot(aes(x = x, y = y, fill = value)) + geom_tile() + theme_bw() +
        xlab(NULL) + ylab(NULL) +
        (if(min(M) < 0){
            scale_fill_gradient2(low = "red", mid = "white", high = "blue")
        } else {
            scale_fill_gradient(low = "white", high = "red")
        }) +
        guides(fill = "none") + geom_text(aes(label = pct_value)) -> G1


    if(!is.null(true_names)){
        G1 <- G1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }


    G2 <- M |>
        as_tbl_graph() |>
        activate(edges) |>
        filter(abs(weight) < 1,  weight > threshold) |>
        ggraph() +
            geom_node_point(color = "blue", size = 2) +
            geom_node_text(aes(label = name), repel = TRUE) +
            geom_edge_link(aes(width = weight), alpha = 0.5) +
            scale_edge_width(range = c(0, 1), guide = "none") +
            theme_graph()


    G1 + G2
}


canonicalize_names <- function(s){
    case_when(s == "Rehnquist" ~ "CJ",
              s == "Stevens" ~ "AJ1",
              s == "O’Connor" ~ "AJ2",
              s == "Scalia" ~ "AJ3",
              s == "Kennedy" ~ "AJ4",
              s == "Souter" ~ "AJ5",
              s == "Thomas" ~ "AJ6",
              s == "Ginsburg" ~ "AJ7",
              s == "Breyer" ~ "AJ8",
              s == "Roberts" ~ "CJ",
              s == "Alito" ~ "AJ2",
              s == "Sotomayor" ~ "AJ5",
              s == "Kagan" ~ "AJ1",
              s == "Gorsuch" ~ "AJ3",
              s == "Kavanaugh" ~ "AJ4",
              s == "Barrett" ~ "AJ7",
              TRUE ~ s)
}

party <- function(s){
    ifelse(s %in% c("Breyer", "Ginsburg", "Sotomayor", "Kagan"), "D", "R")
}

get_network <- function(start_line_no) {
    M1 <- read_csv("scotus_networks.csv", n_max = 9, skip = start_line_no - 1) |>
        select(-1) |>
        mutate_all(function(x) as.numeric(gsub("%", "", x)) / 100) |>
        as.matrix()
    diag(M1) <- 1
    rownames(M1) <- colnames(M1) <- canonicalize_names(colnames(M1))
    ORDER <- c("CJ", paste0("AJ", 1:8))
    M1[ORDER, ORDER]
}


STARTING_LINES <- c(1, 12, 23, 36, 47, 58, 69, 80, 91, 102, 113, 124, 135, 146, 157, 168, 179, 190, 201, 212, 223, 234, 245, 256, 267)
TERMS <- lapply(STARTING_LINES, function(n){
    read_csv("scotus_networks.csv", n_max = 1, skip = n - 1) |> colnames() |> _[[1]]
})

AGREEMENTS <- lapply(STARTING_LINES, get_network)
AGREEMENT_TENSOR <- abind::abind(AGREEMENTS, along = 3)

P0 <- visualize_court_no_abs(AGREEMENT_TENSOR[,,6], threshold = 0.05)
P1 <- visualize_court_no_abs(ss_tpm(AGREEMENT_TENSOR)$V_hat, threshold = 0)
P2 <- visualize_court_no_abs(ss_tpm(AGREEMENT_TENSOR - mean(AGREEMENT_TENSOR))$V_hat, threshold = 0.0)
P3 <- visualize_court(ss_tpm(tensor_cusum(AGREEMENT_TENSOR))$V_hat, threshold = 0.05) # Absolute value here helps interpretation


G <- ((P0[[1]] + ggtitle("October Term 2001")) +
          (P1[[1]] + ggtitle("Baseline (Mean) Court Behavior")) +
          (P2[[1]] + ggtitle("First Variance Component")) +
          (P3[[1]] + ggtitle("First PC of Tensor CUSUM Statistic")) + plot_layout(nrow = 1)) /
    (P0[[2]] + P1[[2]] + P2[[2]] + P3[[2]] + plot_layout(nrow = 1)) + plot_annotation(caption = "Data from the SCOTUSBlog Stat Pack")


ggsave("figure_A5.pdf", G & theme(plot.title = element_text(size = 16),
                                  plot.caption = element_text(size = 16)), width = 18, height = 8)

browseURL("figure_A5.pdf")
