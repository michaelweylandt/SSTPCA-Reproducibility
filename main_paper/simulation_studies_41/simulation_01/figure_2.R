library(tidyverse)

custom_labeller <- function(labels, multi_line = FALSE){
    if(multi_line) stop("Not yet implemented: multi_line = TRUE")

    tryCatch(label_parsed(labels), error = function(e) label_value(labels))
}

G1 <- readRDS("simulation_methods_comparison_results.rds") |>
    filter(T == 20) |>
    filter(rank == 5) |>
    group_by(method, T, p, rank, model) |>
    summarize(VV_L2 = mean(VV_L2), .groups = "drop") |>
    mutate(model = factor(model, levels = c("SBM", "Sphere-RDPG", "Positive-Sphere-RDPG", "Ball-RDPG", "Positive-Ball-RDPG", "Dirichlet-RDPG"))) |>
    ggplot(aes(y = VV_L2, x = p, color = method)) +
        geom_point() +
        facet_grid(~model) +
        scale_y_log10() +
        geom_line() +
        xlab("Number of Vertices") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_color_brewer(name = "Estimator", type = "qual", palette = 2) +
        ylab(expression(paste("Subspace Estimation Error: ", "||"*hat(V)*hat(V)^T - V["*"]*V["*"]^T*"||"[F])))
# \u2016 is better than "||" here but doesn't work well with my PDF engine

ggsave("figure_2.pdf", plot = G1, width=11.25, height = 4)


if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_2.pdf")
} else {
    system("xdg-open figure_2.pdf")
}

