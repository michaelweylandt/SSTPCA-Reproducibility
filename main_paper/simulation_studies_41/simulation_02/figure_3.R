library(tidyverse)

custom_labeller <- function(labels, multi_line = FALSE){
    if(multi_line) stop("Not yet implemented: multi_line = TRUE")

    tryCatch(label_parsed(labels), error = function(e) label_value(labels))
}

readRDS("simulation_computational_convergence_results.rds") |>

    unnest(c(U_U_STAR,   U_U_HAT,   V_V_STAR,   V_V_HAT,
             U_U_STAR_F, U_U_HAT_F, V_V_STAR_F, V_V_HAT_F,
             iteration)) |>
    group_by(rank, iteration) |>
    filter(iteration <= 30) |>
    summarize(U_U_STAR   = 1 - mean(U_U_STAR),
              U_U_HAT    = 1 - mean(U_U_HAT),
              V_V_STAR   = 1 - mean(V_V_STAR),
              V_V_HAT    = 1 - mean(V_V_HAT),
              U_U_STAR_F = mean(U_U_STAR_F),
              U_U_HAT_F  = mean(U_U_HAT_F),
              V_V_STAR_F = mean(V_V_STAR_F),
              V_V_HAT_F  = mean(V_V_HAT_F),
              .groups = "drop") |>
    pivot_longer(cols=c(U_U_STAR, U_U_HAT, V_V_STAR, V_V_HAT, U_U_STAR_F, U_U_HAT_F, V_V_STAR_F, V_V_HAT_F),
                 names_to = "measure",
                 values_to = "error") |>
    mutate(component = ifelse(grepl("U", measure), "u^{(k)}", "V^{(k)}")) |>
    mutate(target    = ifelse(grepl("HAT", measure), "Computational - Distance to Final Iterate", "Statistical - Distance to True Parameter")) |>
    filter(grepl("U.*F", measure) | grepl("^V[^F]*$", measure)) |>
    ggplot(aes(x = iteration, color = factor(rank), y = error)) +
        facet_grid(target ~ component, scales = "free_y",
                   labeller = custom_labeller) +
        scale_color_viridis_d(name = "r = Rank(V)") +
        geom_line() +
        geom_point() +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_x_log10() +
        xlab("Iteration No. (Log Scale)") +
        ylab("Distance Between SS-TPCA Iterate and Target") -> G

ggsave("figure_3.pdf", G, height = 6, width = 6)


if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_3.pdf")
} else {
    system("xdg-open figure_3.pdf")
}
