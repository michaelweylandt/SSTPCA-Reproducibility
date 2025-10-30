library(tidyverse)
library(patchwork)

RESULTS1 <- readRDS("simulation_1.rds")
RESULTS2 <- readRDS("simulation_2.rds")

# Figure 4
T_i <- 40; p_grid <- c(10, 30, 50, 70, 90, 110)
rbind(RESULTS1 |> filter(T == T_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(T == T_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label)) |>
    filter(p %in% p_grid) |>
    ggplot(aes(x = d, y = X_err_scaled, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ p, labeller = label_both) + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab("d") +
    ylab(expression(group("||", hat(X) - X["*"], "||")[F] / group("||", X["*"], "||")[F])) + ggtitle("Recovery of Rank-1 Network Series") + scale_y_log10()  +
    theme(legend.position = "bottom")-> G1a

rbind(RESULTS1 |> filter(T == T_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(T == T_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label)) |>
    filter(p %in% p_grid) |>
    ggplot(aes(x = d, y = u_acos, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ p, labeller = label_both) + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab("d") +
    #ylab(expression(acos("\u2220" ~ u["*"] ~ ","~ hat(u)))) + ggtitle("Recovery of Time Factor (u)") +
    ylab(expression(acos(u["*"] ~ ","~ hat(u)))) + ggtitle("Recovery of Time Factor (u)") +
    ylim(c(0, 90)) +
    theme(legend.position = "bottom")-> G1b

rbind(RESULTS1 |> filter(T == T_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(T == T_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label)) |>
    filter(p %in% p_grid) |>
    ggplot(aes(x = d, y = v_acos, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ p, labeller = label_both) + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab("d") +
    scale_shape_discrete(name = expression(u["*"] ~ "Drawn From")) +
    #ylab(expression(acos("\u2220" ~ v["*"] ~ ","~ hat(v))))  + ggtitle("Recovery of Base Network (V)") +
    ylab(expression(acos(v["*"] ~ ","~ hat(v))))  + ggtitle("Recovery of Base Network (V)") +
    ylim(c(0, 90)) +
    theme(legend.position = "bottom")-> G1c

G <- (G1a + guides(linetype = "none", shape = "none", color = "none") + xlab(NULL) +
          ggtitle("SS-TPCA Recovery of Network Series - Fixed Number of Observations (T = 40) + Varying Network Size (p)") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())) /
    (G1b + guides(linetype = "none", shape = "none", color = "none") + xlab(NULL) +
         ggtitle(NULL) +
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_blank())) /
    (G1c + ggtitle(NULL) +  theme(strip.text = element_blank()))


if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_4.pdf")
} else {
    system("xdg-open figure_4.pdf")
}


ggsave("figure_4.pdf", G, width = 12, height = 6)

## Figure 5
p_i <- 40; T_grid <- c(10, 30, 50, 70, 90, 110)

rbind(RESULTS1 |> filter(p == p_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(p == p_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label),
           snr = d / sqrt(p * log(T))) |>
    filter(T %in% T_grid) |>
    ggplot(aes(x = snr, y = X_err_scaled, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ T, labeller = label_both, scales = "free_x") + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab(expression(paste("Signal-to-Noise Ratio: ", d / sqrt(p * " " * log(T))))) +
    ylab(expression(group("||", hat(X) - X["*"], "||")[F] / group("||", X["*"], "||")[F])) + ggtitle("Recovery of Rank-1 Network Series") + scale_y_log10()  +
    theme(legend.position = "bottom")-> H1a

rbind(RESULTS1 |> filter(p == p_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(p == p_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label),
           snr = d / sqrt(p * log(T))) |>
    filter(T %in% T_grid) |>
    ggplot(aes(x = snr, y = u_acos, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ T, labeller = label_both, scales = "free_x") + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab(expression(paste("Signal-to-Noise Ratio: ", d / sqrt(p * " " * log(T))))) +
    #ylab(expression(acos("\u2220" ~ u["*"] ~ ","~ hat(u)))) + ggtitle("Recovery of Time Factor (u)") +
    ylab(expression(acos(u["*"] ~ ","~ hat(u)))) + ggtitle("Recovery of Time Factor (u)") +
    ylim(c(0, 90)) +
    theme(legend.position = "bottom")-> H1b

rbind(RESULTS1 |> filter(p == p_i) |> mutate(u = "Positive Orthant") |> select(-iter),
      RESULTS2 |> filter(p == p_i) |> mutate(u = "Uniform")) |>
    mutate(label = R.utils::capitalize(label),
           snr = d / sqrt(p * log(T))) |>
    filter(T %in% T_grid) |>
    ggplot(aes(x = snr, y = v_acos, color = label, shape = u)) + geom_point() +
    geom_line() + facet_grid( ~ T, labeller = label_both, scales = "free_x") + theme_bw() +
    scale_color_viridis_d(name = "Initialization Method") + xlab(expression(paste("Signal-to-Noise Ratio: ", d / sqrt(p * " " * log(T))))) +
    scale_shape_discrete(name = expression(u["*"] ~ "Drawn From")) +
    #ylab(expression(acos("\u2220" ~ v["*"] ~ ","~ hat(v))))  + ggtitle("Recovery of Base Network (V)") +
    ylab(expression(acos(v["*"] ~ ","~ hat(v))))  + ggtitle("Recovery of Base Network (V)") +
    ylim(c(0, 90)) +
    theme(legend.position = "bottom")-> H1c


H <- (H1a + guides(linetype = "none", shape = "none", color = "none") + xlab(NULL) +
          ggtitle("SS-TPCA Recovery of Network Series - Varying Number of Observations (T) + Fixed Network Size (p = 40)") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())) /
    (H1b + guides(linetype = "none", shape = "none", color = "none") + xlab(NULL) +
         ggtitle(NULL) +
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_blank())) /
    (H1c + ggtitle(NULL) +  theme(strip.text = element_blank()))

ggsave("figure_5.pdf", H, width = 12, height = 6)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_5.pdf")
} else {
    system("xdg-open figure_5.pdf")
}


