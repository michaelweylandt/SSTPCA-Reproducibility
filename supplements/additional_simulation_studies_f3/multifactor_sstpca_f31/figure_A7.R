if(!file.exists("../multi_rank_error_simul.rds")){
    stop("Please generate multi_rank_error_simul.rds before proceeding.")
}

library(tidyverse)
library(reshape2)

readRDS("../multi_rank_error_simul.rds") |>
    group_by(T, init, sd, dmin, scenario) |>
    summarize(cosU = mean(cosU),
              cosV = mean(cosV),
              X_RERR = mean(X_RERR)) |>
    reshape2::melt(id.vars=c("T", "init", "sd", "dmin", "scenario")) |>
    filter(variable != "X_RERR") |>
    mutate(scenario = paste("Scenario", scenario),
           variable = case_when(
               variable == "cosU" ~ "cos(U) == sigma[min](hat(U)^T * U['*']) / sigma[min](hat(U)) * sigma[min](U['*'])",
               variable == "cosV" ~ "cos(V) == sigma[min](hat(V)^T * V['*']) / sigma[min](hat(V)) * sigma[min](V['*'])",
               variable == "X_RERR" ~ "abs(X['*'] - hat(X))[F]^2 / abs(X['*'])[F]^2"
           ),
           T = factor(T,
                      levels=c(50, 100, 150, 200),
                      labels=c("T == 50",
                               "T == 100",
                               "T == 150",
                               "T == 200"))) |>
    filter(init==2) |>
    ggplot(aes(x=dmin, y=value, color=factor(scenario))) +
    geom_point() +
    geom_line() +
    facet_grid(variable ~ T, scales="free_y", labeller=label_parsed) +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position="bottom") +
    xlab(expression(d[3] == d[min])) +
    ylab('Estimation Error') +
    scale_color_brewer(type="qual", palette=2, name="Orthogonality of True Factors") +
    theme(strip.text = element_text(size=14),
          strip.text.y = element_text(size=12) ,
          axis.text = element_text(size=20),
          axis.title = element_text(size=24),
          legend.text = element_text(size=22),
          legend.title = element_text(size=24)) +
    guides(color = guide_legend(override.aes = list(size = 6))) -> G

ggsave("figure_A7.pdf", G, width=20, height=8)
browseURL("figure_A7.pdf")
