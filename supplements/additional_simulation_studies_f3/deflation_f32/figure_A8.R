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
    filter(variable != "X_RERR", init != 0) |>
    mutate(scenario = paste("Scenario", scenario),
           variable = case_when(
               variable == "cosU" ~ "cos(U) == sigma[min](hat(U)^T * U['*']) / sigma[min](hat(U)) * sigma[min](U['*'])",
               variable == "cosV" ~ "cos(V) == sigma[min](hat(V)^T * V['*']) / sigma[min](hat(V)) * sigma[min](V['*'])",
               variable == "X_RERR" ~ "abs(X['*'] - hat(X))[F]^2 / abs(X['*'])[F]^2"
           ),
           init = case_when(
               init == -1 ~ "'Deflation SST-PCA'",
               init == 0 ~ "'Simultaneous SST-PCA + Stable Init'",
               init == 1 ~ "'Simultaneous SST-PCA + Random Init'",
               init == 2 ~ "'Simultaneous SST-PCA + Oracle Init'"
           ),
           T = factor(T, levels=c(50, 100, 150, 200), labels=c("T == 50", "T == 100", "T == 150", "T == 200"))) |>
    filter(T == "T == 100") |>
    ggplot(aes(x=dmin, y=value, color=factor(scenario))) +
    geom_point() +
    geom_line() +
    facet_grid(variable ~ init, scales="free_y", labeller=label_parsed) +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position="bottom") +
    xlab(expression(d[3])) +
    ylab('Estimation Error') +
    scale_color_brewer(type="qual", palette=2, name="Orthogonality of True Factors") +
    theme(strip.text = element_text(size=18),
          strip.text.y = element_text(size=12),
          axis.text = element_text(size=18),
          axis.title=element_text(size=24),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size = 6))) -> G

ggsave("figure_A8.pdf", G, width=20, height=8)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A8.pdf")
} else {
    system("xdg-open figure_A8.pdf")
}



