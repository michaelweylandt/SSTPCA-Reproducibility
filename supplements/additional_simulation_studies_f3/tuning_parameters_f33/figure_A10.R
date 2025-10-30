if(!file.exists("select_rank_BIC_error.rds")){
    stop("Please run simulation_select_rank_BIC.R to generate data before proceeding.")
}


suppressPackageStartupMessages(library(tidyverse))


readRDS("select_rank_BIC_error.rds") |>
    group_by(T, p, sd, dmin, strategy) |>
    summarize(E=sqrt(mean(X_RERR/T))) |>
    filter(strategy != "K_MSE_R_BIC",
           sd %in% c(0.5, 1, 2),
           T <= 200) |>
    mutate(T = factor(T,
                      levels = unique(T),
                      labels = paste("T ==", unique(T))),
           strategy = case_when(
               strategy == "DOUBLE_ORACLE" ~ "Oracle Ranks and Oracle Initialization",
               strategy == "INIT_ORACLE" ~ "BIC-Selected Rank and Oracle Initialization",
               strategy == "R_BIC" ~ "BIC-Selected Rank and Random Initialization",
               TRUE ~ "SVT Selected K and BIC-Selected Rank"
           ),
           sd = factor(sd,
                       levels = unique(sd),
                       labels = paste("sigma==", unique(sd)))) |>
    ggplot(aes(x=dmin, y=E, color=strategy)) +
    geom_point() + geom_line() +
    facet_grid(sd~T, labeller = label_parsed) +
    theme_bw() +
    scale_y_log10() +
    xlab(expression(d[min])) +
    ylab(expression('Reconstruction Error: ' * abs(X['*'] - hat(X))[F]^2 / abs(X['*'])[F]^2)) +
    theme(legend.position='bottom') +
    scale_color_brewer(type="qual", palette=6, name="Estimation Strategy") +
    theme(axis.text = element_text(size = 22),
          axis.title = element_text(size = 26),
          legend.title= element_text(size=24),
          legend.text = element_text(size = 24),
          strip.text = element_text(size = 18)) +
    guides(color=guide_legend(nrow=2,
                              override.aes = list(size = 6)))-> G


ggsave("figure_A10.pdf", G, width=18, height=9)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A10.pdf")
} else {
    system("xdg-open figure_A10.pdf")
}

