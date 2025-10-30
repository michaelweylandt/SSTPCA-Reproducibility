if(!file.exists("select_K_error.rds")){
    stop("Please run simulation_select_K.R to generate data before proceeding.")
}

suppressPackageStartupMessages(library(tidyverse))

# Selection of K
readRDS("select_K_error.rds") |>
    mutate(K_err = K - K_hat,
           SNR=min_normX / normE * K) |>
    group_by(dmin, T, p, problem, sd) |>
    filter(sd == 0.2, T <= 200) |>
    summarize(E=mean(K_err),
              SNR=mean(SNR)) |>
    ggplot(aes(x=SNR, y=E, color=factor(problem))) +
    geom_point() +
    geom_line() +
    facet_grid(~T, scales="free_x", labeller = label_both) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_color_brewer(type="qual", palette=3, name="Problem") +
    xlab(expression(SNR == abs(X)/abs(E))) +
    ylab(expression(K['*'] - hat(K)['SVT']))  +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 24),
          legend.title= element_text(size=24),
          legend.text = element_text(size = 20),
          strip.text = element_text(size = 22)) +
    guides(color = guide_legend(override.aes = list(size = 6))) -> G

ggsave("figure_A9.pdf", G, width=16, height=6)


if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A9.pdf")
} else {
    system("xdg-open figure_A9.pdf")
}

