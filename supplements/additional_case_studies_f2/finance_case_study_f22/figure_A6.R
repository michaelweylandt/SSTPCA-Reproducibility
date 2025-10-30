library(quantmod)
library(Rcpp)
library(tidyverse)
library(patchwork)
library(ggraph)
library(tidygraph)

Rcpp::sourceCpp("../../../tensor_factorizations.cpp")

TICKERS <- c("EWJ", "EWT", "EWY", "FXI", "EWZ", "EWC", "EWU", "EWG", "EWL",
             "EWA", "EWW", "EWH", "INDY", "EWD", "EWQ", "EWP", "ERUS", "EWS",
             "EWI", "EIDO", "ECH", "THD", "EWN", "EPOL", "EZA", "EWM", "TUR",
             "EIS", "ENZL", "EPU", "EPHE", "EWO", "EIRL", "EWK", "IVV")

LOOKUPS <- c(Japan = "EWJ", Taiwan = "EWT", `South Korea` = "EWY", `China` = "FXI",
             Brazil = "EWZ", Canada = "EWC", `United Kingdom` = "EWU", Germany = "EWG",
             Switzerland = "EWL", Australia = "EWA", Mexico = "EWW", `Hong Kong` = "EWH",
             `India 50` = "INDY", Sweden = "EWD", France = "EWQ", Spain = "EWP",
             Russia = "ERUS", Singapore = "EWS", Italy = "EWI", Indonesia = "EIDO",
             Chile = "ECH", Thailand = "THD", Netherlands = "EWN", Poland = "EPOL",
             `South Africa` = "EZA", Malaysia = "EWM", Turkey = "TUR", Israel = "EIS",
             `New Zealand` = "ENZL", Peru = "EPU", Philippines = "EPHE", Austria = "EWO",
             `Brazil Small-Cap` = "EWZS", Ireland = "EIRL", Belgium = "EWK", `United States` = "IVV")

REVERSE_LOOKUPS <- c(EWJ = "Japan", EWT = "Taiwan", EWY = "South Korea", FXI = "China",
                     EWZ = "Brazil", EWC = "Canada", EWU = "United Kingdom", EWG = "Germany",
                     EWL = "Switzerland", EWA = "Australia", EWW = "Mexico", EWH = "Hong Kong",
                     INDY = "India 50", EWD = "Sweden", EWQ = "France", EWP = "Spain",
                     ERUS = "Russia", EWS = "Singapore", EWI = "Italy", EIDO = "Indonesia",
                     ECH = "Chile", THD = "Thailand", EWN = "Netherlands", EPOL = "Poland",
                     EZA = "South Africa", EWM = "Malaysia", TUR = "Turkey", EIS = "Israel",
                     ENZL = "New Zealand", EPU = "Peru", EPHE = "Philippines", EWO = "Austria",
                     EWZS = "Brazil Small-Cap", EIRL = "Ireland", EWK = "Belgium", IVV = "United States")

CONTINENTS <- c(EWJ = "Asia", EWT = "Asia", EWY = "Asia", FXI = "Asia",
                EWZ = "South America", EWC = "North America", EWU = "Europe", EWG = "Europe",
                EWL = "Europe", EWA = "Oceania", EWW = "North America", EWH = "Asia",
                INDY = "Asia", EWD = "Europe", EWQ = "Europe", EWP = "Europe",
                ERUS = "Europe", EWS = "Asia", EWI = "Europe", EIDO = "Asia",
                ECH = "South America", THD = "Asia", EWN = "Europe", EPOL = "Europe",
                EZA = "Africa", EWM = "Oceania", TUR = "Europe", EIS = "Israel",
                ENZL = "Oceania", EPU = "South America", EPHE = "Asia", EWO = "Europe",
                EWZS = "South America", EIRL = "Europe", EWK = "Europe", IVV = "North America")

STOCK_ENV <- new.env()

getSymbols(TICKERS, env = STOCK_ENV, from = "2000-01-01")

DATE_RANGE <- "2005::2024"

RETURNS <- do.call(merge, eapply(STOCK_ENV, function(x) na.omit(ROC(Ad(x)))))
RETURNS <- na.omit(RETURNS[DATE_RANGE, !is.na(head(RETURNS[DATE_RANGE], 1))])
colnames(RETURNS) <- gsub(".Adjusted", "", colnames(RETURNS))
colnames(RETURNS) <- REVERSE_LOOKUPS[colnames(RETURNS)]
months <- unique(substr(index(RETURNS), 1, 7))

TENSOR <- array(NA,
                dim = c(NCOL(RETURNS), NCOL(RETURNS), length(months)),
                dimnames = list(colnames(RETURNS), colnames(RETURNS), months))


for(mn_ix in seq_along(months)){
    mn <- months[mn_ix]
    TENSOR[,,mn_ix] <- cor(RETURNS[mn,])
}

LOOKUP_DF <- inner_join(data.frame(Ticker = names(CONTINENTS),
                                   Continent = CONTINENTS),
                        data.frame(Ticker = names(REVERSE_LOOKUPS),
                                   Name = REVERSE_LOOKUPS))

viz_network <- function(M, threshold = 0.0){
    M|> as.data.frame() |>
        rownames_to_column("other") |>
        pivot_longer(-other,
                     names_to="self",
                     values_to="corr") |>
        mutate(self=factor(self),
               other=factor(other)) |>
        filter(self != other) |>
        tbl_graph(edges=_) |>
        activate(nodes) |>
        inner_join(LOOKUP_DF, join_by(name == Name)) |>
        arrange(desc(Continent)) |>
        activate(edges) |>
        filter(abs(corr) < 1, abs(corr) > threshold) |>
        ggraph(layout = "linear", circular = TRUE) +
            geom_node_point(aes(color = Continent), size = 2) +
            geom_node_text(aes(label = name), repel = TRUE) +
            geom_edge_link(aes(width = corr), alpha = 0.15) +
            scale_edge_width(range = c(0, 1), guide = "none") +
            theme_graph()
}


G1 <- viz_network(ss_tpm(TENSOR, rank = 4)$V_hat, 0.03) +
        guides(color = "none")
G2 <- viz_network(ss_tpm(TENSOR - ss_tpm(TENSOR, rank = 1)$X_hat, rank = 4)$V_hat, 0.07) +
        guides(color = "none")
G3 <- viz_network(ss_tpm(tensor_cusum(TENSOR), rank = 4)$V_hat, 0.10) +
        guides(color = "none")

T1 <- ggplot(data.frame(u = ss_tpm(TENSOR, rank = 4)$u_hat,
                        date = as.Date(as.yearmon(dimnames(TENSOR)[[3]]))),
             aes(x = date, y = u)) + geom_point() + geom_line() + theme_bw() +
    xlab("Date") + ylab(expression(hat(u))) + scale_x_date() + ggtitle("Market Baseline Covariance")


T2 <- ggplot(data.frame(u = ss_tpm(TENSOR - ss_tpm(TENSOR, rank = 1)$X_hat, rank = 4)$u_hat,
                        date = as.Date(as.yearmon(dimnames(TENSOR)[[3]]))),
             aes(x = date, y = u)) + geom_point() + geom_line() + theme_bw() +
    xlab("Date") + ylab(expression(hat(u))) + scale_x_date() + ggtitle("Regional Modes of Variation")


T3 <- ggplot(data.frame(u = ss_tpm(tensor_cusum(TENSOR), rank = 4)$u_hat,
                        date = as.Date(as.yearmon(dimnames(TENSOR)[[3]]))[-1]),
             aes(x = date, y = u)) + geom_point() + geom_line() + theme_bw() +
    xlab("Date") + ylab(expression(hat(u))) + scale_x_date() + ggtitle("CUSUM Analysis: European Debt Crisis")


G <- (T1 + (T2 + ylab(NULL)) + (T3 + ylab(NULL))) / (G1 + G2 + G3) + plot_layout(heights = c(2,5))


ggsave("figure_A6.pdf", G & theme(text = element_text(size = 18)), width = 20, height = 6)
browseURL("figure_A6.pdf")

## Cusum analysis:
#- Look for European debt crisis, should be low-loaded on Asia
#- 2008 <-- CUSUM
