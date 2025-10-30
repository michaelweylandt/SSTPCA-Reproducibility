G_EMBEDDINGS <- G2_NO_CD + M2_G2_NO_CD +  TWIST_UU + COSIE_VV +JEG_HH +
    plot_layout(guides="collect", ncol = 2) & theme(legend.position="bottom")

ggsave("figure_A4.pdf", G_EMBEDDINGS, width=14, height=14)

if(grepl("darwin", version$os, ignore.case=TRUE)){
    system("open figure_A4.pdf")
} else {
    system("xdg-open figure_A4.pdf")
}
