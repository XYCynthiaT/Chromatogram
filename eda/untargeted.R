source("eda/global.R")

lipid <- readRDS("data/untargeted_singletp.rds")

edata.prop <- apply(lipid$edata, 2, function(col){
        col/sum(col)*100
})

plotPie <- function(sample){
        df <- data.frame(
                class = rownames(edata.prop),
                prop = edata.prop[,sample]
        )
        n <- nrow(df)
        if (n<=10) {
                pal <- pal_npg()(n)
        } else {
                pal <- colorRampPalette(pal_npg()(10))(n)
        }
        
        ggplot(df, aes(x="", y=prop, fill=class)) +
                geom_bar(stat="identity", width=1, color = "white") +
                coord_polar("y", start=0) +
                scale_fill_manual(values = pal) +
                labs(title = sample) +
                theme_cynthia()
}

# Categories --------------------------------------------------------------


lipid.class <- summarize_feature(lipid, "Categories")
edata.prop <- apply(lipid.class$edata, 2, function(col){
        col/sum(col)*100
})

pies.class1 <- lapply(sampleNames(lipid), plotPie)
legend.class1 <- get_legend(pies.class1[[1]])
pies.class1 <- lapply(pies.class1, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class1, nrow = 3) %>%
        plot_grid(legend.class1, nrow = 1, rel_widths = c(1, 0.2))


# Main class --------------------------------------------------------------

lipid.class2 <- summarize_feature(lipid, "Main.class")
edata.prop <- apply(lipid.class2$edata, 2, function(col){
        col/sum(col)*100
})

pies.class2 <- lapply(sampleNames(lipid), plotPie)
legend.class2 <- get_legend(pies.class2[[1]])
pies.class2 <- lapply(pies.class2, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class2, nrow = 3) %>%
        plot_grid(legend.class2, nrow = 1, rel_widths = c(1, 0.2))
