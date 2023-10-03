source("eda/global.R")
library(ggrepel)

lipid <- readRDS("data/targeted_singletp.rds")


# CV distribution ---------------------------------------------------------

hist(lipid$fdata$qc_cv)
filter(lipid$fdata, qc_cv>30) %>% select(Annotation)


# blank read --------------------------------------------------------------

table(lipid$fdata$sample_mean < lipid$fdata$blank_mean)
filter(lipid$fdata, sample_mean < blank_mean) %>% select(Annotation)


# PCA ---------------------------------------------------------------------

edata.prop <- apply(lipid$edata, 2, function(col){
        col/sum(col)*100
})
edata.scaled <- scale(t(edata.prop), center = T)
res.pca <- PCA(edata.scaled, scale.unit = T, graph = FALSE)
sree <- fviz_eig(res.pca, addlabels = TRUE)
ind <- fviz_pca_ind(res.pca, 
                    axes = c(1,2),
                    repel = TRUE,
                    # col.ind = rownames(edata.scaled), 
                    # palette = colorRampPalette(pal_npg()(10))(15),
                    # addEllipses = TRUE, # Concentration ellipses
                    # ellipse.type = "confidence",
                    title = NULL
) +
        theme_cynthia_bw() +
        theme(
                legend.position = "none"
        ) +
        labs(title = NULL)

fviz_cos2(res.pca, choice = "var", axes = 1:2, top = 20)
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             select.var = list(contrib = 20),
             repel = TRUE # Avoid text overlapping
)


plotPie <- function(sample){
        group <- rownames(edata.prop)
        df <- data.frame(
                class = factor(group, levels = group),
                prop = edata.prop[,sample]
        ) %>% 
                mutate(ypos = 100 - cumsum(prop) + 0.5*prop)
        
        n <- nrow(df)
        if (n<=10) {
                pal <- pal_npg()(n)
        } else {
                pal <- colorRampPalette(pal_npg()(10))(n)
        }
        
        ggplot(df, aes(x="", y=prop)) +
                geom_bar(aes(fill=class), stat="identity", width=1, color = "white") +
                coord_polar("y", start=0) +
                scale_fill_manual(values = pal) +
                # geom_text(aes(y = ypos, label = lab), color = "white", size=3, nudge_x = 0.3) +
                labs(title = sample) +
                theme_cynthia() +
                geom_label_repel(aes(label = class, y = ypos), data = filter(df, prop>1),
                                 nudge_x = 0.6,
                                 size = 3, show.legend = F)
}

plotBar <- function(edata){
        group <- rownames(edata)
        sample_lv <- colnames(edata)
        df <- as.data.frame(edata) %>%
                rownames_to_column("class") %>%
                pivot_longer(-c("class"), names_to = "sample", values_to = "prop") %>%
                mutate(class = factor(class, levels = group),
                       sample = factor(sample, levels = sample_lv),
                       ypos = prop+2)
        n <- length(group)
        if (n<=10) {
                pal <- pal_npg()(n)
        } else {
                pal <- colorRampPalette(pal_npg()(10))(n)
        }

        ggplot(df, aes(x=class, y=prop)) +
                geom_bar(aes(fill = class),stat="identity", width=0.7, color = "black") +
                facet_wrap(~sample, nrow = 5) +
                scale_fill_manual(values = pal) +
                geom_text(aes(y = ypos, label = paste0(round(prop, 1), "%")), size=3, nudge_y = 10) +
                labs(x = "", y = "Proportion, %") +
                theme_cynthia_bw()
}

# Categories --------------------------------------------------------------

# edata
lipid.class <- summarize_feature(lipid, "Categories")
edata.prop <- apply(lipid.class$edata, 2, function(col){
        col/sum(col)*100
})

# pie chart
pies.class1 <- lapply(sampleNames(lipid), plotPie)
legend.class1 <- get_legend(pies.class1[[1]])
pies.class1 <- lapply(pies.class1, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class1, nrow = 3) 

# bar plot
plotBar(edata.prop)


# Main class --------------------------------------------------------------

lipid.class2 <- summarize_feature(lipid, "Main.class")
edata.prop <- apply(lipid.class2$edata, 2, function(col){
        col/sum(col)*100
})

# pie chart
pies.class2 <- lapply(sampleNames(lipid), plotPie)
legend.class2 <- get_legend(pies.class2[[1]])
pies.class2 <- lapply(pies.class2, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class2, nrow = 3) %>%
        plot_grid(legend.class2, nrow = 1, rel_widths = c(1, 0.2))

# bar plot
plotBar(edata.prop)

# Subclass --------------------------------------------------------------

lipid.class3 <- summarize_feature(lipid, "Lipid.subclass")
edata.prop <- apply(lipid.class3$edata, 2, function(col){
        col/sum(col)*100
})
# plotPie("Plasma_001")

# pie chart
pies.class3 <- lapply(sampleNames(lipid), plotPie)
legend.class3 <- get_legend(pies.class3[[1]]+theme(legend.position = "bottom"))
pies.class3 <- lapply(pies.class3, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class3, nrow = 3) %>%
        plot_grid(legend.class3, nrow = 2, rel_heights = c(1, 0.2))

# bar plot
plotBar(edata.prop)
