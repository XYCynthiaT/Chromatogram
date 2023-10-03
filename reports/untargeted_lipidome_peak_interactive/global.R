library(HTSet)
library(tidyverse)
library(ggsci)
library(cowplot)

# Data
lipid <- readRDS("../../data/untargeted_singletp_fin1.rds")
lipid.quant.name <- featureNames(lipid$quant)
lipid <- lipid$peak

# New data
lipid2 <- subset_features(lipid, lipid$fdata$qc_cv<30)
lipid2 <- subset_features(lipid2, !(featureNames(lipid2) %in% lipid.quant.name))

# Protein amount
lipid2 <- transform_by_feature(lipid2, function(x){
        ifelse(is.na(lipid2$pdata$Protein_amt),
               x, x/lipid2$pdata$Protein_amt)
})
lipidNames <- lipid2$fdata$name

# theme_bw italic x -------------------------------------------------------

theme_cynthia_bw <- function(){ 
        font <- "Helvetica"
        
        theme_bw(base_size = 14) %+replace%    #replace elements we want to change
                
                theme(
                        
                        ##grid elements
                        # panel.border = element_rect(colour = "black", fill = NA),
                        legend.position = "none",
                        panel.background = element_blank(), 
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.y = element_line(linetype = "dashed"),
                        axis.text.y = element_text(colour = "black", family = font),
                        axis.text.x = element_text(colour = "black", family = font, face = "italic", angle = 45, vjust = 1, hjust = 1),
                        axis.title.x = element_text(colour = "black", family = font, face = "bold"),
                        axis.title.y = element_text(colour = "black", family = font, face = "bold", angle = 90, vjust = 3), #change hjust to 0.5
                        strip.text.x = element_text(family = font, size = 10),
                        strip.background.x = element_rect(fill = NA),
                        # plot.margin = unit(c(0.2,0.5,0,0.5), "cm"),
                        plot.title = ggplot2::element_text(family=font,
                                                           size=12,
                                                           face="bold",
                                                           color="black", hjust = 0.5, vjust = 0),
                        plot.subtitle = ggplot2::element_text(family=font,
                                                              size=10,
                                                              margin=ggplot2::margin(9,0,9,0), hjust = 0.5)
                        
                )
}


# plotBar -----------------------------------------------------------------

plotBars <- function(feature){
        name <- lipid2$fdata[feature,"name"]
        df <- cbind(feature = lipid2$edata[feature,], lipid2$pdata) %>%
                mutate(Treatment = factor(Treatment, levels = unique(Treatment)))
        print(feature)
        print(lipid2$edata[1:2,])
        print(lipid2)
        y.plasma <- df["Plasma_001", "feature"]
        y.Inj.Lp <- df["Inj Lp_003", "feature"]
        y.Inj.DLp <- df["Inj DLp_011", "feature"]
        ypos <- max(c(y.plasma, y.Inj.Lp, y.Inj.DLp))
        df <- df[c(2, 4:10, 12:14),] # D-Albumin looks weired
        ggplot(df, aes(x=Treatment, y=feature)) +
                geom_bar(aes(fill = Treatment),stat="identity", width=0.7, color = "black") +
                geom_hline(yintercept = c(y.plasma, y.Inj.Lp, y.Inj.DLp)) +
                scale_fill_manual(values = colorRampPalette(pal_npg()(10))(12)) +
                annotate("text",label = c("Plasma", "Inj Lp", "Inj DLp"),
                         x = 11,
                         y = c(y.plasma, y.Inj.Lp, y.Inj.DLp)+ypos*0.05) +
                labs(x = "", title = name, y = "Normalized Peak Height") +
                theme_cynthia_bw()
}
