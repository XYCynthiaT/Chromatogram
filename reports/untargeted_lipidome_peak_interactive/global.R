library(HTSet)
library(tidyverse)
library(ggsci)
library(cowplot)

# Data
# lipid <- readRDS("../../data/untargeted_singletp_fin1.rds")
# lipid <- readRDS("../../data/untargeted_singletp_fin3_001.rds")
lipid <- readRDS("../../data/untargeted_singletp_fin4_all.rds")
lipid.quant.name <- featureNames(lipid$quant)
lipid <- lipid$peak

# Normalize to total lipid peak height
uq.rel <- apply(lipid$edata, 2, function(x)x/sum(x)*100)
if(identical(rownames(uq.rel), featureNames(lipid))){
        lipid$edata <- uq.rel
}

# New data
lipid2 <- subset_features(lipid, lipid$fdata$qc_cv<30)
lipid2 <- subset_features(lipid2, !(featureNames(lipid2) %in% lipid.quant.name))

# Protein amount
# lipid2 <- transform_by_feature(lipid2, function(x){
#         ifelse(is.na(lipid2$pdata$Protein_amt),
#                x, x/lipid2$pdata$Protein_amt)
# })
lipidNames <- lipid2$fdata$name
# saveRDS(lipidNames, "lipidNames.rds")

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
        # y.plasma <- df["Plasma_001", "feature"]
        # y.Inj.Lp <- df["Inj Lp_003", "feature"]
        # y.Inj.DLp <- df["Inj DLp_011", "feature"]
        # ypos <- max(c(y.plasma, y.Inj.Lp, y.Inj.DLp))
        # df <- df[c(2, 4:10, 12:14),] # D-Albumin looks weired
        ggplot(df, aes(x=Treatment, y=feature)) +
                geom_bar(aes(fill = Treatment),stat="identity", width=0.7, color = "black") +
                # geom_hline(yintercept = c(y.plasma, y.Inj.Lp, y.Inj.DLp)) +
                scale_fill_manual(values = colorRampPalette(pal_npg()(10))(nrow(df))) +
                # annotate("text",label = c("Plasma", "Inj Lp", "Inj DLp"),
                #          x = 11,
                #          y = c(y.plasma, y.Inj.Lp, y.Inj.DLp)+ypos*0.05) +
                labs(x = "", title = name, y = "Relative to total peak height, %") +
                theme_cynthia_bw()
}


# plot pie ----------------------------------------------------------------

plotPie <- function(){
        lipid$fdata$major <- ifelse(featureNames(lipid) %in% lipid.quant.name, "Major", "Minor")
        lipid.major <- summarize_feature(lipid, "major")$edata
        df <- as.data.frame(lipid.major) %>%
                rownames_to_column("Type") %>%
                pivot_longer(cols = contains("_0"), names_to = "fraction", values_to = "Perc") %>%
                mutate(fraction = factor(fraction, levels = colnames(lipid.major))) %>%
                group_by(fraction) %>%
                mutate(ypos = 100 - cumsum(Perc) + 0.5 * Perc)
        ggplot(df, aes(x = "", y = Perc, fill = Type)) +
                geom_bar(stat = "identity", color = "white") +
                geom_label(aes(y = ypos, label = round(Perc, 1)), fill = "white") +
                coord_polar("y", start = 0) +
                theme_void() +
                scale_fill_npg() +
                facet_wrap(~fraction, nrow = 3)
}
