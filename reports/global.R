library(HTSet)
library(tidyverse)
library(ggsci)
library(cowplot)
library(factoextra)
library(FactoMineR)

# plotPie -----------------------------------------------------------------

plotPie <- function(sample, edata){
        df <- data.frame(
                class = factor(rownames(edata), 
                               levels = rownames(edata)),
                value = edata[,sample]
        ) %>%
                mutate(group = sub(".* \\[(.*)\\]", "\\1", .$class),
                       ypos = 100-cumsum(.$value)+ .$value/2,
                       group = ifelse(value>1, group, ""))
        n <- nrow(edata)
        if (n>10) {
                pal <- colorRampPalette(pal_npg()(10))(n)
        } else {
                pal <- pal_npg()(n)
        }
        ggplot(df, aes(x = "", y = value, fill = class)) +
                geom_bar(stat = "identity", color = "white") +
                geom_label_repel(aes(y = ypos, label = group), fill = "white", max.overlaps = Inf) +
                coord_polar("y", start = 0) +
                theme_cynthia_void() +
                scale_fill_manual(values = pal) +
                labs(title = sample)
}


# plotBar -----------------------------------------------------------------

plotBar <- function(edata){
        df <- as.data.frame(edata) %>%
                rownames_to_column("class") %>%
                mutate(class = factor(class, levels = rownames(edata))) %>%
                pivot_longer(1:ncol(edata)+1, 
                             names_to = "sample", values_to = "prop") %>%
                mutate(label = paste0(round(prop, 1), "%"),
                       sample = factor(sample, levels = colnames(edata)))
        n <- nrow(edata)
        if (n>10) {
                pal <- colorRampPalette(pal_npg()(10))(n)
        } else {
                pal <- pal_npg()(n)
        }
        ggplot(df, aes(class, prop, fill = class)) +
                geom_bar(stat = "identity", color = "black") +
                geom_text(aes(label = label), nudge_y = 7) +
                facet_wrap(~sample, ncol = 3) +
                theme_cynthia_bw() +
                labs(x="", y = "Proportion, %") +
                scale_fill_manual(values = pal)
}
# theme_void --------------------------------------------------------------

theme_cynthia_void <- function(){ 
        font <- "Helvetica"
        
        theme_void(base_size = 14) %+replace%    #replace elements we want to change
                
                theme(
                        
                        ##grid elements
                        # panel.border = element_rect(colour = "black", fill = NA),
                        # legend.position = "none", 
                        panel.background = element_blank(), 
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        # panel.grid.major.y = element_line(linetype = "dashed"), 
                        # axis.text.y = element_text(colour = "black", family = font, face = "bold"),
                        # axis.text.x = element_text(colour = "black", family = font, face = "italic"), 
                        # axis.title.x = element_text(colour = "black", family = font, face = "bold"),
                        # axis.title.y = element_text(colour = "black", family = font, face = "bold", angle = 90, vjust = 3, size = 12), #change hjust to 0.5
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


# theme_bw normal ---------------------------------------------------------

theme_cynthia_bw2 <- function(){ 
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
                        axis.text.x = element_text(colour = "black", family = font, angle = 45, vjust = 1, hjust = 1),
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


# empty plot --------------------------------------------------------------

emptyPlot <- ggplot() + geom_blank() + theme_void()

