library(HTSet)
library(tidyverse)
library(ggsci)
library(FactoMineR)
library(factoextra)
library(cowplot)

theme_cynthia <- function(){ 
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

emptyPlot <- ggplot() + geom_blank() + theme_void()

bar_color <- ggsci::pal_npg()(3)
bar_color[3] <- "#FCA510"
