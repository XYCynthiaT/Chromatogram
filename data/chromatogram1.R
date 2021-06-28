setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("tidyverse", "readxl")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}


# Clean up data
data1 <- read_xlsx(
        "../raw-data/RR\ Triplicate\ Same Processing\ Day\ UV\ Values.xlsx",
        sheet = 1,
        range = "A2:F845"
) 
colnames(data1) <- data1[1,]
data1 <- data1[-1,] 

# Tidy data
data1_tidy1 <- data1[, 1:2] %>%
        mutate(sample = 1) 
data1_tidy2 <- data1[, 3:4] %>%
        mutate(sample = 2) 
data1_tidy3 <- data1[, 5:6] %>%
        mutate(sample = 3) 
data1_tidy_all <- rbind(data1_tidy1, data1_tidy2) %>%
        rbind(data1_tidy3) %>%
        mutate_all(as.numeric) %>%
        mutate(sample = as.factor(sample))


# quick plot
library(ggplot2)
ggplot(data1_tidy_all, aes(x = ml, y = mAU)) +
        geom_line(aes(group = sample, color = sample)) +
        labs(x = "Volume eluted (ml)", y = "UV readings (mAu)") +
        theme_bw()

saveRDS(data1_tidy_all, "chromatogram1.rds")