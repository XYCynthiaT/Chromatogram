setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("tidyverse", "readxl")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}


# Clean up data--------------------#
data <- read_xlsx(
        "../raw-data/ERCC 6 std UV vs ml.xlsx",
        sheet = 1,
        range = "A2:L913"
) 

# For example:
# data <- read_xlsx(
#         "../raw-data/RR\ Triplicate\ Same Processing\ Day\ UV\ Values.xlsx",
#         sheet = 1,
#         range = "A2:F845"
# ) 


sample <- colnames(data)[seq(1, ncol(data), 2)] %>%
        str_split(" ", simplify = T) %>%
        apply(1, function(x)str_c(x, collapse = "_"))

colnames(data) <- data[1,]
data <- data[-1,] 
# Clean up data--------------------#



# Tidy data------------------------#
data_tidy_all <- data.frame()
for(i in 1:length(sample)){
        data_split <- data[, (2*i-1):(2*i)] %>%
                mutate_all(as.numeric) %>%
                mutate(sample = sample[i]) 
        data_tidy_all <- rbind(data_tidy_all, data_split) %>%
                mutate(sample = as.factor(sample))
}
# Tidy data------------------------#


# quick plot
library(ggplot2)
ggplot(data_tidy_all, aes(x = ml, y = mAU)) +
        geom_line(aes(group = sample, color = sample)) +
        labs(x = "Volume eluted (ml)", y = "UV readings (mAu)") +
        theme_bw()

# save
saveRDS(data_tidy_all, "chromatogram2.rds")