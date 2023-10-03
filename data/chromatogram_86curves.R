setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("tidyverse", "readxl")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

file <- c("../raw-data/ERCC 2 Col 8 Curves.xls",
          "../raw-data/ERCC 2 Col 8 missing data.xls",
          "../raw-data/ERCC1 Col 6 Curves.xls",
          "../raw-data/ERCC1 Col 6 missing data.xls")
# file <- "../raw-data/ERCC1 Col 6 Curves.xls"
data <- lapply(file, function(x){
        read_xls(
                x,
                sheet = 1,
                skip = 1
        ) 
})
names(data) <- c("eight", "eight_sup", "six", "six_sup")

data_tidy_all <- vector("list", 4)
names(data_tidy_all) <- names(data)
for (nm in names(data)) {
        # Clean up data--------------------#
        dt <- data[[nm]]
        # Remove column of injection and fractions  
        fracCols <- grep("fractions", colnames(dt), ignore.case = T)
        injCols <- grep("inject", colnames(dt), ignore.case = T)
        frac_inj <- c(fracCols, fracCols+1, injCols, injCols+1)
        dt <- select(dt, !c(frac_inj))
        # Replace space to underline
        sample <- colnames(dt)[seq(1, ncol(dt), 2)] %>%
                gsub(" ", "_", .)
        # rename columns
        colnames(dt) <- dt[1,]
        dt <- dt[-1,]
        # Clean up data end--------------------#
        
        # Tidy data--------------------#
        data_tidy_all[[nm]] <- data.frame()
        for (i in 1:length(sample)) {
                data_split <- dt[,(2*i-1):(2*i)] %>%
                        mutate_all(as.numeric) %>%
                        mutate(sample = sample[i])
                data_tidy_all[[nm]] <- rbind(data_tidy_all[[nm]], data_split) %>%
                        mutate(sample=as.factor(sample))
                data_tidy_all[[nm]] <- data_tidy_all[[nm]][complete.cases(data_tidy_all[[nm]]),]
        }
        # Tidy data end-------------------#
}

# Replace wrong samples
wrong1 <- unique(data_tidy_all$eight$sample[which(data_tidy_all$eight$mAU< -20)])
wrong2 <- unique(data_tidy_all$six$sample[which(data_tidy_all$six$mAU>100)])
data_fin <- data_tidy_all[c("eight", "six")]
data_fin$eight <- filter(data_tidy_all$eight, sample != wrong1) %>%
        mutate(sample = droplevels(sample)) %>%
        rbind(data_tidy_all$eight_sup)
data_fin$six <- filter(data_tidy_all$six, sample != wrong2) %>%
        mutate(sample = droplevels(sample)) %>%
        rbind(data_tidy_all$six_sup)

# quick plot
library(ggplot2)
ggplot(data_fin$eight, aes(x = ml, y = mAU)) +
        geom_line(aes(group = sample, color = sample)) + 
        labs(x = "Volume eluted (ml)", y = "UV readings (mAU)") +
        theme_bw()
ggplot(data_fin$six, aes(x = ml, y = mAU)) +
        geom_line(aes(group = sample, color = sample)) + 
        labs(x = "Volume eluted (ml)", y = "UV readings (mAU)") +
        theme_bw()

# save
saveRDS(data_fin, "chromatogram_86curves.rds")
