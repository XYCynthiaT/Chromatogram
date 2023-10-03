library(tidyverse)

filenames_sec <- list.files("raw-data/SEC_chromatograms", full.names = T)
filenames_tem <- list.files("raw-data/TEM_size_data", full.names = T)

## Extract subject ID
subject_sec <- substr(basename(filenames_sec), 1, 3)
subject_tem <- str_split_fixed(basename(filenames_tem), "_", 2)[,1] %>%
        ifelse(nchar(.)==2, paste0(0, .), .)
## Subjects in common
subject <- intersect(subject_sec, subject_tem)
subject <- subject[-175] # incorrect SEC format

## Import data
names(filenames_sec) <- subject_sec
names(filenames_tem) <- subject_tem

sec <- lapply(filenames_sec[subject], function(x)read.csv(x, sep = "\t", skip = 2))
tem <- lapply(filenames_tem[subject], function(x)read.csv(x))

## Mean for TEM (mm)
mean_tem <- sapply(tem, function(x)mean(x$Diameter))

## Weighted mean for HDL fraction of SEC (Volume eluted (ml))
mean_sec <- sapply(sec, function(x){
        lower <- filter(x, X.Fractions.==" 7")$ml.1
        upper <- filter(x, X.Fractions.==" 9")$ml.1
        temp <- filter(x, ml>=lower&ml<=upper)
        weighted.mean(temp$ml, temp$mAU)
})

identical(names(mean_tem), names(mean_sec))

(corres <- cor.test(mean_tem, mean_sec))

data.frame(
        subject = names(mean_tem),
        TEM = mean_tem,
        SEC = mean_sec
) %>%
        ggplot(aes(TEM, SEC)) +
        geom_point() +
        geom_smooth(method = "lm") +
        labs(x = "TEM, Mean of Diameter (nm)", y = "SEC, Mean of Volume Eluted (ml)") +
        theme_bw() +
        annotate("text", x = 11, y = 12, label = paste0("p = ", signif(corres$p.value, 2))) +
        annotate("text", x = 11, y = 12.05, label = paste0("R = ", round(corres$estimate, 2)))


