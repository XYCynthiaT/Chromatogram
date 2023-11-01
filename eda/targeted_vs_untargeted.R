setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(HTSet)
library(tidyverse)
library(plotly)

targeted <- readRDS("../data/targeted_singletp_fin3_all.rds")
untargeted <- readRDS("../data/untargeted_singletp_fin4_all(single_istd).rds")$quant
istd <- readRDS("../data/istd(abundant).rds")

targeted
untargeted

# Compare annontation
table(targeted$fdata$Annotation %in% untargeted$fdata$name)
setdiff(targeted$fdata$Annotation, untargeted$fdata$name) # only in targeted
table(untargeted$fdata$name %in% targeted$fdata$Annotation)
setdiff(untargeted$fdata$name, targeted$fdata$Annotation) # only in untargeted

# Compare InChI key
table(duplicated(targeted$fdata$InChIKey)) # no duplicated InChIKey
table(duplicated(untargeted$fdata$InChIKey)) # no duplicated InChIKey
table(targeted$fdata$InChIKey %in% untargeted$fdata$InChIKey)
setdiff(targeted$fdata$InChIKey, untargeted$fdata$InChIKey) # only in targeted
table(untargeted$fdata$InChIKey %in% targeted$fdata$InChIKey)
setdiff(untargeted$fdata$InChIKey, targeted$fdata$InChIKey) # only in untargeted

# For common InChI Key
common <- intersect(untargeted$fdata$InChIKey, targeted$fdata$InChIKey)


featureNames(targeted) <- targeted$fdata$InChIKey
targeted.part1 <- targeted
targeted.part1[common[1:3]]$edata


featureNames(untargeted) <- untargeted$fdata$InChIKey
untargeted.part1 <- untargeted
untargeted.part1[common[1:3]]$edata

df <- data.frame(
        targeted = targeted.part1[common]$edata[2,],
        untargeted = untargeted.part1[common]$edata[2,]
)
ggplotly(
        ggplot(df, aes(targeted, untargeted)) +
                geom_point() +
                # geom_smooth(method = "lm") +
                theme_bw()
) # different istd cause different results

# For FAs
untargeted.FA <- subset_features(untargeted, untargeted$fdata$class == "FA")
featureNames(untargeted.FA) <- untargeted.FA$fdata$name

targeted.FA <- subset_features(targeted, targeted$fdata$class == "FA")
featureNames(targeted.FA) <- targeted.FA$fdata$Annotation

df <- data.frame(untargeted = untargeted.FA$edata["FA 18:0 (stearic acid)",], 
            targeted = targeted.FA$edata["FA 18:0 (stearic acid)",]) %>%
        rownames_to_column("sample") %>%
        mutate(sample = factor(sample, levels = sampleNames(targeted.FA)))

ggplotly(
        ggplot(df, aes(targeted, untargeted)) +
                geom_point() +
                # geom_smooth(method = "lm") +
                theme_bw()
) # the conc was amplified in untargeted data

# IDL-targeted
154014/34562*1000/450*14.96/0.05
# IDL-untargeted
192576/6391*1000/450*44.44/0.05
