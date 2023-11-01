setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(HTSet)
library(tidyverse)
library(ComplexHeatmap)

# Which data?
lipid <- readRDS("../data/untargeted_singletp_fin4_all.rds")$peak

lipid <- subset_features(lipid, lipid$fdata$qc_cv<30)
lipid.rel <- apply(lipid$edata, 2, function(x){
        x/sum(x)
})

# background fractions
bg.fraction <- c(2, 4:10, 12:15)
plasma.fraction <- 1
# target.fraction <- 11 #D-HDL-L_013
# other.fraction <- setdiff(bg.fraction, target.fraction)
target.sp <- "Feature0109"
# other.sp <- setdiff(featureNames(lipid), target.sp)


# Visualization
all.lipid <- colSums(lipid.rel)
target.sp.conc <- lipid.rel[target.sp, ]

df <- data.frame(
        all = all.lipid,
        target.sp = target.sp.conc,
        fraction = factor(names(target.sp.conc), levels = names(target.sp.conc))
)
ggplot(df)+
        # geom_bar(aes(fraction, all), stat = "identity") +
        geom_bar(aes(fraction, target.sp), stat = "identity", fill = "red")

target.fraction <- 9
ggplot(df, aes(x="", y=target.sp)) +
        geom_jitter(width = 0.1) +
        geom_hline(yintercept = target.sp.conc[target.fraction]) +
        scale_y_log10() +
        theme_bw() +
        annotate("text", x = 0.8, y = target.sp.conc[target.fraction]*1.1, label = sampleNames(lipid)[target.fraction])+
        labs(x= "", y = "Normalized Peak Height", title = lipid$fdata[target.sp, "name"])

# is the mean of other fractions significantly different from the targeted fraction
# boxplot(lipid.rel[target.sp, other.fraction])
enrich_one_sp <- function(target.sp, target.fraction, bg.fraction, log = FALSE){
        other.fraction <- setdiff(bg.fraction, target.fraction)
        # other.sp <- setdiff(featureNames(lipid), target.sp)
        other.fraction.rel <- lipid.rel[target.sp, other.fraction]
        target.fraction.rel <- lipid.rel[target.sp, target.fraction]
        if (log) {
                other.fraction.rel <- log(other.fraction.rel)
                target.fraction.rel <- log(target.fraction.rel)
        }
        norm <- shapiro.test(other.fraction.rel)$p.value
        tres <- t.test(other.fraction.rel, 
                       mu = target.fraction.rel, 
                       alternative = "less")
        return(c(target = target.fraction.rel,
                 mean.other = mean(other.fraction.rel), shapiro = norm, 
                 pval = tres$p.value))
}
enrich_one_sp("Feature0109", 9, bg.fraction)
enrich_one_sp("Feature0109", 9, bg.fraction, TRUE)

target.fraction <- 9 # set target fraction
logp <- lapply(1:dim(lipid.rel)[1], function(sp){
        enrich_one_sp(sp, target.fraction, bg.fraction, TRUE)
}) %>%
        do.call(rbind, .)
p <- lapply(1:dim(lipid.rel)[1], function(sp){
        enrich_one_sp(sp, target.fraction, bg.fraction)
}) %>%
        do.call(rbind, .)
plot(x = logp[,4], y = p[,4])
lines(x = logp[,4], y = rep(0.05, nrow(p)), col = "red")
lines(x = rep(0.05, nrow(logp)), y = p[,4], col = "red") # for all p<0.05, logp<0.05

enrich_all_sp <- function(target.fraction, bg.fraction, log = FALSE){
        res <- lapply(1:dim(lipid.rel)[1], function(sp){
                enrich_one_sp(sp, target.fraction, bg.fraction, log)
        }) %>%
                do.call(rbind, .) %>%
                as.data.frame() %>%
                mutate(padj = p.adjust(.$pval, method = "BH"))
        rownames(res) <- featureNames(lipid)
        res <- cbind(lipid$fdata[,c("name", "class")], res)
        return(res)
}
p.one.frac <- enrich_all_sp(9, bg.fraction)
p.one.frac.log <- enrich_all_sp(9, bg.fraction, log = TRUE)


enrich_all_fraction <- function(bg.fraction, log = FALSE){
        res <- lapply(bg.fraction, function(target.fraction){
                res <- enrich_all_sp(target.fraction, bg.fraction, log)
                return(res$padj)
        }) %>%
                do.call(cbind, .)
        colnames(res) <- sampleNames(lipid)[bg.fraction]
        rownames(res) <- featureNames(lipid)
        return(res)
}
p <- enrich_all_fraction(bg.fraction, TRUE)


# Heatmap of enrich species in each fraction
mat <- as.matrix(p)
mat[mat<0.01] <- 0.01
mat <- -log(mat)

major.classes <- c("Cer", "CAR", "CE", "DG", "FA", "LPC", "LPE", "PC", "PE", "SM", "TG")
classes <- lipid$fdata$class
classes[!classes %in% major.classes] <- "Other"
cols <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
cols <- setNames(cols, unique(classes))

Heatmap(mat, name = "-log(p)",
        right_annotation = rowAnnotation(class = classes, col = list(class = cols)),
        show_row_names = FALSE)

# Enrichment of lipid classes
lipid$fdata$class[is.na(lipid$fdata$class)] <- "UA"

lapply(bg.fraction, function(target.fraction){
        p.one.frac <- enrich_all_sp(target.fraction, bg.fraction)
        featureID <- rownames(filter(p.one.frac, padj < 0.05))
        groups <- lipid$fdata[featureID,"class"]
        print(sampleNames(lipid)[target.fraction])
        print(table(groups))
        
        pre.contigency <- data.frame(
                feature = featureNames(lipid),
                isSp = factor("N", levels = c("Y", "N")),
                inFrac = factor("N", levels = c("Y", "N"))
        ) %>%
                column_to_rownames("feature")
        
        out <- sapply(unique(groups), function(x){
                pre.contigency[featureID, "inFrac"] <- "Y"
                pre.contigency[lipid$fdata$class==x, "isSp"] <- "Y"
                contigency <- table(pre.contigency)
                fisher.test(contigency)$p.value
        })
        return(out)
}) %>%
        `names<-`(sampleNames(lipid)[bg.fraction])


