source("data/lipidome_functions.R")
library(webchem)
library(lipidmapsR)
library(cowplot)

# 1 Import data
filename <- "raw-data/HDL_lipidome/mx733759_Agus_Lipids_Single-point quant_with cholesterol_Submit_10-17-23 combine.xlsx"

htset <- import_wcmc_excel(filename, 
                           sheet = "positive quant", 
                           edata_range = "H9:AB676", 
                           pdata_range = "G1:AB7", 
                           fdata_range = "A8:G676", 
                           "InChI Key", "na") #667

# 2 exclude features without InChI Key
htset <- subset_features(htset, !is.na(htset$fdata$InChIKey)) #635

# 3 Compare blank to pool
# blankvspool <- apply(htset$edata, 1, function(x){
#         x[is.na(x)] <- 0
#         t.test(x[16:18], x[19:21], alternative = "less")$p.value # blank < pool
# })
# htset$fdata$blankvspool <- round(blankvspool, 2)
# htset <- subset_features(htset, htset$fdata$blankvspool<0.05)
htset #454

# 4 Exclude InchI Key without a CID 
# cid <- get_cid(htset$fdata$InChIKey,
#                from = "inchikey",
#                match = "first",
#                verbose = FALSE)
# cid <- cid[!is.na(cid$cid),]
# cid <- cid[!duplicated(cid$query),] #483
# saveRDS(cid, "data/targeted_singletp_fin_cid.rds")

cid <- readRDS("data/targeted_singletp_fin_cid.rds")
htset <- subset_features(htset, htset$fdata$InChIKey %in% cid$query)
htset #391

# 5 Calculate the CVs of QC samples
htset <- collapse_qc(htset, paste0("Pool00", 1:3))

# 6 Calculate the min of blank samples
htset <- collapse_blank(htset, paste0("MB00", 1:3))

# 7 View the distribution of missing values
apply(htset$edata, 1, function(row){
        sum(is.na(row))/length(row)*100
}) %>% 
        hist(main = "Histogram of missingness", 
             xlab = "% of samples with missing data", 
             ylab = "Number of lipid species")

# 8 Filter by CVs
htset <- filter_by_cv(htset, cv = "qc_cv") #353
htset
View(htset$fdata)

# 7 Assign lipid classes
htset$fdata <- rownames_to_column(htset$fdata) %>%
        left_join(cid, by = c("InChIKey" = "query")) %>%
        column_to_rownames()

View(htset$fdata)

# From lipid maps
lipid.class <- lapply(htset$fdata$cid, function(x){
        classi <- compound_search(input_item = "pubchem_cid",
                                  input_value = x,
                                  output_item = "classification")
        if (length(classi)==0) {
                classi <- c("input" = x, "core" = NA, "main_class" = NA, "sub_class" = NA)
                return(classi)
        } else {
                return(classi[c(1, 5:7)])
        }
}) %>%
        do.call(rbind, .)
colnames(lipid.class)[2:4] <- paste0("LM_", colnames(lipid.class)[2:4])
if(!identical(htset$fdata$cid, lipid.class[,1]))
        stop("Check the order of CID in both tables.")
htset$fdata <- cbind(htset$fdata, lipid.class[,2:4])

# From manually curated table
lipid.class2 <- read.csv("raw-data/HDL_lipidome/Lipid subclass nomenclature.csv")
htset <- assignClass(htset, "Annotation", lipid.class2$Abbreviation)

# 7.1 Add additional feature info
lipid.class2 <- lipid.class2[!duplicated(lipid.class2[,"Abbreviation"]),]
htset$fdata <- htset$fdata %>%
        rownames_to_column() %>%
        left_join(lipid.class2[,1:4], by = c("class"="Abbreviation")) %>%
        column_to_rownames()

# 7.2 Manually modify some classes
htset$fdata$class[grepl("Cer-.*", htset$fdata$class)] <- "Cer"


# 9 Add additional sample quantities
filename <- "raw-data/HDL_lipidome/Sample-list-zivkovic_agus_ERCC2.xlsx"
sample.list <- read_xlsx(filename, range = "A10:H25", na = "N/A") %>%
        select(5,7,8) %>%
        `colnames<-`(c("Protein_conc", "Sample_volume", "ID")) %>%
        mutate(Protein_amt = Protein_conc * Sample_volume) %>%
        column_to_rownames("ID")
if(!identical(rownames(sample.list), sampleNames(htset)))
        stop("Rownames of sample.list didn't match sample names of the htset object.")
htset$pdata <- cbind(htset$pdata, sample.list)

# 10 Exclude EtherXX
htset <- subset_features(htset, !grepl("^Ether", htset$fdata$class)) #342

# 13 Finalize

saveRDS(htset, "data/targeted_singletp_fin3_all.rds")



library(ggrepel)
library(ggsci)
htset.quant.class <- summarize_feature(htset, "class")
edata.prop3 <- apply(htset.quant.class$edata, 2, function(col){
        col/sum(col)*100
})

# pie chart
pies.class3 <- lapply(sampleNames(htset), plotPie, edata.prop3)
legend.class3 <- get_legend(pies.class3[[1]]+theme(legend.position = "bottom"))
pies.class3 <- lapply(pies.class3, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class3, nrow = 3)