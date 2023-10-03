source("data/lipidome_functions.R")
library(webchem)
library(lipidmapsR)
library(cowplot)

# 1 import data
filename <- "raw-data/HDL_lipidome/mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx"

htset <- import_wcmc_excel(filename, 
                           sheet = "Data", 
                           edata_range = "H10:AB3503", 
                           pdata_range = "G1:AB8", 
                           fdata_range = "A9:G3503", 
                           "InChI Key", "na")

# 2 exclude features without InChI Key
htset <- subset_features(htset, !is.na(htset$fdata$InChIKey))
htset

# 3 Compare blank to pool
blankvspool <- apply(htset$edata, 1, function(x){
        x[is.na(x)] <- 0
        t.test(x[16:18], x[19:21], alternative = "less")$p.value # blank < pool
})
htset$fdata$blankvspool <- round(blankvspool, 2)
htset <- subset_features(htset, htset$fdata$blankvspool<0.01)
htset

# 3 Exclude InchI Key without a CID (Use these codes for the first run)
# cid <- get_cid(htset$fdata$InChIKey, 
#                from = "inchikey", 
#                match = "first", 
#                verbose = FALSE)
# cid <- cid[!is.na(cid$cid),]
# cid <- cid[!duplicated(cid$query),]
# saveRDS(cid, "data/untargeted_singletp_fin_cid.rds")

cid <- readRDS("data/untargeted_singletp_fin_cid.rds")
htset <- subset_features(htset, htset$fdata$InChIKey %in% cid$query)
htset

# 3 Calculate the CVs of QC samples
htset <- collapse_qc(htset, paste0("Pool00", 1:3))

# 4 Calculate the min of blank samples
htset <- collapse_blank(htset, paste0("MB00", 1:3))

# 5 View the distribution of missing values
apply(htset$edata, 1, function(row){
        sum(is.na(row))/length(row)*100
}) %>% 
        hist(main = "Histogram of missingness", 
            xlab = "% of samples with missing data", 
            ylab = "Number of lipid species")

# 6 Impute missing values
htset <- fillUpNAs(htset)

# 8 Filter by CVs
htset <- filter_by_cv(htset, cv = "qc_cv")
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

# From manually curated list
lipid.class2 <- read.csv("raw-data/HDL_lipidome/Lipid subclass nomenclature.csv")
htset <- assignClass(htset, lipid.class2$Abbreviation)

# 7.1 Manually add some classes
htset$fdata$class[htset$fdata$name == "GlcCerd14:14E/20:02OH"] <- "GlcCer"
htset$fdata$class[htset$fdata$name == "lactosylceramide d18:1/24:1 15Z"] <- "LacCer"
htset$fdata$class[htset$fdata$name %in%  c("Phosphatidylethanolamine alkenyl 16:0-20:4", "Phosphatidylethanolamine alkenyl 18:0-18:2")] <- "EtherPE"

# 7.2 Add additional feature info
lipid.class2 <- lipid.class2[!duplicated(lipid.class2[,"Abbreviation"]),]
htset$fdata <- htset$fdata %>%
        rownames_to_column() %>%
        left_join(lipid.class2[,1:4], by = c("class"="Abbreviation")) %>%
        column_to_rownames()

# 7.3 chain length and double bonds

format.name <- sub("(.*)\\|.*", "\\1", htset$fdata$name)
format.name <- sub("(.*) or .*", "\\1", format.name)

chain.bond <- lapply(format.name, function(x){
        parts <- str_split(x, ":") %>%
                unlist()
        n <- length(parts)
        chainL <- sum(as.integer(sub(".*\\D(\\d{1,2}$)", "\\1", parts[2:n-1])))
        nbonds <- sum(as.integer(sub("^(\\d{1,2})\\D.*$", "\\1", parts[2:n])))
        return(c(chainL, nbonds))
}) %>%
        do.call(rbind, .) %>%
        `colnames<-`(c("chain.length", "double.bond"))
htset$fdata <- cbind(htset$fdata, chain.bond)

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

# 10 Import istd
istd <- readRDS("data/istd(all).rds")
htset$assay$sample_volumn_ul <- 50
htset$fdata$class <- sub("(Cer)-.*", "\\1", htset$fdata$class)

# 11 Prepare data for quantification
# (Not all classes had an istd)
htset$fdata$class[grepl("Cholesterol", htset$fdata$name)] <- "CE"
htset$fdata[is.na(htset$fdata$chain.length), c("chain.length", "double.bond")] <- c(0,0)

quant.idx <- lapply(unique(istd$fdata$class), function(x){
        grep(paste0("^",x,"$"), htset$fdata$class, ignore.case = T) # PC-O etc. was excluded
}) %>% 
        unlist() 

htset.quant <- htset[quant.idx,]

# 12 Quantification
htset.quant <- calibrate_lipidomics_wcmc(htset.quant, istd, 
                                         class_name = "class", 
                                         ESI_name = "ESI mode", 
                                         chain_name = "chain.length")

htset.quant$fdata$class[grepl("Cholesterol", htset.quant$fdata$name)] <- "Cholesterol"

htset.quant.class <- summarize_feature(htset.quant, "class")
# 13 Finalize
data <- list(quant = htset.quant,
             peak = htset)

saveRDS(data, "data/untargeted_singletp_fin3_001.rds")





edata.prop3 <- apply(htset.quant.class$edata, 2, function(col){
        col/sum(col)*100
})

# pie chart
pies.class3 <- lapply(sampleNames(htset), plotPie, edata.prop3)
legend.class3 <- get_legend(pies.class3[[1]]+theme(legend.position = "bottom"))
pies.class3 <- lapply(pies.class3, function(x)x+theme(legend.position = "none"))
plot_grid(plotlist = pies.class3, nrow = 3)
