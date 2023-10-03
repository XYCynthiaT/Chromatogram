source("data/lipidome_functions.R")

# 1 Import data
filename <- "raw-data/HDL_lipidome/mx733759_Agus_Lipids_Single-point quant_Submit_CS_08-14-23_combine.xlsx"

htset <- import_wcmc_excel(filename, 
                           sheet = "pos quant submit", 
                           edata_range = "H9:AB675", 
                           pdata_range = "G1:AB7", 
                           fdata_range = "A8:G675", 
                           "InChI Key", "na")

# 2 exclude features without InChI Key
htset <- subset_features(htset, !is.na(htset$fdata$InChIKey))

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

# 7 Assign lipid classes
lipid.class <- read.csv("raw-data/HDL_lipidome/Lipid subclass nomenclature.csv")
htset <- assignClass(htset, lipid.class$Abbreviation)

# 7.2 Add additional feature info
lipid.class <- lipid.class[!duplicated(lipid.class[,"Abbreviation"]),]
htset$fdata <- htset$fdata %>%
        rownames_to_column() %>%
        left_join(lipid.class[,1:4], by = c("class"="Abbreviation")) %>%
        column_to_rownames()

# 8 Filter by CVs
htset <- filter_by_cv(htset, cv = "qc_cv")

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

# 13 Finalize

saveRDS(htset, "data/targeted_singletp_fin1.rds")