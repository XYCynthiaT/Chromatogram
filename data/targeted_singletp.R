library(HTSet)
library(readxl)
library(tidyverse)


# Read data ---------------------------------------------------------------

filename <- "raw-data/HDL_lipidome/mx733759_Agus_Lipids_Single-point quant_Submit_CS_08-14-23.xlsx"

# lipid species
fdata.pos <- read_xlsx(filename, sheet = "pos quant submit", range = "B8:F477") %>%
        mutate(ESI.mode = "pos")
fdata.neg <- read_xlsx(filename, sheet = "neg quant submit", range = "B8:F206") %>%
        mutate(ESI.mode = "neg")
fdata <- rbind(fdata.pos, fdata.neg)
fdata <- as.data.frame(fdata)

# sample inform
pdata <- read_xlsx(filename, sheet = "pos quant submit", range = "F1:AA7", col_names = FALSE) %>% t()
colnames(pdata) <- pdata[1,]
rownames(pdata) <- NULL
pdata <- as.data.frame(pdata[-1,])
rownames(pdata) <- gsub(" ", "_", gsub("-", "_", pdata$Label))

# concentration
edata.pos <- read_xlsx(filename, sheet = "pos quant submit", range = "G9:AA477", col_names = FALSE, na = "na") %>%
        as.matrix()
edata.neg <- read_xlsx(filename, sheet = "neg quant submit", range = "G9:AA206", col_names = FALSE, na = "na") %>%
        as.matrix()
edata <- rbind(edata.pos, edata.neg)
colnames(edata) <- rownames(pdata)

# keep only annotated features
keep.features <- !is.na(fdata$Annotation)
fdata <- fdata[keep.features,]
edata <- edata[keep.features,]

# check duplicated feature name
sum(duplicated(fdata$Annotation))
fdata$Annotation[duplicated(fdata$Annotation)]

## select ones with the highest peak
dup.features <- fdata$Annotation[duplicated(fdata$Annotation)]
dup.features.id <- which(fdata$Annotation %in% dup.features)

conc.sum <- rowSums(edata, na.rm = TRUE)
conc.sum.uni <- tapply(conc.sum, fdata$Annotation, max)
keep.features <- which(conc.sum %in% conc.sum.uni)
fdata <- fdata[keep.features,]
edata <- edata[keep.features,]

rownames(fdata) <- fdata$Annotation
rownames(edata) <- fdata$Annotation

# separate qc and blank samples
edata_qc <- edata[,grepl("Pool", colnames(edata))]
edata <- edata[,!grepl("Pool", colnames(edata))]
edata_blank <- edata[,grepl("MB", colnames(edata))]
edata <- edata[,!grepl("MB", colnames(edata))]


# Assign lipid class ------------------------------------------------------

lipid.class <- read.csv("raw-data/HDL_lipidome/Lipid subclass nomenclature.csv")

assignAbbreviation <- function(fdata, lipid.class){
        fdata$Abbreviation <- NA
        abbrs <- unique(lipid.class$Abbreviation)
        for (abbr in abbrs) {
                pattern <- paste0("^", abbr, " ")
                fdata$Abbreviation[grepl(pattern, fdata$Annotation)] <- abbr
        }
        return(fdata)
}

fdata.abbrev <- assignAbbreviation(fdata, lipid.class)
filter(fdata.abbrev, is.na(Abbreviation)) %>% View()

# Unique abbreviation
lipid.class2 <- lipid.class[!duplicated(lipid.class$Abbreviation),]
fdata.abbrev <- left_join(fdata.abbrev, lipid.class2, by = "Abbreviation") 

# Use QC to filter features -----------------------------------------------

rownames(fdata.abbrev) <- fdata.abbrev$Annotation
identical(rownames(fdata.abbrev), rownames(edata))
identical(rownames(fdata.abbrev), rownames(edata_qc))

fdata.abbrev$qc_mean <- apply(edata_qc, 1, mean)
fdata.abbrev$qc_sd <- apply(edata_qc, 1, sd)
fdata.abbrev$qc_cv <- fdata.abbrev$qc_sd/fdata.abbrev$qc_mean*100

# QC distribution
hist(fdata.abbrev$qc_cv)

# Use blank as LOD --------------------------------------------------------

fdata.abbrev$blank_mean <- apply(edata_blank, 1, mean, na.rm = T)
fdata.abbrev$blank_sd <- apply(edata_blank, 1, sd, na.rm = T)
fdata.abbrev$blank_min <- apply(edata_blank, 1, min, na.rm = T)

fdata.abbrev$sample_mean <- apply(edata, 1, mean, na.rm = T)
fdata.abbrev$sample_min <- apply(edata, 1, min, na.rm = T)


# NAs ---------------------------------------------------------------------

# replace NA with the min value of blank samples
fillUpNAs <- function(data){
        pos <- which(is.na(data), arr.ind = T)
        data[pos] <- fdata.abbrev[rownames(pos), "blank_min"]
        pos <- which(is.infinite(data), arr.ind = T)
        data[pos] <- fdata.abbrev[rownames(pos), "sample_min"]/2
        return(data)
}
edata <- fillUpNAs(edata)

lpd.targeted <- HTSet(
        edata = edata,
        fdata = fdata.abbrev,
        pdata = pdata[1:15,]
)

saveRDS(lpd.targeted, "data/targeted_singletp.rds")
