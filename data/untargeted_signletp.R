library(HTSet)
library(readxl)
library(tidyverse)


# Read data ---------------------------------------------------------------

filename <- "raw-data/HDL_lipidome/mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx"

# lipid species
fdata <- read_xlsx(filename, sheet = "Data", range = "A9:G3503")
fdata <- as.data.frame(fdata)

# sample inform
pdata <- read_xlsx(filename, sheet = "Data", range = "G1:AB8", col_names = FALSE) %>% t()
colnames(pdata) <- pdata[1,]
rownames(pdata) <- NULL
pdata <- as.data.frame(pdata[-1,])
rownames(pdata) <- gsub(" ", "_", gsub("-", "_", pdata$Label))

# peak height
edata <- read_xlsx(filename, sheet = "Data", range = "H10:AB3503", col_names = FALSE, na = "na") %>%
        as.matrix()
colnames(edata) <- rownames(pdata)

# keep only annotated features
keep.features <- !is.na(fdata$name)
fdata <- fdata[keep.features,]
edata <- edata[keep.features,]

# check duplicated feature name
sum(duplicated(fdata$name))
fdata$name[duplicated(fdata$name)]

## select ones with the highest peak
dup.features <- fdata$name[duplicated(fdata$name)]
dup.features.id <- which(fdata$name %in% dup.features)

peak.sum <- rowSums(edata, na.rm = TRUE)
peak.sum.uni <- tapply(peak.sum, fdata$name, max)
keep.features <- which(peak.sum %in% peak.sum.uni)
fdata <- fdata[keep.features,]
edata <- edata[keep.features,]

rownames(fdata) <- fdata$name
rownames(edata) <- fdata$name

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
                fdata$Abbreviation[grepl(pattern, fdata$name)] <- abbr
        }
        return(fdata)
}

fdata.abbrev <- assignAbbreviation(fdata, lipid.class)
filter(fdata.abbrev, is.na(Abbreviation)) %>% View()

# add some abbreviations manually
fdata.abbrev$Abbreviation[fdata.abbrev$name%in% c("13-HODE", "13-HpOTrE", "13-KODE", "16(17)-EpDPE", 
                                           "8-HETE", "cis-13-docosenoic acid")] <- "FA"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("3-beta-hydroxy-5-cholestenoic acid")] <- "ST"
fdata.abbrev$Abbreviation[fdata.abbrev$name == "13-HODE cholesteryl ester"] <- "CE"
fdata.abbrev$Abbreviation[grepl("Cholesterol", fdata.abbrev$name)] <- "CE"
fdata.abbrev$Abbreviation[grepl("Ceramide", fdata.abbrev$name)] <- "Cer"
fdata.abbrev$Abbreviation[grepl("CoQ10", fdata.abbrev$name)] <- "CoQ"
fdata.abbrev$Abbreviation[grepl("lactosylceramide d18:1/24:1 15Z", fdata.abbrev$name)] <- "GlcCer"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("N-(eicosanoyl)sphingosine")] <- "Cer"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("docosanamide", "palmitamide")] <- "PAm"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("GlcCerd14:14E/20:02OH")] <- "GlcCer"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("N-acetyldihydrosphingosine")] <- "Cer"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("palmitoyl ethanolamide")] <- "NAE"
fdata.abbrev$Abbreviation[grepl("PC-O|PC O-|PC p-", fdata.abbrev$name, ignore.case = T)] <- "EtherPC"
fdata.abbrev$Abbreviation[grepl("LPC-O|LPC O-|LPC p-", fdata.abbrev$name, ignore.case = T)] <- "EtherLPC"
fdata.abbrev$Abbreviation[grepl("PE-O|PE O-|PE p-", fdata.abbrev$name, ignore.case = T)] <- "EtherPE"
fdata.abbrev$Abbreviation[grepl("LPE-O|LPE O-|LPE p-", fdata.abbrev$name, ignore.case = T)] <- "EtherLPE"
fdata.abbrev$Abbreviation[grepl("TG-O|TG O-", fdata.abbrev$name, ignore.case = T)] <- "EtherTG"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("PC-p 18:0_22:6")] <- "EtherPC"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("Phosphatidylethanolamine alkenyl 16:0-20:4", "Phosphatidylethanolamine alkenyl 18:0-18:2")] <- "EtherPE"
fdata.abbrev$Abbreviation[fdata.abbrev$name %in% c("prostaglandin F2.alpha. dimethylamine", "PGE2")] <- "PGE"

filter(fdata.abbrev, is.na(Abbreviation))

# Unique abbreviation
lipid.class2 <- lipid.class[!duplicated(lipid.class$Abbreviation),]
fdata.abbrev <- left_join(fdata.abbrev, lipid.class2, by = "Abbreviation") 

# Clean their names
name1 <- str_split_fixed(fdata.abbrev$name, "\\|", 2)[,1]
name2 <- str_split_fixed(name1, " or ", 2)[,1]

name2[name2=="PC-O 16:0_18:1"] <- "PC O-16:0_18:1"
name2[name2=="PE-O 18:1_22:6"] <- "PE O-18:1_22:6"
name2[name2=="PC-p 18:0_22:6"] <- "PC P-18:0_22:6"
name2[name2=="GlcCerd14:14E/20:02OH"] <- "GlcCer d14:1/20:0;(2OH)"
name2[name2=="lactosylceramide d18:1/24:1 15Z"] <- "GlcCer d18:1/24:1"
name2[name2=="C16 Lactosyl Ceramide (d18:1/16:0)"] <- "GlcCer d18:1/16:0"

extractCbond <- function(name){
        parts <- unlist(str_split(name, ":"))
        nparts <- length(parts)
        if (nparts > 1) {
                carb <- as.integer(sub(".*\\D(\\d{1,2}$)", "\\1", parts[1:(nparts-1)]))
                dbond <- as.integer(sub("(^\\d{1,2})\\D.*", "\\1", parts[2:nparts]))
        } else {
                carb <- 0
                dbond <- 0
        }
        
        totalbond <- paste0(sum(carb), ":", sum(dbond))
        return(totalbond)
}
res <- sapply(name2, extractCbond)
data.frame(name = names(res), value = res) %>% View()

fdata.abbrev <- mutate(fdata.abbrev,
       cbond = res,
       shortName = paste(fdata.abbrev$Abbreviation, fdata.abbrev$cbond),
       isomer = ifelse(grepl("A$", name), "A", ifelse(grepl("B$", name), "B", "")))

# Use blank as LOD --------------------------------------------------------

fdata.abbrev$blank_mean <- apply(edata_blank, 1, mean, na.rm = T)
fdata.abbrev$blank_sd <- apply(edata_blank, 1, sd, na.rm = T)
fdata.abbrev$blank_min <- apply(edata_blank, 1, min, na.rm = T)
fdata.abbrev$blank_mean[is.na(fdata.abbrev$blank_mean)] <- 0

fdata.abbrev$sample_mean <- apply(edata, 1, mean, na.rm = T)
fdata.abbrev$sample_min <- apply(edata, 1, min, na.rm = T)

# Use QC to filter features -----------------------------------------------

rownames(fdata.abbrev) <- fdata.abbrev$name
identical(rownames(fdata.abbrev), rownames(edata))
identical(rownames(fdata.abbrev), rownames(edata_qc))

fdata.abbrev$qc_mean <- apply(edata_qc, 1, mean, na.rm = T)
fdata.abbrev$qc_sd <- apply(edata_qc, 1, sd, na.rm = T)
fdata.abbrev$qc_cv <- fdata.abbrev$qc_sd/fdata.abbrev$qc_mean*100

# QC distribution
hist(fdata.abbrev$qc_cv)

# filter features
fdata.abbrev <- filter(fdata.abbrev, !is.na(`InChI Key`)) # missing InChI Key
duplicatedShortName <- fdata.abbrev[duplicated(fdata.abbrev$shortName),"shortName"]
dupTable <- filter(fdata.abbrev, shortName %in% duplicatedShortName) %>%
        filter(!grepl(" 0:0", shortName)) %>%
        select(name, shortName, `ESI mode`, qc_cv, sample_mean) %>% 
        group_by(shortName) %>%
        mutate(rank = rank(qc_cv))
View(dupTable)
write.csv(dupTable, file = "tbl/duplicated species.csv")

discard <- filter(fdata.abbrev, isomer != "") %>% 
        select(name, shortName, `ESI mode`, qc_cv, sample_mean) %>% 
        mutate(newName = sub("A$|B$", "", name)) %>%
        group_by(newName) %>%
        mutate(rank = rank(qc_cv)) %>%
        filter(rank == 2) %>% 
        pull(name)
fdata.abbrev <- filter(fdata.abbrev, !(name %in% discard))

# NAs ---------------------------------------------------------------------

# replace NA with the min value of blank samples
fillUpNAs <- function(data){
        pos <- which(is.na(data), arr.ind = T)
        data[pos] <- fdata.abbrev[rownames(pos), "blank_min"]
        pos <- which(is.infinite(data), arr.ind = T)
        data[pos] <- fdata.abbrev[rownames(pos), "sample_min"]/2
        return(data)
}
edata <- edata[fdata.abbrev$name,]
edata <- fillUpNAs(edata)

lpd.untargeted <- HTSet(
        edata = edata,
        fdata = fdata.abbrev,
        pdata = pdata[1:15,]
)

# Calculate concentration -------------------------------------------------

istd <- readRDS("data/istd.rds")

lpd.untargeted$fdata$Abbreviation <- sub("(Cer)-.*", "\\1", lpd.untargeted$fdata$Abbreviation)
quant.idx <- lapply(unique(istd$fdata$Abbreviation), function(x){
        grep(paste0("^",x,"$"), lpd.untargeted$fdata$Abbreviation)
}) %>% 
        unlist() 

lpd.untargeted.quant <- lpd.untargeted[quant.idx,]
lpd.untargeted.quant <- subset_features(lpd.untargeted.quant, !(lpd.untargeted.quant$fdata$Abbreviation=="FA"&lpd.untargeted.quant$fdata$`ESI mode`=="ESI pos"))

# peak height ratio
edata.conc <- lapply(1:nfeatures(lpd.untargeted.quant), function(i){
        abbrev <- lpd.untargeted.quant$fdata$Abbreviation[i]
        esi <- lpd.untargeted.quant$fdata$`ESI mode`[i]
        match.abbrev <- grepl(paste0("^",abbrev,"$"), istd$fdata$Abbreviation)
        match.esi <- grepl(esi, istd$fdata$`ESI mode`)
        istd.mass <- istd$fdata$`ng in 110 uL resus`[(match.abbrev&match.esi)]
        ratio <- lpd.untargeted.quant$edata[i,]/istd$edata[(match.abbrev&match.esi),] # divide by rows
        conc <- ratio*istd.mass*1000/450/50
        return(conc)
}) %>% Reduce(rbind, .)
rownames(edata.conc) <- featureNames(lpd.untargeted.quant)
lpd.untargeted.quant$edata <- edata.conc

saveRDS(lpd.untargeted, "data/untargeted_singletp.rds")
saveRDS(lpd.untargeted.quant, "data/untargeted_singletp_quant.rds")

lpd.untargeted.quant$fdata %>%
        group_by(Abbreviation) %>%
        summarise(n=n())

summarize_feature(lpd.untargeted.quant, "Abbreviation")$edata
