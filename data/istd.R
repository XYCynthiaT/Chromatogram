source("data/lipidome_functions.R")


# Read data ---------------------------------------------------------------

filename <- "raw-data/HDL_lipidome/mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx"

istd.peakheight <- import_wcmc_excel(file = filename, sheet = "Internal standards", 
                  edata_range = "H10:AB129",
                  pdata_range = "G1:AB8",
                  fdata_range = "A9:G129", 
                  InChIKey = "InChI Key", na = "na")

# Internal standard concentration
itsd.conc <- "raw-data/HDL_lipidome/Lipidomics_Single_Point_Quant_Concentrations_JW_CT_revised.xlsx"
conc <- read_xlsx(itsd.conc, sheet = "Sheet1", range = "A1:E79") %>%
        as.data.frame(conc)

# Calculate spike_amt
conc <- mutate(conc, spike_amt = `ng in 110 uL resus`*1000/450,
               spike_unit = "ug") %>%
        select(1, 6:7)

# Format istd.peakheight$fdata$name
istd.peakheight$fdata$name <- sub("1_(.*) iSTD", "\\1", istd.peakheight$fdata$name)
istd.peakheight$fdata$name <- sub(" \\(.*\\)", "", istd.peakheight$fdata$name)

str_split_fixed(istd.peakheight$fdata$name, " ", 2)

conc$iSTD <- sub("Lyso ", "L", conc$iSTD)

correctOrder <- function(x){
        split <- str_split_fixed(x, " ", 2)
        head <- grepl("^[A-Z]", split)
        tail <- gsub("-", "\\/", split[!head])
        if (grepl("-d\\d$", split[head])) {
                deuterium <- sub("[A-Z].+(-d\\d$)", "\\1", split[head])
                split[head] <- sub("-d\\d", "", split[head])
        } else {
                deuterium <- ""
        }
        
        paste0(split[head], " ", tail, deuterium)
}
conc$iSTD <- sapply(conc$iSTD, correctOrder)
conc$iSTD <- sub("/d", "-d", conc$iSTD)

# Check names
table(conc$iSTD %in% istd.peakheight$fdata$name)
conc$iSTD[!(conc$iSTD %in% istd.peakheight$fdata$name)]
table(istd.peakheight$fdata$name %in% conc$iSTD)
istd.peakheight$fdata$name[!(istd.peakheight$fdata$name %in% conc$iSTD)]

istd.peakheight$fdata$name[istd.peakheight$fdata$name=="SM 24:1 d18:1/24:1-d95"] <- "SM 24:1 d18:1/24:1-d9"
istd.peakheight <- subset_features(istd.peakheight, istd.peakheight$fdata$name != "CUDA")
table(istd.peakheight$fdata$name %in% conc$iSTD)

temp.fdata <- left_join(istd.peakheight$fdata, conc, by = c("name"="iSTD"))
if(!identical(istd.peakheight$fdata$name, temp.fdata$name))
        stop("The row order in the two tables doesn't match. Check the row order of both table.")
istd.peakheight$fdata <- cbind(istd.peakheight$fdata, temp.fdata[,c("spike_amt", "spike_unit")])

# Add classes
istd.peakheight$fdata$class <- sub("(^[A-Za-z]*) .*", "\\1", istd.peakheight$fdata$name)


# Summarise chain length and double bonds
format.name <- sub("(^[A-Za-z]* )\\d*:\\d.* (.*)", "\\1\\2", istd.peakheight@fdata$name)
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
istd.peakheight$fdata <- cbind(istd.peakheight$fdata, chain.bond)

# QC and blank
istd.peakheight <- collapse_qc(istd.peakheight, qc_names = paste0("Pool00", 1:3))
istd.peakheight <- collapse_blank(istd.peakheight, blank_names = paste0("MB00", 1:3))

# Remove  got neg values?
# istd.peakheight <- removeNoise(istd.peakheight)

# Compare the mean of blank samples and pooled samples
qcBlank <- cbind(istd.peakheight$assay$blank$edata, istd.peakheight$assay$qc$edata)
blankvspool <- apply(qcBlank, 1, function(x){
        x[is.na(x)] <- 0
        t.test(x[1:3], x[4:6], alternative = "less")$p.value # blank < pool
})
istd.peakheight$fdata$blank.vs.pool <- round(blankvspool, 1)

blankvssample <- apply(qcBlank, 1, function(x){
        t.test(x[1:3], x[4:6], alternative = "less")$p.value # blank < sample
})
istd.peakheight$fdata$blank.vs.sample <- round(blankvssample, 1)


# Select the istd with lowest p.value in each class category and ESI mode
istd.peakheight$fdata <- istd.peakheight$fdata %>% 
        rownames_to_column() %>%
        group_by(class, `ESI mode`) %>%
        mutate(min = min(blank.vs.pool),
               keep = blank.vs.pool==min) %>%
        column_to_rownames()
istd.peakheight <- subset_features(istd.peakheight, istd.peakheight$fdata$keep)

saveRDS(istd.peakheight, "data/istd(all).rds")


# Select the istd with lowest CV
istd.peakheight$fdata <- istd.peakheight$fdata %>% 
        rownames_to_column() %>%
        group_by(class, `ESI mode`) %>%
        mutate(rank = rank(qc_cv)) %>%
        column_to_rownames()
istd.peakheight <- subset_features(istd.peakheight, istd.peakheight$fdata$rank == 1)

saveRDS(istd.peakheight, "data/istd(pval).rds")
