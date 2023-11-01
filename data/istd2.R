source("data/lipidome_functions.R")


# Read data ---------------------------------------------------------------

filename <- "raw-data/HDL_lipidome/mx733759_Agus_Lipids_Single-point quant_Submit_CS_08-14-23 (1).xlsx"

pos <- read_xlsx(filename,
                 sheet = 3)
neg <- read_xlsx(filename,
                 sheet = 6)
pos.neg <- list(
        pos = pos,
        neg = neg
)
pos.neg <- lapply(pos.neg, function(x){
        colnames(x) <- sub("Agus.*_.*_.*_(.*)", "\\1", colnames(x))
        colnames(x)[28:33] <- c(paste0("MB00", 1:3), paste0("Pool00", 1:3))
        x <- select(x, 
                      1:12, 
                      ends_with("-001"),
                      ends_with("-002"),
                      ends_with("-003"),
                      ends_with("-004"),
                      ends_with("-005"),
                      ends_with("-006"),
                      ends_with("-007"),
                      ends_with("-008"),
                      ends_with("-009"),
                      ends_with("-010"),
                      ends_with("-011"),
                      ends_with("-012"),
                      ends_with("-013"),
                      ends_with("-014"),
                      ends_with("-015"),
                      28:33)
})

identical(colnames(pos.neg$pos), colnames(pos.neg$neg))
pos.neg <- lapply(pos.neg, function(x){
        x <- filter(x, !is.na(`Metabolite name`))
})
pos.neg$pos <- mutate(pos.neg$pos, ESI = "pos")
pos.neg$neg <- mutate(pos.neg$neg, ESI = "neg")
pos.neg <- Reduce(rbind, pos.neg)

istd.peakheight <- filter(pos.neg, grepl("istd", `Metabolite name`, ignore.case = T))

# Format istd.peakheight$fdata$name
istd.peakheight$`Metabolite name` <- sub("1_(.*) iSTD", "\\1", istd.peakheight$`Metabolite name`)
istd.peakheight$`Metabolite name` <- sub(" \\(.*\\)", "", istd.peakheight$`Metabolite name`)

str_split_fixed(istd.peakheight$`Metabolite name`, " ", 2)

# Internal standard concentration
itsd.conc <- "raw-data/HDL_lipidome/Lipidomics_Single_Point_Quant_Concentrations_JW_CT_revised.xlsx"
conc <- read_xlsx(itsd.conc, sheet = "Sheet1", range = "A1:E79") %>%
        as.data.frame(conc)

# Calculate spike_amt
conc <- mutate(conc, spike_amt = `ng in 110 uL resus`*1000/450,
               spike_unit = "ug") %>%
        select(1, 6:7)

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
table(conc$iSTD %in% istd.peakheight$`Metabolite name`)
conc$iSTD[!(conc$iSTD %in% istd.peakheight$`Metabolite name`)]
table(istd.peakheight$`Metabolite name` %in% conc$iSTD)
istd.peakheight$`Metabolite name`[!(istd.peakheight$`Metabolite name` %in% conc$iSTD)]

istd.peakheight$`Metabolite name`[istd.peakheight$`Metabolite name`=="SM 24:1 d18:1/24:1-d95"] <- "SM 24:1 d18:1/24:1-d9"
istd.peakheight <- filter(istd.peakheight, istd.peakheight$`Metabolite name` != "CUDA")
table(istd.peakheight$`Metabolite name` %in% conc$iSTD)

temp.fdata <- left_join(istd.peakheight, conc, by = c("Metabolite name"="iSTD"))
if(!identical(istd.peakheight$`Metabolite name`, temp.fdata$`Metabolite name`))
        stop("The row order in the two tables doesn't match. Check the row order of both table.")
istd.peakheight <- temp.fdata

# Add classes
istd.peakheight$class <- sub("(^[A-Za-z]*) .*", "\\1", istd.peakheight$`Metabolite name`)


# Summarise chain length and double bonds
format.name <- sub("(^[A-Za-z]* )\\d*:\\d.* (.*)", "\\1\\2", istd.peakheight$`Metabolite name`)
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
istd.peakheight <- cbind(istd.peakheight, chain.bond)

# Create a HTSet
edata <- as.matrix(istd.peakheight[,grepl("\\d", colnames(istd.peakheight))])
rownames(edata) <- istd.peakheight$identifier
pdata <- data.frame(sample = colnames(edata))
rownames(pdata) <- colnames(edata)
fdata <- istd.peakheight[,!grepl("\\d", colnames(istd.peakheight))]
rownames(fdata) <- istd.peakheight$identifier
istd <- HTSet(edata=edata, fdata=fdata, pdata=pdata)

# QC and blank
istd <- collapse_qc(istd, qc_names = paste0("Pool00", 1:3))
istd <- collapse_blank(istd, blank_names = paste0("MB00", 1:3))

# Remove  got neg values?
# istd.peakheight <- removeNoise(istd.peakheight)

# Compare the mean of blank samples and pooled samples
qcBlank <- cbind(istd$assay$blank$edata, istd$assay$qc$edata)
blankvspool <- apply(qcBlank, 1, function(x){
        x[is.na(x)] <- 0
        t.test(x[1:3], x[4:6], alternative = "less")$p.value # blank < pool
})
istd$fdata$blank.vs.pool <- round(blankvspool, 1)

blankvssample <- apply(qcBlank, 1, function(x){
        t.test(x[1:3], x[4:6], alternative = "less")$p.value # blank < sample
})
istd$fdata$blank.vs.sample <- round(blankvssample, 1)


# Select the istd with lowest p.value in each class category and ESI mode
istd$fdata <- istd$fdata %>% 
        rownames_to_column() %>%
        group_by(class, ESI) %>%
        mutate(min = max(min(blank.vs.pool), 0.3),
               keep = blank.vs.pool <= min) %>% # min or < 0.3
        column_to_rownames()
istd2 <- subset_features(istd, istd$fdata$keep)
istd2$fdata <- istd2$fdata %>%
        rownames_to_column() %>%
        group_by(`Metabolite name`, ESI) %>%
        mutate(max = max(`sample ave peak height`),
               keep2 = `sample ave peak height`==max) %>%
        column_to_rownames()
istd2 <- subset_features(istd2, istd2$fdata$keep2)

saveRDS(istd2, "data/istd2(all).rds")

inter <- intersect(istd.old$fdata$name, istd2$fdata$`Metabolite name`)
istd.old$fdata$name[!(istd.old$fdata$name %in% inter)]
istd2$fdata$`Metabolite name`[!(istd2$fdata$`Metabolite name` %in% inter)]
