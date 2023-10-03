library(HTSet)
library(readxl)
library(tidyverse)

# Workflow

# 1 Read data

import_wcmc_excel <- function(file,
                             sheet,
                             edata_range,
                             pdata_range,
                             fdata_range,
                             InChIKey = NULL,
                             na = "na"){
        if(!requireNamespace("readxl"))
                stop("The package 'readxl' is required for this funciton. Please install it.")
        
        # read edata
        edata <- readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = edata_range,
                col_names = FALSE,
                na        = na
        ) %>% as.matrix
        
        # read pdata
        pdata <- readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = pdata_range
        ) %>%
                as.data.frame
        if(any(duplicated(pdata[,1])))
                pdata <- pdata[!duplicated(pdata[,1]),]
                warning("Within the pdata_range, there are duplicated names in the first column. Rows with redundent names were removed.")
        pdata <- column_to_rownames(pdata, 
                                    colnames(pdata)[1]) %>%
                t %>% as.data.frame
        # read fdata
        fdata <- readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = fdata_range
        ) %>% as.data.frame
        if(!is.null(InChIKey))
                colnames(fdata)[colnames(fdata) == InChIKey] <- "InChIKey"
        
        # assign names
        rownames(fdata) <- paste0("Feature",
                                  str_pad(1:nrow(fdata), side = "left",
                                          width = ceiling(log10(nrow(fdata))), 
                                          pad = "0"))
        rownames(edata) <- rownames(fdata)
        colnames(edata) <- rownames(pdata)
        
        htset <- HTSet(edata = edata, pdata = pdata, fdata = fdata)
        
        return(htset)
}
# 2 Remove NAs
        ## subset_features()

# 3 Calculate CV from QC samples and remove QC from peak height table
collapse_qc <- function(HTobj, qc_names){
        qc <- subset_samples(HTobj, samples = sampleNames(HTobj) %in% qc_names)
        qc_mean <- apply(qc$edata, 1, mean, na.rm = T)
        qc_sd <- apply(qc$edata, 1, sd, na.rm = T)
        HTobj$fdata$qc_mean <- qc_mean
        HTobj$fdata$qc_cv <- qc_sd/qc_mean*100
        HTobj$assay$qc <- qc
        HTobj <- subset_samples(HTobj, !(sampleNames(HTobj) %in% qc_names))
        return(HTobj)
}

# 4 Calulate minimum and mean from blank sampels
collapse_blank <- function(HTobj, blank_names){
        blank <- subset_samples(HTobj, samples = sampleNames(HTobj) %in% blank_names)
        blank_min <- apply(blank$edata, 1, min, na.rm = T)
        HTobj$fdata$blank_mean <- apply(blank$edata, 1, mean, na.rm = T)
        HTobj$fdata$blank_min <- blank_min
        HTobj$assay$blank <- blank
        HTobj <- subset_samples(HTobj, !(sampleNames(HTobj) %in% blank_names))
        return(HTobj)
}

# 5 Impute NA values
fillUpNAs <- function(htset){
        data <- htset$edata
        pos <- which(is.na(data), arr.ind = T)
        data[pos] <- htset$fdata[rownames(pos), "blank_min"]
        pos <- which(is.infinite(data), arr.ind = T)
        data[pos] <- apply(data[rownames(pos), ], 1, min, na.rm = TRUE)/2
        htset$edata <- data
        return(htset)
}

# 6 Subtract background noise
removeNoise <- function(htset){
        noise <- htset$fdata$blank_mean
        htset$edata <- apply(htset$edata, 2, function(x)x-noise)
        return(htset)
}

# 6 Assign lipid classes
assignClass <- function(htset, lipid.class){
        htset$fdata$class <- NA
        abbrs <- unique(lipid.class)
        for (abbr in abbrs) {
                pattern <- paste0("^", abbr, " ")
                htset$fdata$class[grepl(pattern, htset$fdata$name)] <- abbr
        }
        htset$fdata$class[grepl("Cholesterol", htset$fdata$name)] <- "Cholesterol"
        htset$fdata$class[grepl("Ceramide", htset$fdata$name)] <- "Cer"
        htset$fdata$class[grepl("PC-O|PC O-|PC p-|PC-p", htset$fdata$name, ignore.case = T)] <- "EtherPC"
        htset$fdata$class[grepl("LPC-O|LPC O-|LPC p-", htset$fdata$name, ignore.case = T)] <- "EtherLPC"
        htset$fdata$class[grepl("PE-O|PE O-|PE p-", htset$fdata$name, ignore.case = T)] <- "EtherPE"
        htset$fdata$class[grepl("LPE-O|LPE O-|LPE p-", htset$fdata$name, ignore.case = T)] <- "EtherLPE"
        htset$fdata$class[grepl("TG-O|TG O-", htset$fdata$name, ignore.case = T)] <- "EtherTG"
        htset$fdata$class[grepl("sphingosine", htset$fdata$name, ignore.case = T)] <- "Sphingosine"
        return(htset)
}

# 7 Import iSTD concentration and peak height
        ## readRDS()

# 8 Select isomers by CVs
filter_by_cv <- function(htset, cv){
        if(!cv %in% colnames(htset$fdata))
                stop("The argument cv not found in feature data",
                     call. = FALSE)
        options(warn = -1)
        htset$fdata <- htset$fdata %>%
                rownames_to_column("feature_id") %>%
                group_by(InChIKey) %>%
                mutate(keep = eval(parse(text = paste0(cv, "== min(", cv, ")")))) %>%
                ungroup %>%
                # mutate(prefix = sub(" Isomer [A-Z]$", "", htset$fdata$name)) %>%
                # mutate(prefix = sub(" [A-Z]$", "", .$prefix)) %>%
                # group_by(prefix) %>%
                # mutate(keep = eval(parse(text = paste0(cv, "== min(", cv, ")")))) %>%
                # ungroup %>%
                as.data.frame() %>%
                column_to_rownames("feature_id") 
        
        htset <- subset_features(htset, htset$fdata$keep)
        return(htset)
}


# 9 Calculate concentration
calibrate_lipidomics_wcmc <- function(htset, istd, class_name, ESI_name, chain_name){
        if(missing(class_name))
                stop("[ HTSet ] Must have 'class' that specifies the feature variable of lipid class")
        
        if(!class_name %in% colnames(htset$fdata))
                stop(paste0("[ HTSet ] The class variable ", class_name, "not found in the fdata"))
        
        if(missing(ESI_name))
                stop("[ HTSet ] Must have 'ESI' that specifies the feature variable of ESI mode")
        
        if(!ESI_name %in% colnames(htset$fdata))
                stop(c("[ HTSet ] The ESI variable '", ESI, "' not found in the fdata"),
                     call. = FALSE)
        if(!all(htset$fdata$class %in% istd$fdata$class)) {
                idx <- which(!(htset$fdata$class %in% istd$fdata$class))
                stop(paste0("Feature #", str_c(idx, sep = ", "), " doesn't have corresponding istds."))
        }

        sample_vol <- htset$assay$sample_volumn_ul
        edata <- lapply(1:nfeatures(htset), function(i){
                class <- htset$fdata[i, class_name]
                esi <- htset$fdata[i, ESI_name]
                chain <- htset$fdata[i, chain_name]
                match.class <- grepl(paste0("^",class,"$"), istd$fdata[,class_name])
                match.esi <- grepl(esi, istd$fdata[,ESI_name])
                istd.refined <- istd[match.class&match.esi]
                match.chain <- which.min(abs(chain - istd.refined$fdata[,chain_name]))
                spike_amt <- istd.refined$fdata$spike_amt[match.chain]
                ratio <- htset$edata[i,]/istd.refined$edata[match.chain,] # divide by rows
                conc <- ratio*spike_amt/sample_vol * 1000
        }) %>% 
                do.call(rbind, .) 
        rownames(edata) <- featureNames(htset)
        htset$edata <- edata
        htset$assay$edata_unit = "ug/mL"
        return(htset)
}
