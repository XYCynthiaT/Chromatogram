library(Metabase)
library(readxl)
library(tidyverse)
library(HTSet)

# Read data ---------------------------------------------------------------

filename <- "raw-data/HDL_lipidome/mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx"

mset = import_wcmc_excel(
        file            = filename, 
        sheet           = "Data",
        conc_range      = "H10:AB3503",
        sample_range    = "G1:AB7",
        feature_range   = "A9:G3503",
        InChIKey        = "InChI Key",
        experiment_type = "Lipidomics"
)
mset

mset = subset_features(mset, !is.na(feature_data(mset)$InChIKey))

mset = collapse_QC(mset, qc_names = paste0("Pool00", 1:3))
mset = subset_samples(mset, !grepl("^MB", sampleNames(mset)))

plot_hist_NA(mset)

mset = transform_by_feature(
        mset, function(x) ifelse(is.na(x), min(x, na.rm = TRUE)/2, x)
)

feature_data(mset)$class = assign_lipid_class(feature_data(mset)$name)
# feature_data(mset)$`ESI mode` = sub("ESI ", "", feature_data(mset)$`ESI mode`)

internal_standards <- readRDS("data/internal_standards.rds")
internal_standards_peakHeight <- readRDS("data/internal_standards_peakHeight.rds")
experiment_data(mset)$institute = "West Coast Metabolomics Center"
experiment_data(mset)$sample_volumn_ul = 50
experiment_data(mset)$internal_standards = internal_standards
experiment_data(mset)

# Exclude those class == NA
mset = subset_features(
        mset,
        !is.na(mset@feature_data$class)
)
mset = calibrate_lipidomics_wcmc(mset, cid = "InChIKey", 
                                 class = "class", ESI = "ESI mode",
                                 iSTDPeakHeight = internal_standards_peakHeight)
mset

# Species didn't fit the ESI mode
sum(rowSums(conc_table(mset))==0) # 11
mset <- subset_features(
        mset,
        rowSums(mset@conc_table) !=0
)

mset = filter_by_cv(mset, cv = "qc_cv", cid = "InChIKey")
mset

lpd.untargeted.quant <- HTSet(edata = as(conc_table(mset), "matrix"), 
                              fdata = as(feature_data(mset), "data.frame"), 
                              pdata = as(sample_table(mset), "data.frame"))
saveRDS(lpd.untargeted.quant, "data/untargeted_singletp_quant_mset.rds")

# Modified functions ------------------------------------------------------

import_wcmc_excel = function(file,
                             sheet,
                             conc_range,
                             sample_range,
                             feature_range,
                             InChIKey = NULL,
                             experiment_type = "metabolomics",
                             na = "na"){
        if(!requireNamespace("readxl"))
                stop("The package 'readxl' is required for this funciton. Please install it.")
        
        experiment_type = tolower(experiment_type)
        
        # read conc_table
        conc_table = readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = conc_range,
                col_names = FALSE,
                na        = "na"
        ) %>% as.data.frame %>% as.matrix
        # read sample_table
        sample_table = readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = sample_range,
                col_names = TRUE,
                na        = ""
        ) %>%
                as.data.frame
        sample_table = column_to_rownames(sample_table,
                                          colnames(sample_table)[1]) %>%
                t %>% as.data.frame
        # read feature_data
        feature_data = readxl::read_excel(
                path      = file,
                sheet     = sheet,
                range     = feature_range,
                col_names = TRUE,
                na        = ""
        ) %>% as.data.frame
        if(!is.null(InChIKey))
                colnames(feature_data)[colnames(feature_data) == InChIKey] = "InChIKey"
        
        # assign names
        rownames(feature_data) = paste0("Feature",
                str_pad(1:nrow(feature_data), side = "left",
                        width = ceiling(log10(nrow(feature_data))), pad = "0"))
        rownames(conc_table) = rownames(feature_data)
        colnames(conc_table) = rownames(sample_table)
        
        if(experiment_type == "lipidomics"){
                object = LipidomicsSet(
                        conc_table   = conc_table(conc_table),
                        sample_table = sample_table(sample_table),
                        feature_data = feature_data(feature_data)
                )
        }else{
                object = MetabolomicsSet(
                        conc_table   = conc_table(conc_table),
                        sample_table = sample_table(sample_table),
                        feature_data = feature_data(feature_data)
                )
        }
        
        return(object)
}

calibrate_lipidomics_wcmc = function(object, class, cid, ESI, iSTDPeakHeight){
        if(!isClass(object, Class = "LipidomicsSet"))
                stop("This function only supports LipidomicsSet data.", call. = FALSE)
        
        if(is.null(object@experiment_data))
                stop("[ Metabase ] The experiment_data slot must not be null")
        
        if(is.null(object@experiment_data$internal_standards))
                stop("[ Metabase ] The experiment_data slot must contain a list item named internal_standards")
        
        if(is.null(object@experiment_data$sample_volumn_ul))
                stop("[ Metabase ] The experiment_data slot must contain `sample_volumn_ul`")
        
        if(missing(class))
                stop("[ Metabase ] Must have 'class' that specifies the feature variable of lipid class")
        
        if(!class %in% colnames(object@feature_data))
                stop("[ Metabase ] The class variable '" %+% class %+% "' not found in the feature_data")
        
        if(missing(cid))
                stop("[ Metabase ] Must have 'cid' that specifies the compound ID for each feature such as InChIKey")
        
        if(!cid %in% colnames(object@feature_data))
                stop("[ Metabase ] The cid variable '" %+% cid %+% "' not found in the feature_data")
        
        if(missing(ESI))
                stop("[ Metabase ] Must have 'ESI' that specifies the feature variable of ESI mode")
        
        if(!ESI %in% colnames(object@feature_data))
                stop("[ Metabase ] The ESI variable '" %+% ESI %+% "' not found in the feature_data",
                     call. = FALSE)
        
        is_df = object@experiment_data$internal_standards
        
        sample_vol = object@experiment_data$sample_volumn_ul
        conc_table = lapply(1:nfeatures(object), function(i){
                int = object@conc_table[i,]
                lipid.class = object@feature_data[i, class]
                if(is.na(lipid.class))
                        stop(paste0("Unknown lipid class. Feature ID: ", featureNames(object)[i]),
                             call. = FALSE)
                ESI.mode = object@feature_data[i, ESI]
                int.is = internal_standards_peakHeight[internal_standards[,ESI] == ESI.mode &
                                                        internal_standards[,class] == lipid.class,]
                if (length(int.is) == 0) {
                        conc <- vector("numeric", 15)
                        names(conc) <- names(int)
                        return(conc)
                }
                spike.amt = is_df$spike_amt[is_df[,ESI] == ESI.mode &
                                                    is_df[,class] == lipid.class]
                return(int/int.is * spike.amt / sample_vol * 1000)
        }) %>% 
                do.call(rbind, .) 
        rownames(conc_table) = featureNames(object)
        object@conc_table = conc_table(conc_table)
        object@experiment_data$conc_table_unit = "ug/ml"
        return(object)
}


filter_by_cv = function(object, cv, cid){
        if(!inherits(object, "MetabolomicsSet"))
                stop("The function only works for MetabolomicsSet data",
                     call. = FALSE)
        if(!cv %in% colnames(object@feature_data))
                stop("The argument cv not found in feature data",
                     call. = FALSE)
        if(!cid %in% colnames(object@feature_data))
                stop("The argument cid not found in feature data",
                     call. = FALSE)
        options(warn = -1)
        feature_data(object) = as(feature_data(object), "data.frame") %>%
                rownames_to_column("feature_id") %>%
                group_by(!!sym(cid)) %>%
                mutate(keep = eval(parse(text = paste0(cv, "== min(", cv, ")")))) %>%
                ungroup %>%
                as.data.frame() %>%
                column_to_rownames("feature_id") %>%
                feature_data()
        
        object = subset_features(object, feature_data(object)$keep)
        return(object)
}


file = file.path(system.file("extdata", package = "Metabase"), "wcmc_lipidomics_standards.csv")
internal_standards_sample = read.csv(file)