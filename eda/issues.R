source("data/lipidome_functions.R")
library(webchem)
library(lipidmapsR)
library(cowplot)
library(plotly)

# 1 import data
filename <- "raw-data/HDL_lipidome/mx733759_Agus_Lipids_Single-point quant_Submit_CS_08-14-23 (1).xlsx"

pos <- read_xlsx(filename,
          sheet = 3)
# istd <- readRDS("data/istd(all).rds")

filename <- "raw-data/HDL_lipidome/mx 733759_Agus_lipidomics_human plasma lipoproteins_07-2023 submit.xlsx"

istd.peakheight <- read_xlsx(path = filename, 
                             sheet = "Internal standards", 
                             range = "A9:AB129",
                             na = "na")

# 2 process pos
colnames(pos) <- sub("Agus.*_.*_.*_(.*)", "\\1", colnames(pos))
colnames(pos)[28:33] <- c(paste0("MB00", 1:3), paste0("Pool00", 1:3))
pos <- select(pos, 
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

dfs <- lapply(istd.peakheight$Identifier[1:10], function(id){
        id.peak <- which(istd.peakheight$Identifier == id)
        id.pos <- which(pos$identifier == id)
        df <- data.frame(
                raw = as.numeric(istd.peakheight[id.peak,8:28]),
                quant = as.numeric(pos[id.pos, 13:33]),
                label = factor(colnames(pos)[13:33]),
                id = istd.peakheight$name[id.peak]
        )
}) %>%
        Reduce(rbind, .)

dfs <- filter(dfs, !is.na(quant))

ggplotly(ggplot(dfs,aes(raw, quant, color = id)) +
        geom_point() +
        theme_bw())

id <- "0.27_319.27"
id.peak <- which(istd.peakheight$Identifier == id)
id.pos <- which(pos$identifier == id)
df <- data.frame(
        raw = as.numeric(istd.peakheight[id.peak,8:28]),
        quant = as.numeric(pos[id.pos, 13:33]),
        label = factor(colnames(pos)[13:33]),
        id = id
)

ggplot(df,aes(raw, quant)) +
        geom_point() +
        theme_bw()

# combined ids
filter(pos, `Metabolite name` == "1_CE 18:1-d7 iSTD") %>%
        select(contains("-"))

