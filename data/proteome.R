
pkgs <- c("tidyverse", "readxl")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

file <- list.files("raw-data/Proteome_Byonic/", full.names = T)

f1 <- read_xlsx(file[1], sheet = 3) %>%
        `colnames<-`(sub("\r\n", " ", colnames(.)))

ggplot(f1, aes(`# of spectra`, `Total Intensity`)) + 
        geom_point()

ggplot(filter(f1, `# of spectra`<200), aes(`# of spectra`, `Total Intensity`)) + 
        geom_point()
   