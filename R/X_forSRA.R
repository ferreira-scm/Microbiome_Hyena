## The FORMAT for SRA submission
format <- read.table("Data/Metagenome.environmental.1.0.tsv", header=TRUE)
## only the mandatory fields
BF <- format[, grepl("X\\.", colnames(format))]
colnames(BF) <-  gsub("X\\.", "", colnames(BF))

## The metadata for the hyena INDIVIDUALS
metadata <- read.csv("Data/microbiome_tagged.csv", header=TRUE)

## Some IDs used somehow for naming the files
IDs1 <- read.csv("Data/Index_Pool_1.csv")
IDs2 <- read.csv("Data/Index_Pool_2.csv")

colnames(IDs1) <- gsub("_1", "", colnames(IDs1))
IDs1$Pool <- "ONE_1"

colnames(IDs2) <- gsub("_2", "", colnames(IDs2))
IDs2$Pool <- "Two_2"

IDs <- rbind(IDs1[, colnames(IDs2)], IDs2)


### THE FILENAMES 
path <- c(
    ## Hyena Pool 2 (Single amplicon run) 2nd Full sequencing Run (Good run)
    "2018_22_hyena", 
    ## Hyena Pool 1 (Multiamplicon run) Preliminary test sequencing Run
    "2018_22_Hyena",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 1
    "2018_22_hyena1_extra3_part1",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 2
    "2018_22_hyena1_extra3_part2",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 3
    "2018_22_hyena1_extra3_part3",
    ## Hyena Pool 2 (Single amplicon run) Preliminary test sequencing Run
    "2018_22_hyena2_main_run_1",
    ## Hyena Pool 2 (Single amplicon run) 1st Full sequencing Run (Susan report some problems and few reads)
    "2018_22_hyena2_run2",
    ## Hyena Pool 1 (Multiamplicon run) Full sequencing Run
    "2018_22_hyena_main_run",
    ## Hyena Pool 2 (Single amplicon run) Extra sequencing Run for more reads 
    "2018_22_P2_extra") 

fullpath <- paste0("/SAN/Victors_playground/Metabarcoding/AA_Hyena/", path)
names(fullpath) <- path

fastqFiles <- list.files(fullpath, pattern=".fastq.gz$", full.names=FALSE) 

### WTF is going on here... duplicated filenaems!
fastqFiles <- unique(fastqFiles)

## this is what we have to fill
colnames(BF)


submission <- data.frame(sample_name = fastqFiles,
                         organism = "metagenomic",
                         host = "Crocuta crocuta",
                         collection_date = "missing",
                         geo_loc_name = "Serengeti National Park (Tansania)",
                         lat_lon = "-2.32_34.83")

colnames(sumbission) <- paste0("*", colnames(submission))

write.table(submission,
            "Data/Metagenome.environmental.1.0.HYENA.tsv",
            row.names=FALSE)
