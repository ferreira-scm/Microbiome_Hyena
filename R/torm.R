setwd('/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/')

df <- read.csv("Data/microbiome_tagged_fitness_2021-06-30.csv")
head(df)

df2 <- read.csv("Data/microbiome_tagged.csv")

all(df$sample_ID.x==df2$sample_ID.x)



df$sample_ID.x

    ## now some more data cleanig and prep we have some NA's at the
    ## social rank, for juveniles whose genetic mums were dead at hte
    ## time of sampling. We will give them the social rank of their
    ## surrogate mums



#df$sample_ID.x[df$sample_ID.x=="C47"]  <- "C0047"
#df$sample_ID.x[df$sample_ID.x=="C52"]  <- "C0052"
#df$sample_ID.x[df$sample_ID.x=="C129"]  <- "C0129"

df$CSocialRank <- rowMeans(df[, c("social_rank_hyena_ID",
                                      "social_rank_genetic_mum")],
                             na.rm=TRUE)

mytable <- read.csv("Data/Covariates_int_biomes_23-06-29.csv")

mytable <- mytable[match(df$sample_ID.x, mytable$sample_ID.x),]



all(df$sample_ID.x==mytable$sample_ID.x) # looks good

    df$ID_sib <- mytable$ID_sib
    df$genetic_mum <- mytable$genetic_mum
    df$surrogate_mum <- mytable$surrogate_mum


## B8423 and B8180 have 2 surrogate mums. Let's give them a mean
    ## social rank.  B8427, B8738 and B8833 have 1 surrogate mum and
    ## we will give them that social rank.
df[df$sample_ID.x%in%c("B8427", "B8738", "B8833", "X4350"),
           "CSocialRank"] <-
        df[df$sample_ID.x%in%c("B8427", "B8738", "B8833", "X4350"),
           "social_rank_surrogate_mum1"]

    df[df$sample_ID.x%in%c("B8423"),]$CSocialRank <-
        mean(c(df[df$sample_ID.x%in%c("B8423"),]$social_rank_surrogate_mum1,
               df[df$sample_ID.x%in%c("B8423"),]$social_rank_surrogate_mum2))

    df[df$sample_ID.x%in%c("B8180"),]$CSocialRank <-
        mean(c(df[df$sample_ID.x%in%c("B8180"),]$social_rank_surrogate_mum1,
               df[df$sample_ID.x%in%c("B8180"),]$social_rank_surrogate_mum2))

df$CSocialRank
    
    ## now let's see if we can complete immune and parasite data with
    ## Susana's dataset
comp.I <- read.csv("Data/Sample_database_Soares_20230525.csv")

id_NA <- df[is.na(df$IgA),]$hyena_ID

tb <- subset(comp.I, comp.I$ID%in%id_NA)
#tb <- tb[1:19,]

tb[, c("ID", "sample_date", "frozen_sample", "formalin_sample")]

df[is.na(df$IgA),c("hyena_ID", "date_sampling")]

    #B6493 instead of B6492 (one day later)
    
    df[df$Sample=="B6492",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- unique(
        tb[tb$formalin_sample=="B6493",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")])

#ES2736 instead of B9482 (8 days later)
    df[df$Sample=="B9482",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- tb[tb$frozen_sample=="ES2736",
                                 c("Ancylostoma_egg_load",
                                   "Cystoisospora_oocyst_load",
                                   "Spirometra_egg_load",
                                   "Dipylidium_egg_load",
                                   "Trichuris_egg_load",
                                   "Taeniidae_egg_load",
                                   "Spirurida_egg_load",
                                   "Mesocestoides_egg_load", "mucin",
                                   "IgA", "IgG", "lysozyme",
                                   "neopterin")]
#ES887 instead of X6643 (one day later)
    metadt[metadt$Sample=="X6643",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- tb[tb$frozen_sample=="ES887",
                                 c("Ancylostoma_egg_load",
                                   "Cystoisospora_oocyst_load",
                                   "Spirometra_egg_load",
                                   "Dipylidium_egg_load",
                                   "Trichuris_egg_load",
                                   "Taeniidae_egg_load",
                                   "Spirurida_egg_load",
                                   "Mesocestoides_egg_load",
                                   "mucin",
                                   "IgA",
                                   "IgG",
                                   "lysozyme",
                                   "neopterin")]
    



head(df)

names(df)

df <-df[,1:52]

df$module <- NULL
df$sub_module <- NULL
df$litter_size <- NULL
df$twin_status <- NULL
df$sex_comp <- NULL
df$date_reversal1 <- NULL
df$date_reversal2 <- NULL
df$type_adoption <- NULL
df$genetic_mum_alive <- NULL
df$social_rank_genetic_mum <- NULL
df$surrogate_mum1_alive <- NULL
df$social_rank_genetic_mum <- NULL
df$social_rank_surrogate_mum <- NULL
df$social_rank_hyena_ID <- NULL
df$social_rank_surrogate_mum2 <- NULL
df$date_adoption <- NULL
df$prey_level <- NULL
df$prey_level_imputed <- NULL
df$sampling_month <- NULL
df$season <- NULL
df$total_parasite_load <- NULL

names(df)

df <- df[,-c(11:30)]


write.csv(df, "/SAN/Susanas_den/gitProj/Microbiome_Hyena/Data/microbiome_tagged.csv")

head(df)
