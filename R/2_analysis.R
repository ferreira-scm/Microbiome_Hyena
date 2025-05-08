## this script has all the statistical analysis and data visualisation

library(ggplot2)
library(vegan)
#library(tidyverse)
library(cowplot)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(vegan)
library(phyloseq, lib.loc="/usr/local/lib/R/site-library")
library(igraph)
library(Hmisc, lib.loc="/usr/local/lib/R/site-library")
library(Matrix)
library(SpiecEasi, lib.loc="/usr/local/lib/R/site-library")
library(microbiome)
library(ggpmisc)
#library(cmdstanr)
library(brms, lib.loc="/usr/local/lib/R/site-library/")
library(dplyr)

PMS <- readRDS("tmp/fPMS.rds")

### some basic data preparation ot imput a few missing values for IgA and mucin
do_DataPrep <- FALSE
if (do_DataPrep){
    # imputing IgA and mucin
    library(mice)
    pMiss <- function(x){sum(is.na(x))/length(x)*100}
    apply(PMS@sam_data,2,pMiss)
    # to input: mucin and IgA
    names(data.frame(PMS@sam_data))
    df <- data.frame(PMS@sam_data)[29:39]
    set.seed(135)
    dfI <- mice(df, m=10, maxit=50, meth="pmm", seed=500)
    dfI <- complete(dfI, 1)
    PMS@sam_data$mucin_imputed<- dfI$mucin
    PMS@sam_data$IgA_imputed <- dfI$IgA
    PMS@sam_data$Clan <- substr(PMS@sam_data$hyena_ID, 1,1)
    saveRDS(PMS, "tmp/PMS_imputed.rds")
} else
    PMS <- readRDS("tmp/PMS_imputed.rds")


### Testing for batch effects
PMS@sam_data$batch <- (paste(PMS@sam_data$run, PMS@sam_data$Chip_Pool_1, sep="_"))
unique(paste(PMS@sam_data$run, PMS@sam_data$Chip_Pool_1, sep="_"))

adonis2(distance(PMS, method="bray") ~
            PMS@sam_data$batch,
            strata=PMS@sam_data$hyena_ID) ## We have batch effects, so we will need to account for that in the models

# We consider parasites of the gastrointestinal tract of hyenas:
Parasite <- subset_taxa(PMS, Genus%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))

##  only bacteria
Bacteria <- subset_taxa(PMS, Kingdom %in%"Bacteria") 
Bacteria

##  only Eukaryotes
Eukarya <- subset_taxa(PMS, Kingdom %in%"Eukarya") 
Eukarya

# only fungi
Fungi <- subset_taxa(PMS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))
Fungi

###################################### plotting Figure 1
## small adjustment here
Eukarya@tax_table[which(is.na(Eukarya@tax_table[,2])),2] <- "Unknown_phylum_in_Eukarya"
Eukarya@tax_table[which(Eukarya@tax_table[,2]=="Unknown_phylum"), 2] <- "Unknown_phylum_in_Eukarya"
Bacteria@tax_table[which(Bacteria@tax_table[,2]=="Unknown_phylum"), 2] <- "Unknown_phylum_in_Bacteria"
Eukarya.p <- tax_glom(Eukarya, "Phylum") # have to do this separate here otherwise the figure is too heavy

# counting the number of ASVs per phylum within eukaryotes
Euk.df <- data.frame(tax_table(Eukarya))
Euk.df %>% dplyr::count(Phylum) -> Euk.df

Eukarya.t <- microbiome::transform(Eukarya.p, "compositional")
euk.df <- psmelt(Eukarya.t)

# ordering the x axis for visualization purposes
Euk.t2 <- tax_glom(Eukarya, taxrank="Phylum")
di.m <- vegan::vegdist(Euk.t2@otu_table,
                       method="bray")
di.m[is.na(di.m)]<- 0 # defining those as 0 distances
clustering <- hclust(di.m, method="complete")
euk.df$Sample2 <- factor(euk.df$Sample, levels=clustering$labels[clustering$order])
euk.df$Sample <- as.factor(euk.df$Sample)
euk.df$Phylum <- as.factor(euk.df$Phylum)

# now defining groups as in the models
euk.df$Group<-"Other"
euk.df$Group[which(euk.df$Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))] <- "Fungi"
euk.df$Group[which(euk.df$Phylum%in%unique(Parasite@tax_table[,2]))] <- "Parasite"

levels(euk.df$Phylum) <- c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota", "Apicomplexa", "Nematozoa", "Platyhelminthes", as.vector(unique(euk.df$Phylum[euk.df$Group=="Other"])))

Euk.df$Phylum <- as.factor(Euk.df$Phylum)

levels(Euk.df$Phylum) <- levels(euk.df$Phylum)

coul <- c(c("#99760f", "#cc9d14", "#ffcb32", "#ffbf00", "#e5ab00", "#ffdc75",
            "#202c6a", "#293885", "#535f9d"), c("#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c", "#8a8a8c"))

E1 <- ggplot(Euk.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
    coord_polar()+
        scale_fill_manual(values=coul)+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")


E <- ggplot(euk.df, aes(x=Sample2, y=Abundance, fill=Phylum))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=coul)+
#    scale_fill_viridis(discrete=T)+
    labs(x="Samples", y="Relative abundance")+
    theme_classic()+
    guides(fill="none")+
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

Bac.p <- tax_glom(Bacteria, "Phylum") # have to do this separate here otherwise the figure is too heavy
Bac.t <- microbiome::transform(Bac.p, "compositional")
Bac.t2 <- tax_glom(Bacteria, taxrank="Phylum")

Bac.df <- psmelt(Bac.t)
di.m <- vegan::vegdist(Bac.t2@otu_table,
                       method="bray")
di.m[is.na(di.m)]<- 0 # defining those as 0 distances

clustering <- hclust(di.m, method="complete")
Bac.df$Sample2 <- factor(Bac.df$Sample, levels=clustering$labels[clustering$order])
Bac.df$Sample <- as.factor(Bac.df$Sample)
bac.c <- c("#d63232", "#da4646", "#de5a5a", "#e26f6f", "#e68484", "#ea9898", "#eeadad", "#f2c1c1", "#c02d2d", "#ab2828", "#952323", "#801e1e")

Bac.df$Phylum <- as.factor(Bac.df$Phylum)

bac.df <- data.frame(tax_table(Bacteria))
bac.df %>% dplyr::count(Phylum) -> bac.df
bac.df$Phylum <- as.factor(bac.df$Phylum)
levels(bac.df$Phylum) <- levels(Bac.df$Phylum)


B1 <- ggplot(bac.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
    coord_polar()+
        scale_fill_manual(values=bac.c)+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")

B <- ggplot(Bac.df, aes(x=Sample2, y=Abundance, fill=Phylum))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=bac.c)+
    labs(x="Samples", y="Relative abundance")+
    theme_classic()+
    guides(fill="none")+
    theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


age.df <- PMS@sam_data %>%
    group_by(hyena_ID) %>%
    mutate(min_age=min(age_sampling),
           n= n(),
           max_age=max(age_sampling))

# age at sampling and repeated samples
PMS@sam_data$TimeP <- 0
PMS@sam_data$TimeP[age.df$age_sampling==age.df$min_age] <- 1
PMS@sam_data$TimeP[age.df$age_sampling>age.df$min_age] <- 3
PMS@sam_data$TimeP[age.df$age_sampling==age.df$max_age] <- 2

repeatedS <- ggplot(data=sample_data(PMS),
                    aes(y=age_sampling/365,
                        x=reorder(hyena_ID, age_sampling),
                        group=hyena_ID,
                        fill=age_sampling_cat))+
    geom_point(aes(colour=age_sampling_cat), pch=16, size=2, alpha=0.8)+
        geom_line(size=0.3,alpha=0.5)+
    geom_hline(yintercept=2, colour="firebrick", size=1, alpha=0.5)+
    xlab("Hyena ID")+
    ylab("Age at sampling (years)")+
    scale_color_manual(values=c("#a5c3ab", "#eac161"))+
    scale_y_continuous(breaks=0:16)+
    theme_classic()+
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black"))+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


Fig1_1 <-  cowplot::plot_grid(E1, B1, ncol=2, labels=c("B", "C"))
Fig1_2 <-  cowplot::plot_grid(E, B, ncol=1, labels=c("D", "E"))

Fig1 <- cowplot::plot_grid(repeatedS, Fig1_1,Fig1_2, ncol=1, rel_heights=c(0.8,1, 1))

#Fig1

ggplot2::ggsave(file="fig/Figure1.pdf", Fig1, width = 180, height = 280, dpi = 300, units="mm")

#############################################################################
############# Modeling ß-diversity
# 3 different ß-diversity measures: Bray-Curtis, Jaccard and Aitchison
# we model 1) all observed taxa - overall model
# and 2) a multivariate model with the dependent variables:
## a) Bacteria component - all taxa from bacteria domain
## b) Fungi component - all taxa from Fungi kingdom
## c) Parasite component - all known eukaryotic parasites of hyenas
 do_Models <- FALSE # all models are saved, set TRUE to rerun all models. Warning: takes a few days
if(do_Models){

    sample_data(PMS)$key <- paste(sample_data(PMS)$hyena_ID,
                              sample_data(PMS)$Sample, sep="_")
    key <- data.frame(ID=sample_data(PMS)$key)
    metadt <- sample_data(PMS)
    metadt$IgAP <- log(metadt$IgA_imputed)
    metadt$MucinP <- log(metadt$mucin_imputed)
    # renaming this now
    sample_names(PMS) <- key$ID

#################### we first prepare the data
    # the overall ß-diversity measures
    ## 1) Jaccard distance
    JACM <- as.matrix(phyloseq::distance(PMS,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
    JACM <- 1-JACM
    jac<-c(as.dist(JACM))

    ## 2) Aitchison distance
    AITM <- as.matrix(vegan::vegdist(PMS@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
    AITM <- 1-AITM
    ait<-c(as.dist(AITM))

    ## 3) Bray distance
    BRAM <- as.matrix(vegan::vegdist(PMS@otu_table,
                                 method="bray"))
    BRAM <- 1-BRAM
    bra<-c(as.dist(BRAM))
    
    #the parasite ß-diversity measures
    ## 4) Jaccard distance Parasites
    JACP <- as.matrix(phyloseq::distance(Parasite,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
    JACP[is.na(JACP)]<- 0 # defining those as 0 distances
    JACP <- 1-JACP
    jacP<-c(as.dist(JACP))

    ## 5) Aitchison distance Parasites
    AITP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
    AITP[is.na(AITP)]<- 0 # defining those as 0 distances
    AITP <- 1-AITP
    aitP<-c(as.dist(AITP))

    ## 6) Bray distance Parasites
    BRAP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                     method="bray"))
    BRAP[is.na(BRAP)]<- 0 # defining those as 0 distances
    BRAP <- 1-BRAP
    braP<-c(as.dist(BRAP))
    
    # the bacteria ß-diversity
    ## 7) Jaccard distance Bacteria
    JACB <- as.matrix(phyloseq::distance(Bacteria,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
    JACB[is.na(JACB)]<- 0 # defining those as 0 distances
    JACB <- 1-JACB
    jacB<-c(as.dist(JACB))

    ## 8) Aitchison distance Bacteria
    AITB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
    AITB[is.na(AITB)]<- 0 # defining those as 0 distances
    AITB <- 1-AITB
    aitB<-c(as.dist(AITB))

    ## 9) Bray distance Bacteria
    BRAB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="bray"))
    BRAB[is.na(BRAB)]<- 0 # defining those as 0 distances
    BRAB <- 1-BRAB
    braB<-c(as.dist(BRAB))

    # ß-diversity estimates for fungi
    ## 7) Jaccard distance Fungi
    JACF <- as.matrix(phyloseq::distance(Fungi,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
    JACF[is.na(JACF)]<- 0 # defining those as 0 distances
    JACF <- 1-JACF
    jacF<-c(as.dist(JACF))

    ## 8) Aitchison distance Fungi
    AITF <- as.matrix(vegan::vegdist(Fungi@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
    AITF[is.na(AITF)]<- 0 # defining those as 0 distances
    AITF <- 1-AITF
    aitF<-c(as.dist(AITF))

    ## 8.1) Bray distance Fungi
    BRAF <- as.matrix(vegan::vegdist(Fungi@otu_table,
                                 method="bray"))
    BRAF[is.na(BRAF)]<- 0 # defining those as 0 distances
    BRAF <- 1-BRAF
    braF<-c(as.dist(BRAF))

    
    # pairwise comparisons of clan (shared or not)
    ClanM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$clan[i]==metadt$clan[j]){
                ClanM[i,j]= "1"
            } else{
                ClanM[i,j]= "0"
            }
        }
    }
    ClanM<-c(ClanM[lower.tri(ClanM)])  


    ## now other distances
    AGEM <- c(dist(metadt[,c("age_sampling")]))
    metadt$date_sampling <- as.Date(metadt$date_sampling,
                                    format='%m/%d/%Y')
SamplingM <- c(dist(metadt$date_sampling))


    RANKM<-c(dist(metadt[,c("CSocialRank")]))
    IGAPM<-c(dist(metadt[,c("IgAP")]))
    MUCPM<-c(dist(metadt[,c("MucinP")]))

metadt$Season <- as.numeric(format(metadt$date_sampling, "%m"))
metadt$season <- "wet"
metadt$season[metadt$Season>5&metadt$Season<11] <- "dry"
#table(as.factor(metadt$Season), metadt$season)

    # pairwise comparisons of clan (shared or not)
    SeasonM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$season[i]==metadt$season[j]){
                SeasonM[i,j]= "1"
            } else{
                SeasonM[i,j]= "0"
            }
        }
    }
    SeasonM<-c(SeasonM[lower.tri(SeasonM)])  


    ## genetic mum
    MUM_frame<-metadt[,c("key","genetic_mum")]
    MUMM<-array(0,c(nrow(MUM_frame),nrow(MUM_frame)))
    for(i in 1:nrow(MUM_frame)){
        for(j in 1:nrow(MUM_frame)){
            if(MUM_frame$genetic_mum[i]==MUM_frame$genetic_mum[j]){
                MUMM[i,j]= "1"
            } else{
                MUMM[i,j]= "0"
            }
        }
    }
    MUMM <- as.character(c(as.dist(MUMM)))

    # and we control for batch effects
    BatchM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$batch[i]==metadt$batch[j]){
                BatchM[i,j]= "1"
            } else{
                BatchM[i,j]= "0"
            }
        }
    }
    BatchM<-c(BatchM[lower.tri(BatchM)])  


    
    ##Combine these vectors into a data frame
    data.dyad<-data.frame(MS_J=jac,
                          MS_A=ait,
                          MS_B=bra,
                          Parasite_J=jacP,
                          Parasite_A=aitP,
                          Parasite_B=braP,
                          Bacteria_J=jacB,
                          Bacteria_A=aitB,
                          Bacteria_B=braB,
                          Fungi_J=jacF,
                          Fungi_A=aitF,
                          Fungi_B=braF,
                          Age=AGEM,
                          IgAP=IGAPM,
                          MucinP=MUCPM,
                          Rank=RANKM,
                          Gen_mum=MUMM,
                          Temporal=SamplingM,
                          Clan=ClanM,
                          Season=SeasonM,
                          Batch=BatchM)

    list<-expand.grid(key$ID, key$ID)
    ## This created individual-to-same-individual pairs as well. Get rid of these:
    list<-list[which(list$Var1!=list$Var2),]
    ## this still has both quantiles in--> add 'unique' key
    list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
    list<-subset(list, !duplicated(list$key))
    ## sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
    i=nrow(key)
    JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
    ## add the names of both individuals participating in each dyad into the data frame
    data.dyad$IDA_s<-list$Var2
    data.dyad$IDB_s<-list$Var1
    data.dyad$IDA <- gsub("_.*", "", data.dyad$IDA_s)
    data.dyad$IDB <- gsub("_.*", "", data.dyad$IDB_s)

    ## We keep intra individual comparisons (not the same sample!!)
    ##data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]


    ##scale all predictors to range between 0-1 if they are not already naturally on that scale
    ##define scaling function:
    range.use <- function(x, min.use, max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    scalecols<-c("Rank",
                 "Age",
                 "IgAP",
                 "MucinP",
                 "Parasite_A",
                 "Bacteria_A",
                 "Parasite_A",
                 "Fungi_A",
                 "MS_A",
                 "MS_J",
                 "Temporal")
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-
            range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
    
    saveRDS(data.dyad, "tmp/data.dyad.rds")
    

#### Modeling all taxa
    modelJ<-brm(MS_J~ 1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                    Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
    saveRDS(modelJ, "tmp//modelJ.rds")

    modelA<-brm(MS_A~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                    Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
    saveRDS(modelA, "tmp/modelA.rds")

    modelB<-brm(MS_B~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                    Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
    saveRDS(modelB, "tmp/modelB.rds")
    
######################### Now we model the data ####################
### modeling bacteria, fungi and parasites separately

    modelB_P<-brm(Parasite_B~ 1+ Fungi_B+Bacteria_B+IgAP+ MucinP + Age+  Rank+ Gen_mum+
                   Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                      (1|mm(IDA,IDB)),
                  data = data.dyad,
                  family= "gaussian",
                  warmup = 1000, iter = 3000,
                  cores = 20, chains = 4,
               inits=0)

saveRDS(modelB_P, "tmp/modelB_P.rds")

modelB_F<-brm(Fungi_B~ 1+ Bacteria_B+Parasite_B+IgAP+ MucinP + Age+  Rank+ Gen_mum+
                  Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                      (1|mm(IDA,IDB)),
                  data = data.dyad,
                  family= "gaussian",
                  warmup = 1000, iter = 3000,
                  cores = 20, chains = 4,
              inits=0)

saveRDS(modelB_F, "tmp/modelB_F.rds")

modelB_B<-brm(Bacteria_B~ 1+ Fungi_B+Parasite_B+IgAP+ MucinP + Age+  Rank+ Gen_mum+
                  Temporal+ Clan+Season+Batch+Age:IgAP+Age:MucinP+
                      (1|mm(IDA,IDB)),
                  data = data.dyad,
                  family= "gaussian",
                  warmup = 1000, iter = 3000,
                  cores = 20, chains = 4,
              inits=0)

saveRDS(modelB_B, "tmp/modelB_B.rds")

saveRDS(modelB, "tmp/modelB.rds")
} else 
    data.dyad <- readRDS("tmp/data.dyad.rds")


########################################################################
##### now we test sex effects
##### because all adults are females we do this only for juveniles
doJuv <- FALSE # models are saved, set to TRUE for rerunning the models
if(doJuv){
    Juvenile <- subset_samples(PMS, age_sampling_cat=="sampled_as_juvenile")
    Para <- subset_taxa(Juvenile, Genus%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))
    
##  only bacteria
    Bac <- subset_taxa(Juvenile, Kingdom %in%"Bacteria") 

# only fungi
    Fun <- subset_taxa(Juvenile, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))

    sample_data(Juvenile)$key <- paste(sample_data(Juvenile)$hyena_ID,
                                   sample_data(Juvenile)$Sample, sep="_")

    key <- data.frame(ID=sample_data(Juvenile)$key)
    metadt <- sample_data(Juvenile)
    metadt$IgAP <- log(metadt$IgA_imputed)
    metadt$MucinP <- log(metadt$mucin_imputed)
    # renaming this now
    sample_names(Juvenile) <- key$ID
    
    ## 1) Bray distance for all taxa
    BRAM <- as.matrix(vegan::vegdist(Juvenile@otu_table,
                                 method="bray"))
    BRAM <- 1-BRAM
    bra<-c(as.dist(BRAM))
    
    ## 2) Bray distance Parasites
    BRAP <- as.matrix(vegan::vegdist(Para@otu_table,
                                     method="bray"))
    BRAP[is.na(BRAP)]<- 0 # defining those as 0 distances
    BRAP <- 1-BRAP
    braP<-c(as.dist(BRAP))
    
    ## 3) Bray distance Bacteria
    BRAB <- as.matrix(vegan::vegdist(Bac@otu_table,
                                 method="bray"))
    BRAB[is.na(BRAB)]<- 0 # defining those as 0 distances
    BRAB <- 1-BRAB
    braB<-c(as.dist(BRAB))

    ## 4) Bray distance Fungi
    BRAF <- as.matrix(vegan::vegdist(Fun@otu_table,
                                 method="bray"))
    BRAF[is.na(BRAF)]<- 0 # defining those as 0 distances
    BRAF <- 1-BRAF
    braF<-c(as.dist(BRAF))

    ## 5) Sex pairs
    Sex_frame<-metadt[,c("key","sex")]
    Sex_frame$key<-as.character(Sex_frame$key)
    Sex_frame$sex<-as.character(Sex_frame$sex)
    #Create an empty character matrix to fill with characters
    SEXM<-array(as.character(NA),c(nrow(Sex_frame),nrow(Sex_frame)))
    for(i in 1:nrow(Sex_frame)){
        for(j in 1:nrow(Sex_frame)){
            if(Sex_frame$sex[i]=="F" & Sex_frame$sex[i]==Sex_frame$sex[j]){
                SEXM[i,j]= "FF"}
            if(Sex_frame$sex[i]=="M" & Sex_frame$sex[i]==Sex_frame$sex[j]){
                SEXM[i,j]= "MM"}
            if( Sex_frame$sex[i]!=Sex_frame$sex[j]){
                SEXM[i,j]= "FM"}
        }
    }
    sex<-c(SEXM[lower.tri(SEXM)])

    ## 6)  pairwise comparisons of clan (shared or not)
    ClanM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$clan[i]==metadt$clan[j]){
                ClanM[i,j]= "1"
            } else{
                ClanM[i,j]= "0"
            }
        }
    }
    ClanM<-c(ClanM[lower.tri(ClanM)])  

    ## 7) now other distances
    AGEM <- c(dist(metadt[,c("age_sampling")]))
    metadt$date_sampling <- as.Date(metadt$date_sampling,
                                    format='%m/%d/%Y')
    SamplingM <- c(dist(metadt$date_sampling))
    RANKM<-c(dist(metadt[,c("CSocialRank")]))
    IGAPM<-c(dist(metadt[,c("IgAP")]))
    MUCPM<-c(dist(metadt[,c("MucinP")]))

    metadt$Season <- as.numeric(format(metadt$date_sampling, "%m"))
    metadt$season <- "wet"
    metadt$season[metadt$Season>5&metadt$Season<11] <- "dry"
    # pairwise comparisons of clan (shared or not)
    SeasonM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$season[i]==metadt$season[j]){
                SeasonM[i,j]= "1"
            } else{
                SeasonM[i,j]= "0"
            }
        }
    }
    SeasonM<-c(SeasonM[lower.tri(SeasonM)])  

    ## genetic mum
    MUM_frame<-metadt[,c("key","genetic_mum")]
    MUMM<-array(0,c(nrow(MUM_frame),nrow(MUM_frame)))
    for(i in 1:nrow(MUM_frame)){
        for(j in 1:nrow(MUM_frame)){
            if(MUM_frame$genetic_mum[i]==MUM_frame$genetic_mum[j]){
                MUMM[i,j]= "1"
            } else{
                MUMM[i,j]= "0"
            }
        }
    }
    MUMM <- as.character(c(as.dist(MUMM)))

    # and we control for batch effects
    BatchM<-array(0,c(nrow(metadt),nrow(metadt)))
    for(i in 1:nrow(metadt)){
        for(j in 1:nrow(metadt)){
            if(metadt$batch[i]==metadt$batch[j]){
                BatchM[i,j]= "1"
            } else{
                BatchM[i,j]= "0"
            }
        }
    }
    BatchM<-c(BatchM[lower.tri(BatchM)])  

    #####Combine these vectors into a data frame
    data.dyad<-data.frame(MS_B=bra,
                          Parasite_B=braP,
                          Bacteria_B=braB,
                          Fungi_B=braF,
                          Sex=sex,
                          Age=AGEM,
                          IgAP=IGAPM,
                          MucinP=MUCPM,
                          Rank=RANKM,
                          Gen_mum=MUMM,
                          Temporal=SamplingM,
                          Clan=ClanM,
                          Season=SeasonM,
                          Batch=BatchM)

    list<-expand.grid(key$ID, key$ID)
    ## This created individual-to-same-individual pairs as well. Get rid of these:
    list<-list[which(list$Var1!=list$Var2),]
    ## this still has both quantiles in--> add 'unique' key
    list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
    list<-subset(list, !duplicated(list$key))
    ## sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
    i=nrow(key)
    BRAM[which(rownames(BRAM)==list$Var1[i]),which(colnames(BRAM)==list$Var2[i])]==bra[i]
    ## add the names of both individuals participating in each dyad into the data frame
    data.dyad$IDA_s<-list$Var2
    data.dyad$IDB_s<-list$Var1
    data.dyad$IDA <- gsub("_.*", "", data.dyad$IDA_s)
    data.dyad$IDB <- gsub("_.*", "", data.dyad$IDB_s)
    ## We keep intra individual comparisons (not the same sample!!)
    ##data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]
    ##scale all predictors to range between 0-1 if they are not already naturally on that scale
    ##define scaling function:
    range.use <- function(x, min.use, max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    scalecols<-c("Rank",
                 "Age",
                 "IgAP",
                 "MucinP",
                 "Temporal")
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-
            range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

    model_Juv<-brm(MS_B~1+ Sex+ IgAP+ MucinP + Age+
                    Temporal+ Clan+Season+Batch+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
    saveRDS(model_Juv, "tmp/modelJuv.rds")

        model_JuvP<-brm(Parasite_B~1+ Sex+ IgAP+ MucinP + Age+
                    Temporal+ Clan+Season+Batch+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
    saveRDS(model_JuvP, "tmp/modelJuv_Parasite.rds")

        model_JuvB<-brm(Bacteria_B~1+ Sex+ IgAP+ MucinP + Age+
                    Temporal+ Clan+Season+Batch+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
    saveRDS(model_JuvB, "tmp/modelJuv_Bacteria.rds")

        model_JuvF<-brm(Fungi_B~1+ Sex+ IgAP+ MucinP + Age+
                    Temporal+ Clan+Season+Batch+
                    (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
    saveRDS(model_JuvF, "tmp/modelJuv_Fungi.rds")
}

### now we load all the models
model_Juv <- readRDS("tmp/modelJuv.rds")
model_JuvF <- readRDS("tmp/modelJuv_Fungi.rds")
model_JuvB <- readRDS("tmp/modelJuv_Bacteria.rds")
model_JuvP <- readRDS("tmp/modelJuv_Parasite.rds")

modelA <- readRDS("tmp/modelA.rds")
modelJ <- readRDS("tmp/modelJ.rds")
modelB <- readRDS("tmp/modelB.rds")

modelB_P <- readRDS("tmp/modelB_P.rds")
modelB_B <- readRDS("tmp/modelB_B.rds")
modelB_F <- readRDS("tmp/modelB_F.rds")


## summaries of overall models (3 different metrics of beta diversity)
print(summary(modelB), digits=3)
print(summary(modelA), digits=4)
print(summary(modelJ), digits=3)

print(summary(model_Juv), digits=3)

print(summary(model_JuvF), digits=3)

print(summary(model_JuvP), digits=3)

print(summary(model_JuvB), digits=3)

## ther results are largely similar regardless of the distances used.
## we proceed the BC, a widely used distance in microbiome studies

## summaries for bacteria, parasite and fungi models
print(summary(modelB_P), digits=3)

print(summary(modelB_B), digits=3)

print(summary(modelB_F), digits=3)

print(summary(modelB), digits=3)

#bayes_R2(modelB_FPB)

bayes_R2(modelB_B)

bayes_R2(modelB_P)

bayes_R2(modelB_F)

bayes_R2(model_Juv)

bayes_R2(modelB)

library(ggmcmc)
library(ggridges)

### saving juvenile results for all models
para <- summary(model_Juv)$fixed
res.df <- data.frame(Effect=rownames(para),
           Estimate=para$Estimate,
           lCI=para$'l-95% CI',
           uCI=para$'u-95% CI')

res.df$Domain <- "Overall"
res.df <- res.df[!res.df$Effect=="Intercept",] # removing intercepts

res2 <- summary(model_JuvP)$fixed
res2$Domain <- "Parasite"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

res2 <- summary(model_JuvF)$fixed
res2$Domain <- "Fungi"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

res2 <- summary(model_JuvB)$fixed
res2$Domain <- "Bacteria"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

write.csv(res.df, "tmp/summary_BRMSmodelsJUVENILES.csv")

### saving overall results for all ß-diversity metrics
para <- summary(modelA)$fixed
res.df <- data.frame(Effect=rownames(para),
           Estimate=para$Estimate,
           lCI=para$'l-95% CI',
           uCI=para$'u-95% CI')

res.df$Domain <- "Aitchison"
res.df <- res.df[!res.df$Effect=="Intercept",] # removing intercepts

res2 <- summary(modelB)$fixed
res2$Domain <- "Bray"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

res2 <- summary(modelJ)$fixed
res2$Domain <- "Jaccard"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

res.df

write.csv(res.df, "tmp/summary_BRMSmodelsOVERALL.csv")

############################# plot Figure 2
#the effects of all variables on the overall microbiome
modelB_transformed <- ggs(modelB)

cat <- filter(modelB_transformed, Parameter %in% c("b_IgAP", "b_MucinP", "b_Age", "b_Rank", "b_Gen_mum1", "b_Temporal", "b_Clan1", "b_Season1", "b_IgAP:Age", "b_MucinP:Age"))
# let's not forget to subset here
cat2 <- cat[cat$Iteration>1000,]

cat2$Parameter <- droplevels(cat2$Parameter)

unique(cat2$Parameter)

cat2$Parameter <-factor(cat2$Parameter, levels=c("b_Rank", "b_Season1","b_Temporal", "b_Clan1", "b_Gen_mum1", "b_IgAP:Age", "b_MucinP:Age", "b_IgAP", "b_MucinP", "b_Age"))

colours <- c("#7d1a2a", "#74b2b4", "#bf7e7e", "#bd5569", "#60776a", "#e5cbd9", "#c1d2a2", "#b1658d", "#8e9f6f", "#d8b255") 

Fig2 <- ggplot(cat2, aes(x = value, y=Parameter, fill=Parameter))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
    #    facet_wrap(~Parameter)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA, width=0.2)+
 #   scale_fill_manual(values=c("#fc9c24", "#fce69e", "#24acb4", "#b4ccbc"))+
    scale_fill_manual(values=colours)+
    xlab("Posterior probability distribution")+
    ylab("")+
    theme_classic()+
    theme(legend.position="none")

Fig2

ggplot2::ggsave(file="fig/Figure2.pdf", Fig2, width = 180, height = 150, dpi = 300, units="mm")


##################### Figure 3
##################### A
coul= c("#d63232", "#ffbf00","#293885","#374f2f")

### A plotting the effects of IgA and mucin accross microbiome components
modelB_B3 <- ggs(modelB_B)
catB <- filter(modelB_B3, Parameter %in% c("b_IgAP", "b_MucinP"))
# let's not forget to subset here
cat2B <- catB[catB$Iteration>1000,]
cat2B$Parameter <- droplevels(cat2B$Parameter)
unique(cat2B$Parameter)
cat2B$Domain <- "Bacteria"


modelB_P3 <- ggs(modelB_P)
catP <- filter(modelB_P3, Parameter %in% c("b_IgAP", "b_MucinP"))
# let's not forget to subset here
cat2P <- catP[catP$Iteration>1000,]
cat2P$Parameter <- droplevels(cat2P$Parameter)
unique(cat2B$Parameter)
cat2P$Domain <- "Parasite"

modelB_F3 <- ggs(modelB_F)
catF <- filter(modelB_F3, Parameter %in% c("b_IgAP", "b_MucinP"))
# let's not forget to subset here
cat2F <- catF[catF$Iteration>1000,]
cat2F$Parameter <- droplevels(cat2F$Parameter)
unique(cat2F$Parameter)
cat2F$Domain <- "Fungi"

modelB_3 <- ggs(modelB)
cat <- filter(modelB_3, Parameter %in% c("b_IgAP", "b_MucinP"))
# let's not forget to subset here
cat2 <- cat[cat$Iteration>1000,]
cat2$Parameter <- droplevels(cat2$Parameter)
unique(cat2$Parameter)
cat2$Domain <- "Overall"

cat.df <- rbind(cat2B, cat2P, cat2F, cat2)
cat.df$Domain <- factor(cat.df$Domain, levels=c("Bacteria", "Fungi", "Parasite", "Overall"))

cat.df_A <- cat.df[cat.df$Parameter=="b_IgAP",]

IgAB <- ggplot(cat.df_A, aes(x = value, y=Domain, fill=Domain))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
    #    facet_wrap(~Parameter)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA, width=0.2)+
 #   scale_fill_manual(values=c("#fc9c24", "#fce69e", "#24acb4", "#b4ccbc"))+
    scale_fill_manual(values=coul)+
    xlab("effect size estimate of f-IgA distances")+
    ylab("")+
    theme_classic()+
    theme(legend.position="none")

cat.df_M <- cat.df[cat.df$Parameter=="b_MucinP",]

MucinB <- ggplot(cat.df_M, aes(x = value, y=Domain, fill=Domain))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
    #    facet_wrap(~Parameter)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA, width=0.2)+
 #   scale_fill_manual(values=c("#fc9c24", "#fce69e", "#24acb4", "#b4ccbc"))+
    scale_fill_manual(values=coul)+
    xlab("effect size estimate of f-mucin distances")+
    ylab("")+
    theme_classic()+
    theme(legend.position="none")

Fig3A <- cowplot::plot_grid(IgAB, MucinB, ncol=2)

Fig3A

# B predicted effects of the interaction between immune responses and age
library("interactions")


MucinPara <- interact_plot(modelB_P,
                           pred=MucinP,
                           modx=Age,
                           interval=TRUE,
                           int.type=c("confidence"),
                           int.width=0.5
                           )+
    theme_classic()+
    theme(legend.position="none")

IgAInt <- interact_plot(modelB,
                         pred=IgAP,
                         modx=Age,
                         interval=TRUE,
                         int.type=c("confidence"),
                         int.width=0.5
                                        #                         rug=TRUE,
#                         partial.residuals=TRUE
                        )+
    scale_y_continuous(limits=c(0.20, 0.40))+
    theme_classic()+
    theme(legend.position="none")

MucinInt <- interact_plot(modelB,
                           pred=MucinP,
                           modx=Age,
                           interval=TRUE,
                           int.type=c("confidence"),
                           int.width=0.5
                          )+
    scale_y_continuous(limits=c(0.20, 0.40))+
    theme_classic()+
    theme(legend.position="none")



Fig3B <- cowplot::plot_grid(IgAInt, MucinInt, ncol=2)
Figure3 <- cowplot::plot_grid(Fig3A, Fig3B, ncol=1, labels="AUTO")

ggplot2::ggsave(file="fig/Figure3.pdf", Figure3, width = 180, height = 200, dpi = 300, units="mm")

#############################
# random forest regressions
library(caret)
library(ranger)
library(pdp)
library(patchwork)

PMS.t <- transform(PMS, "compositional") # transform to make all taxa from 0 to 1
otu <- PMS.t@otu_table
colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS.t@tax_table[,6], sep="_")
metadt <- as.data.frame(PMS@sam_data)
metadt$logIgA <- log(metadt$IgA_imputed)
metadt$logmucin <- log(metadt$mucin_imputed)


metadt$date_sampling <- as.Date(metadt$date_sampling,
                                    format='%m/%d/%Y')
metadt$Season <- as.numeric(format(metadt$date_sampling, "%m"))
metadt$season <- "wet"
metadt$season[metadt$Season>5&metadt$Season<11] <- "dry"


## first we want to predict IgA
df <- data.frame(otu, ID=metadt$hyena_ID, Clan=metadt$clan, Age=metadt$age_sampling,
                 CSocialRank=metadt$CSocialRank, Mum=metadt$genetic_mum,
                 Season=metadt$season, batch=metadt$batch, IgA=metadt$logIgA)

df$ID <- as.factor(df$ID)
df$Season <- as.factor(df$Season)
df$batch <- as.factor(df$batch)
df$Mum <- as.factor(df$Mum)
df$Clan <- as.factor(df$Clan)

set.seed(123)
trainIndex <- createDataPartition(df$IgA, p=0.8, list=FALSE, times=1)
IgA_df_train <- df[trainIndex,]
IgA_df_test <- df[-trainIndex,]

tgrid <- expand.grid(
    mtry = 1:ncol(otu),
    splitrule = "variance",
    min.node.size = c(5, 10))


set.seed(123)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

doIgARF <- FALSE

if (doIgARF){
set.seed(1111)
rfFit1 <- train(IgA~., data=IgA_df_train,
                method="ranger",
                trControl=fitControl,
                tuneGrid = tgrid,
                importance="permutation")
saveRDS(rfFit1, "tmp/rfFit_IgA.rds")
} else
    rfFit1 <- readRDS("tmp/rfFit_IgA.rds")


print(rfFit1)
print(rfFit1$finalModel)
test_predictions <- predict(rfFit1, newdata=IgA_df_test)

print(postResample(test_predictions, IgA_df_test$IgA))
pred_obs <- data.frame(predicted=test_predictions, observed=IgA_df_test$IgA)
cor.test(pred_obs$predicted, pred_obs$observed, method="spearman")

#corr <-ggplot(pred_obs, aes(x=predicted, y=observed))+
#        geom_point(size=3, color="black", alpha=0.4)+
#    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
#    stat_poly_line(color="black", size=1, alpha=0.2) +
#    stat_poly_eq() +
#    theme_classic()

imp <- varImp(rfFit1)
imp <- imp$importance
imp$taxa <- rownames(imp)
rownames(imp) <- NULL
imp <- imp[order(-imp$Overall),]
imp20 <- imp[1:20,]
imp20 <- droplevels(imp20)
imp20$taxa <- with(imp20, reorder(taxa, Overall))

IgAt <- data.frame(seq=taxa_names(PMS)[as.numeric(rownames(imp20))], name=imp20$taxa, importance=imp20$Overall)

library(seqinr)
write.csv2(IgAt, "tmp/IgAtop20.csv")
write.fasta(as.list(IgAt$seq), IgAt$name, "tmp/IgAtop20.fasta")

## Figure 4A
topImp <-
    ggplot(imp20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()+
    theme(axis.text.y = element_text(angle = 20, vjust = 0.5))


## Figure S5: plotting marginal effects for supplements
top_features <- imp20$taxa
pd_plots <- list(NULL)
top_features <- as.character(top_features)

for (a in 1:length(top_features)) {
    pd_plots[[a]] <-pdp::partial(rfFit1, pred.var=top_features[a], rug=TRUE)%>%
        autoplot()+
        geom_hline(yintercept = mean(IgA_df_train$IgA), linetype = 2, color = "gray") + 
#            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
        theme(panel.border = element_rect(colour = "black", fill = NA),
                      panel.background = element_blank())
    print(paste0("Partial dependence of ", top_features[a]))
}

figS4_2 <- wrap_plots(pd_plots, ncol=4)

#####################################################
## now we want to predict mucin with the taxa abundances
df <- data.frame(otu, ID=metadt$hyena_ID, Clan=metadt$clan, Age=metadt$age_sampling,
                 CSocialRank=metadt$CSocialRank, Mum=metadt$genetic_mum,
                 Season=metadt$season, batch=metadt$batch, Mucin=metadt$logmucin)

df$ID <- as.factor(df$ID)
df$Season <- as.factor(df$Season)
df$batch <- as.factor(df$batch)
df$Mum <- as.factor(df$Mum)
df$Clan <- as.factor(df$Clan)

set.seed(123)
trainIndex <- createDataPartition(df$Mucin, p=0.8, list=FALSE, times=1)
Mucin_df_train <- df[trainIndex,]
Mucin_df_test <- df[-trainIndex,]

doMucinRF <- FALSE

if (doMucinRF){
set.seed(12)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)
set.seed(123)
rfFitM <- train(Mucin~., data=Mucin_df_train,
                method="ranger",
                tuneGrid = tgrid,
                trControl=fitControl,
                importance="permutation")
saveRDS(rfFitM, "tmp/rfFit_mucin.rds")
} else
    rfFitM <- readRDS("tmp/rfFit_mucin.rds")

#print(rfFitM)
print(rfFitM$finalModel)

test_predictions <- predict(rfFitM, newdata=Mucin_df_test)
print(postResample(test_predictions, Mucin_df_test$Mucin))

pred_obs <- data.frame(predicted=test_predictions, observed=Mucin_df_test$Mucin)

cor.test(pred_obs$predicted, pred_obs$observed, method="spearman")

#corrM <- 
#    ggplot(pred_obs, aes(x=predicted, y=observed))+
#    geom_point(size=3, color="black", alpha=0.4)+
#    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
#    stat_poly_line(color="black", size=1, alpha=0.2) +
#        stat_poly_eq() +
#    stat_smooth(method="lm", formula=y~x, geom="smooth")+
#    theme_classic()

impM <- varImp(rfFitM)
impM <- impM$importance
impM$taxa <- rownames(impM)
rownames(impM) <- NULL
impM <- impM[order(-impM$Overall),]
impM20 <- impM[1:20,]
impM20 <- droplevels(impM20)
impM20$taxa <- with(impM20, reorder(taxa, Overall))

#save to table
mucint <- data.frame(seq=taxa_names(PMS)[as.numeric(rownames(impM20))], name=impM20$taxa, importance=impM20$Overall)

write.csv2(mucint, "tmp/mucintop20.csv")
write.fasta(as.list(mucint$seq), mucint$name, "tmp/mucintop20.fasta")

# fig 4B
topMImp <- ggplot(impM20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()+
    theme(axis.text.y = element_text(angle = 20, vjust = 0.5))



## plotting figure S6 - marginal effects for Mucin
top_features <- impM20$taxa
pd_plots <- list(NULL)
top_features <- as.character(top_features)
for (a in 1:length(top_features)) {
    pd_plots[[a]] <-pdp::partial(rfFitM, pred.var=top_features[a], rug=TRUE)%>%
        autoplot()+
        geom_hline(yintercept = mean(Mucin_df_train$Mucin), linetype = 2, color = "gray") + 
#            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
        theme(panel.border = element_rect(colour = "black", fill = NA),
                      panel.background = element_blank())
    print(paste0("Partial dependence of ", top_features[a]))
}

figS5_2 <- wrap_plots(pd_plots, ncol=4)

ggplot2::ggsave(file="fig/FigurePDP_IgA.pdf", figS4_2, width = 200, height = 180, dpi = 300, units="mm")
ggplot2::ggsave(file="fig/FigurePDP_Muc.pdf", figS5_2, width = 200, height = 180, dpi = 300, units="mm")

#####################################################################
###### network ########################################################
taxa_names(PMS) <- colnames(otu) # renaming the taxa ID

#### spiec easi
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)


############## network estimation ########################
DOnetwork <- FALSE # this tmp file is too heavy to submit to github
                   # if running for the first time, set to true

if (DOnetwork) {
    se.net <- spiec.easi(PMS, method="mb", pulsar.params=pargs)
    saveRDS(se.net, "tmp/se.fnet.rds")
} else
se.net <- readRDS("tmp/se.fnet.rds")

se.net$select$stars$summary # lambda path
se.net$select

#coding bacteria/eukaryote nodes
group <- PMS@tax_table[,1]

group[which(PMS@tax_table[,6]%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))] <- "Parasite"

group[which(PMS@tax_table[,2]%in%c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))] <- "Fungi"

all(rownames(group)==taxa_names(PMS)) # sanity check

Genus <- as.vector(PMS@tax_table[,6])
net.ids <- taxa_names(PMS)

#### plotting
bm=symBeta(getOptBeta(se.net), mode="maxabs")
diag(bm) <- 0
#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort
net <- SpiecEasi::adj2igraph(Matrix::drop0(getRefit(se.net)),
                             edge.attr=list(weight=weights),
                             vertex.attr = list(name=net.ids))

E(net)
betaMat=as.matrix(symBeta(getOptBeta(se.net)))
# we want positive edges to be green and negative to be red
edges <- E(net)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(net, edges[e.index])
    xindex=which(net.ids==adj.nodes[1])
    yindex=which(net.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#435E55FF")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#D64161FF")
    }
}
E(net)$color=edge.colors

### defining attributes
V(net)$group <- as.factor(group)
V(net)$genus <- Genus
V(net)$id <- colnames(otu)
# removing negative values from weight. These are coded by colour
E(net)$weight <- abs(E(net)$weight)

levels(V(net)$group)

#levels=c("Fungi", "Parasite", "Bacteria", "Unknown_domain"))

## curating names
#V(net)$phylum[V(net)$phylum=="Unknown_phylum"] <- "Unknown"
V(net)$genus <- gsub("Unknown_genus_in_", "", V(net)$genus)
V(net)$genus <- gsub("_uncultered", "", V(net)$genus)

deg <- igraph::degree(net)
#label <- V(net)$genus
label <- ""

#label[vip$hub_score<0.10] <- "" # keeping only the labels with high degree
# I also want the names of the taxa that are predicting IgA and mucin
label[which(V(net)$id%in%as.vector(imp20$taxa))] <- V(net)$id[which(V(net)$id%in%as.vector(imp20$taxa))]
label[which(V(net)$id%in%as.vector(impM20$taxa))] <- V(net)$id[which(V(net)$id%in%as.vector(impM20$taxa))]

V(net)$size<-1
V(net)$size[which(V(net)$id%in%as.vector(imp20$taxa))] <- 4
V(net)$size[which(V(net)$id%in%as.vector(impM20$taxa))] <- 4

close <- closeness(net)
betw <- betweenness(net)

# node importance in the network
vip <- data.frame(degree=deg, closeness=close, betweness=betw, genus=V(net)$genus, name=V(net)$id)

vip$hub_score <- hub_score(net)$vector

rownames(vip) <- NULL

all(vip$name==V(net)$id)

vip$IgA <- "no"
vip$IgA[which(V(net)$id%in%as.vector(imp20$taxa))] <- "yes"

vip$muc <- "no"
vip$muc[which(V(net)$id%in%as.vector(impM20$taxa))] <- "yes"

vip


## taxa important for either IgA or mucin
vip$Imm <- vip$IgA
vip$Imm[vip$muc=="yes"] <- "yes"

# testing importance of the predictors in the network
summary(aov(hub_score~Imm, data=vip))

summary(aov(closeness~Imm, data=vip))

summary(aov(betweness~Imm, data=vip))

summary(aov(degree~Imm, data=vip))

## the means and sd
mean(vip$hub_score[vip$Imm=="yes"])
mean(vip$hub_score[vip$Imm=="no"])
sd(vip$hub_score[vip$Imm=="yes"])
sd(vip$hub_score[vip$Imm=="no"])

mean(vip$degree[vip$Imm=="yes"])
mean(vip$degree[vip$Imm=="no"])
sd(vip$degree[vip$Imm=="yes"])
sd(vip$degree[vip$Imm=="no"])

mean(vip$closeness[vip$Imm=="yes"])
mean(vip$closeness[vip$Imm=="no"])
sd(vip$closeness[vip$Imm=="yes"])
sd(vip$closeness[vip$Imm=="no"])

mean(vip$betweness[vip$Imm=="yes"])
mean(vip$betweness[vip$Imm=="no"])
sd(vip$betweness[vip$Imm=="yes"])
sd(vip$betweness[vip$Imm=="no"])



IgA_hub <- ggplot(vip, aes(y=hub_score, x= IgA))+
            geom_boxplot()+
    theme_classic()

ggplot(vip, aes(y= closeness, x= IgA))+
            geom_boxplot()+
            theme_classic()

ggplot(vip, aes(y= betweness, x= IgA))+
            geom_boxplot()+
            theme_classic()


ggplot(vip, aes(y= degree, x= IgA))+
            geom_boxplot()+
            theme_classic()


write.csv(vip, "tmp/Network_VIN.csv")


V(net)$domain <- "circle"
V(net)$domain[V(net)$group=="Bacteria"] <- "circle"
V(net)$domain[V(net)$group=="Fungi"] <- "square"
V(net)$domain[V(net)$group=="Eukarya"] <- "square"
V(net)$domain[V(net)$group=="Parasite"] <- "square"
V(net)$domain[V(net)$group=="Unknown_domain"] <- "square"

# we also want the node color to code for phylum
#nb.col <- length(levels(as.factor(V(net)$group)))
#coul <- colorRampPalette(brewer.pal(8, "Set1"))(nb.col)
#mc <- coul[as.numeric(as.factor(V(net)$group))]

coul= c("#d63232", "#5b5b5b", "#ffbf00","#293885","#b7b7b7")
mc <- coul[as.numeric(V(net)$group)]

pdf("fig/Network.pdf",
    width =50, height = 50)
set.seed(100)
plot(net,
     layout=layout_with_drl,
     vertex.shape=V(net)$domain,
     vertex.label="",
#     vertex.label.dist=0.4,
#     vertex.label.degree=-pi/2,
#     vertex.size=as.integer(cut(hub_score(net)$vector, breaks=10))+1,
     vertex.size=V(net)$size,
     vertex.color=adjustcolor(mc,0.8),
#    edge.width=as.integer(cut(E(net)$weight, breaks=6))/6,
     edge.width=0.5,
     margin=c(0,0,0,0))
legend(x=-2, y=1, legend=levels(V(net)$group), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)
dev.off()

### let's see only the important taxa
idx <- match(unique(c(as.vector(imp20$taxa), as.vector(impM20$taxa))), V(net)$id)

## some taxa names have dots in it and messes up with the match function, so we need to manually adjust
grep("ASV_59_.Eubacterium._nodatum_group", V(net)$id)
idx[10] <- 59
grep("ASV_84_.Eubacterium._brachy_group", V(net)$id) 
idx[23] <- 84
grep("ASV_730_Prevotellaceae_UCG.003", V(net)$id) 
idx[25] <- 730

# remove NA from age
idx <- idx[!is.na(idx)]

g2 <- igraph::induced_subgraph(net, c(idx, neighbors(net, idx)))


V(g2)$group

coul= c("#d63232", "#5b5b5b", "#ffbf00","#293885","#b7b7b7")
mc <- coul[as.numeric(V(g2)$group)]

#V(g2)$label <- V(g2)$name

V(g2)$label <- " "
V(g2)$label[which(V(g2)$id%in%as.vector(imp20$taxa))] <- V(g2)$id[which(V(g2)$id%in%as.vector(imp20$taxa))]
V(g2)$label[which(V(g2)$id%in%as.vector(impM20$taxa))] <- V(g2)$id[which(V(g2)$id%in%as.vector(impM20$taxa))]

pdf("fig/NetworkS.pdf",
    width =8, height = 8)
set.seed(100)
plot(g2,
#     layout=layout.circle,
#     layout=layout.fruchterman.reingold(g2, niter=1000, area=5*vcount(net)^2),
     layout=layout_with_kk,
     vertex.shape=V(g2)$domain,
     vertex.label=V(g2)$label,
     vertex.label.dist=0.4,
     vertex.label.degree=-pi/2,
#     vertex.size=as.integer(cut(hub_score(net)$vector, breaks=10))+1,
     vertex.color=adjustcolor(mc,0.8),
     edge.width=0.5,
     margin=c(0,0,0,0))
legend(x=-2, y=1, legend=levels(V(g2)$group), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)
dev.off()

############### compiling Figure 4 ############################
#### heatmap with correlations
# we first selected the top predictive taxa for both IgA and mucin
ab.df <- as.data.frame(otu[,idx]) # idx was created previously in the network section
ab.df$IgA <- metadt$logIgA
ab.df$mucin <- metadt$logmucin
#ab.df$Sample <- rownames(ab.df)
## correlation
cor.matrix <- cor(ab.df, method="spearman")
# only taxa and IgA/mucin columns
cor_df <- as.data.frame(cor.matrix)[1:32,c(33,34)]
cor_df$taxa <- rownames(cor_df)

## we want the exact same order, which is taxa predicting IgA, then taxa predicting both and finally taxa
# predicting only mucin.
cor_df$taxa <- factor(cor_df$taxa, levels=unique(cor_df$taxa))

library("tidyr")
#library(viridis)
library(scico)


cor_long <- pivot_longer(data=cor_df,
                        cols=-taxa,
                        names_to="Immune_Marker",
                        values_to="Correlation")
## order by ASV number, no longer using
#cor_long <- cor_long %>%
#    mutate(taxa_number =as.numeric(gsub("ASV_([0-9]+).*", "\\1", taxa))) %>%
#    arrange(taxa_number) %>%
#    mutate(taxa=factor(taxa, levels=(unique(taxa))))

heatmap <- ggplot(data=cor_long, mapping=aes(y=Immune_Marker, x=taxa, fill= Correlation))+
    geom_tile()+
    #    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)+
                                        #    scale_fill_viridis(option="C", direction=-1)+
    scale_fill_scico(palette = "vik")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 80, vjust = 0.5)) +
    labs(x="", y="", fill="Spearman's rho")

Fig4A <- cowplot::plot_grid(topImp, topMImp, labels=c("A", "B"))

Fig4ABC <- plot_grid(Fig4A, heatmap, ncol=1, labels=c("", "C"))

ggplot2::ggsave(file="fig/Figure4.pdf", Fig4ABC, width = 200, height = 200, dpi = 300, units="mm")

###########################################################
#### plotting figure S3

plot.df <- as.data.frame(sample_data(PMS))
class(plot.df) <- "data.frame"

plot.df$IgAP <- log(plot.df$IgA_imputed)
plot.df$MucinP <- log(plot.df$mucin_imputed)

plot.df$date_sampling <- as.Date(plot.df$date_sampling,
                                    format='%m/%d/%Y')
plot.df$Year <- as.numeric(format(plot.df$date_sampling, "%Y"))
plot.df$Season <- as.numeric(format(plot.df$date_sampling, "%m"))
plot.df$season <- "wet"
plot.df$season[plot.df$Season>5&plot.df$Season<11] <- "dry"

IgA_hist <- ggplot(plot.df, aes(x=IgAP))+
    geom_histogram()+
    labs(x="Faecal IgA (log scale)", y="Sample count")+
    theme_classic()

mucin_hist <- ggplot(plot.df, aes(x=MucinP))+
    geom_histogram()+
    labs(x="Faecal mucin (log scale)", y="Sample count")+
    theme_classic()

Year_hist <- ggplot(plot.df, aes(x=Year))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Year of sampling", y="Sample count")+
    theme_classic()

rank_hist <- ggplot(plot.df, aes(x=CSocialRank))+
    geom_histogram()+
    labs(x="Social rank", y="Sample count")+
    theme_classic()

season_hist <-    ggplot(plot.df, aes(x=season))+
    geom_histogram(stat="count")+
    labs(x="Season", y="Sample count")+
    theme_classic()



### for these plots we remove the multiple sample
ID.df <- plot.df[,c("clan", "genetic_mum", "hyena_ID")]

ID.df <- unique(ID.df)

clan_hist <- ggplot(ID.df, aes(clan))+
    geom_histogram(stat="count")+
    labs(x="Clan", y="Individual count")+
    theme_classic()

sib.df <- as.data.frame(summary(as.factor(ID.df$genetic_mum)))
names(sib.df) <- "sib"
sib.df$sib <- sib.df$sib-1
sib_hist <- ggplot(sib.df, aes(x=sib))+
    geom_histogram()+
    scale_x_continuous(breaks=seq(0, max(sib.df$sib), by=1))+
    labs(x="Number of siblings", y="Individual count")+
    theme_classic()

FigureS3 <- cowplot::plot_grid(IgA_hist, mucin_hist, Year_hist, rank_hist,
                               season_hist, clan_hist, sib_hist, labels="AUTO", ncol=3)

ggplot2::ggsave(file="fig/FigureS3.pdf", FigureS3, width = 180, height = 230, dpi = 300, units="mm")

#####################################################
                                        #Plotting Figure S4
# ordination plots
#####################################################
## all samples
Overall <- vegan::metaMDS(PMS@otu_table, distance="bray", try=250, k=3)

OverallEFit <- envfit(Overall, plot.df[, c("IgAP", "MucinP", "CSocialRank",
                                           "clan", "Season")], na.rm=TRUE)

data.scores <-  as.data.frame(scores(Overall)$sites)
data.scores <- cbind(data.scores, otu_table(PMS)[rownames(data.scores)])

ALL <- ggplot(data=data.scores, aes(x=NMDS1, y=NMDS2))+
    geom_point(aes(colour=plot.df$age_sampling_cat))+
        scale_color_manual(values=c("#a5c3ab", "#eac161"))+
    geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                 data=scores(OverallEFit, "vectors"), size=1, alpha=0.5)+
    geom_text(data = scores(OverallEFit, "vectors"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(scores(OverallEFit, "vectors")), colour = "red", fontface = "bold")+
    theme_classic()+
    theme(legend.position="none")

## for juveniles
Juvenile <- subset_samples(PMS, age_sampling_cat=="sampled_as_juvenile")
otu_juv <- as.data.frame(Juvenile@otu_table)
class(otu_juv) <- "data.frame"

JUVall <- vegan::metaMDS(otu_juv, distance="bray", try=250, k=3)
juv.df <- plot.df[plot.df$age_sampling_cat=="sampled_as_juvenile",]

JUVallEFit <- envfit(JUVall, juv.df[, c("IgAP", "MucinP", "CSocialRank",
                                           "clan", "Season")], na.rm=TRUE)

juv.scores <-  as.data.frame(scores(JUVall)$sites)

JUV <- ggplot(data=juv.scores, aes(x=NMDS1, y=NMDS2))+
    geom_point(colour="#eac161")+
    geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                 data=scores(JUVallEFit, "vectors"), size=1, alpha=0.5)+
    geom_text(data = scores(JUVallEFit, "vectors"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(scores(JUVallEFit, "vectors")), colour = "red", fontface = "bold")+
    theme_classic()


## for Adults
Adult <- subset_samples(PMS, age_sampling_cat=="sampled_as_adult")
otu_ad <- as.data.frame(Adult@otu_table)
class(otu_ad) <- "data.frame"
ADall <- vegan::metaMDS(otu_ad, distance="bray", try=250, k=3)
ad.df <- plot.df[plot.df$age_sampling_cat=="sampled_as_adult",]

ADallEFit <- envfit(ADall, ad.df[, c("IgAP", "MucinP", "CSocialRank",
                                           "clan", "Season")], na.rm=TRUE)
ad.scores <-  as.data.frame(scores(ADall)$sites)

AD <- ggplot(data=ad.scores, aes(x=NMDS1, y=NMDS2))+
    geom_point(colour="#a5c3ab")+
    geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                 data=scores(ADallEFit, "vectors"), size=1, alpha=0.5)+
    geom_text(data = scores(ADallEFit, "vectors"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(scores(ADallEFit, "vectors")), colour = "red", fontface = "bold")+
    theme_classic()

## For all samples and coloured by year
Year <- ggplot(data=data.scores, aes(x=NMDS1, y=NMDS2))+
    geom_point(aes(colour=as.factor(plot.df$Year)))+
#        scale_color_manual(values=c("#a5c3ab", "#eac161"))+
    geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                 data=scores(OverallEFit, "vectors"), size=1, alpha=0.5)+
    geom_text(data = scores(OverallEFit, "vectors"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(scores(OverallEFit, "vectors")), colour = "red", fontface = "bold")+
    theme_classic()
#    theme(legend.position="none")

FigS4BC <-  cowplot::plot_grid(AD, JUV, ncol=2, labels=c("B", "C"))

FigS4 <-  cowplot::plot_grid(ALL, FigS4BC, Year, ncol=1, labels=c("A", "", "D"))

ggplot2::ggsave(file="fig/FigureS4.pdf", FigS4, width = 180, height = 200, dpi = 350, units="mm")


######################################################
### Individuality of microbiome
library(dplyr)
samplecount <- data.frame(sample_data(PMS)) %>%
    dplyr::group_by(hyena_ID) %>%
    dplyr::summarise(NoSamples=n())

sample_data(PMS)$Scount <- "single"
sample_data(PMS)$Scount[sample_data(PMS)$hyena_ID%in%samplecount$hyena_ID[samplecount$NoSamples>1]] <- "rep"


PMS.r <- subset_samples(PMS, Scount=="rep")
PMS.r@sam_data$age_sampling
PMS.r <- prune_taxa(taxa_sums(PMS.r)>0, PMS.r)
PMS.r

test <- PMS.r@sam_data %>%
    group_by(hyena_ID) %>%
    mutate(min_age=min(age_sampling),
           n= n(),
           max_age=max(age_sampling))

all(test$Sample==PMS.r@sam_data$Sample)

PMS.r@sam_data$TimeP <- 0
PMS.r@sam_data$TimeP[test$age_sampling==test$min_age] <- 1
PMS.r@sam_data$TimeP[test$age_sampling>test$min_age] <- 3
PMS.r@sam_data$TimeP[test$age_sampling==test$max_age] <- 2


temp <- as.data.frame(sample_data(PMS) %>% group_by(hyena_ID) %>%
                              sample_n(1)) #subset phyloseq object to one randomly selected sample 
one.ps <- subset_samples(PMS, Sample%in%temp$Sample)
#ps.core <- core(one.ps, detection = 0, prevalence = .3)
#Core<- phyloseq::prune_taxa(taxa_names(PMS)%in%taxa_names(ps.core), PMS)
#Core <- core(PMS.r, detection = 0, prevalence = .3)


library("GUniFrac", lib="/usr/local/lib/R/site-library")

Parasite.r <- subset_taxa(PMS.r, Genus%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))

##  only bacteria
Bacteria.r <- subset_taxa(PMS.r, Kingdom %in%"Bacteria") 

# only fungi
Fungi.r <- subset_taxa(PMS.r, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))


Bray <- vegdist(PMS.r@otu_table, method="bray")

BRAP <- as.matrix(vegan::vegdist(Parasite.r@otu_table,
                                 method="bray"))
BRAP[is.na(BRAP)]<- 0 # defining those as 0 distances

BRAF <- vegan::vegdist(Fungi.r@otu_table,
                       method="bray")
BRAF[is.na(BRAF)]<- 0 # defining those as 0 distances

BRAB <- as.matrix(vegan::vegdist(Bacteria.r@otu_table,
                                 method="bray"))
BRAB[is.na(BRAB)]<- 0 # defining those as 0 distances

# mucin and IgA ICC

dICC.SE.bt(as.matrix(Bray), strata=PMS.r@sam_data$hyena_ID, B=1000)

dICC.SE.bt(as.matrix(BRAF), strata=Fungi.r@sam_data$hyena_ID, B=1000)

dICC.SE.bt(as.matrix(BRAB), strata=Bacteria.r@sam_data$hyena_ID, B=1000)

dICC.SE.bt(as.matrix(BRAP), strata=Parasite.r@sam_data$hyena_ID, B=1000)

############# let's do dyad comparisons now
sample_data(Core)$key <- paste(sample_data(Core)$hyena_ID,
                              sample_data(Core)$Sample, sep="_")
key <- data.frame(ID=sample_data(Core)$key)
mt <- sample_data(Core)
# renaming this now
sample_names(Core) <- key$ID

####################
## 1) Jaccard distance
BM <- as.matrix(phyloseq::distance(Core,
                                     method="bray",
                                     type="samples"))

# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
BM <- 1-BM
ba<-c(as.dist(BM))
#Combine these vectors into a data frame
data.d<-data.frame(MS_B=ba)


list<-expand.grid(key$ID, key$ID)
list<-list[which(list$Var1!=list$Var2),]
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
i=nrow(key)
BM[which(rownames(BM)==list$Var1[i]),which(colnames(BM)==list$Var2[i])]==ba[i]

# add the names of both individuals participating in each dyad into the data frame
data.d$IDA_s<-list$Var2
data.d$IDB_s<-list$Var1
data.d$IDA <- gsub("_.*", "", data.d$IDA_s)
data.d$IDB <- gsub("_.*", "", data.d$IDB_s)

data.d$Var <- "ABC"
data.d$Var[which(data.d$IDA!=data.d$IDB)] <- "Inter"
data.d$Var[which(data.d$IDA==data.d$IDB)] <- "Intra"

res <- (aov(data=data.d, MS_B~Var))
summary(res)

resA <- (aov(data=data.d, MS_A~Var))
summary(resA)

