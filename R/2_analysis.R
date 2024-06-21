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

PMS <- readRDS("tmp/fPMS.rds")

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
    PMS@sam_data$mucin_inputed<- dfI$mucin
    PMS@sam_data$IgA_inputed <- dfI$IgA
    PMS@sam_data$Clan <- substr(PMS@sam_data$hyena_ID, 1,1)
    saveRDS(PMS, "tmp/PMS_imputed.rds")
} else
    PMS <- readRDS("tmp/PMS_imputed.rds")

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

do_Models <- FALSE

if (do_Models){
#############################################################################
############# let's do dyad comparisons now
sample_data(PMS)$key <- paste(sample_data(PMS)$hyena_ID,
                              sample_data(PMS)$Sample, sep="_")
key <- data.frame(ID=sample_data(PMS)$key)
metadt <- sample_data(PMS)
metadt$IgAP <- log(metadt$IgA_inputed)
metadt$MucinP <- log(metadt$mucin_inputed)
# renaming this now
sample_names(PMS) <- key$ID

####################
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

## 2.1) Bray distance
BRAM <- as.matrix(vegan::vegdist(PMS@otu_table,
                                 method="bray"))
BRAM <- 1-BRAM
bra<-c(as.dist(BRAM))

## 3) Jaccard distance Parasites
JACP <- as.matrix(phyloseq::distance(Parasite,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
JACP[is.na(JACP)]<- 0 # defining those as 0 distances
JACP <- 1-JACP
jacP<-c(as.dist(JACP))

## 4) Aitchison distance Parasites
AITP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITP[is.na(AITP)]<- 0 # defining those as 0 distances
AITP <- 1-AITP
aitP<-c(as.dist(AITP))

## 4.1) Bray distance Parasites
BRAP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="bray"))
BRAP[is.na(BRAP)]<- 0 # defining those as 0 distances
BRAP <- 1-BRAP
braP<-c(as.dist(BRAP))

## 5) Jaccard distance Bacteria
JACB <- as.matrix(phyloseq::distance(Bacteria,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
JACB[is.na(JACB)]<- 0 # defining those as 0 distances
JACB <- 1-JACB
jacB<-c(as.dist(JACB))

## 6) Aitchison distance Bacteria
AITB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITB[is.na(AITB)]<- 0 # defining those as 0 distances
AITB <- 1-AITB
aitB<-c(as.dist(AITB))

## 6.1) Bray distance Bacteria
BRAB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="bray"))
BRAB[is.na(BRAB)]<- 0 # defining those as 0 distances
BRAB <- 1-BRAB
braB<-c(as.dist(BRAB))

    
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

    
## 9) Sex pairs
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

# shared clan
ClanM<-array(0,c(nrow(metadt),nrow(metadt)))
for(i in 1:nrow(metadt)){
    for(j in 1:nrow(metadt)){
        if(metadt$Clan[i]==metadt$Clan[j]){
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

# genetic mum
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

#Combine these vectors into a data frame
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
                      Sex=sex,
                      IgAP=IGAPM,
                      MucinP=MUCPM,
                      Rank=RANKM,
                      Gen_mum=MUMM,
                      Temporal=SamplingM,
                      Clan=ClanM)

list<-expand.grid(key$ID, key$ID)
# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]
# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA_s<-list$Var2
data.dyad$IDB_s<-list$Var1
data.dyad$IDA <- gsub("_.*", "", data.dyad$IDA_s)
data.dyad$IDB <- gsub("_.*", "", data.dyad$IDB_s)

# We keep intra individual comparisons (not the same sample!!)
#data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]

######################### Now we model the data ####################
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
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
    
### Multivariate model for Parasite, Fungi and Bacteria
bformJ <- bf(mvbind(Fungi_J, Parasite_J, Bacteria_J)~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)))+ set_rescor(TRUE)
    
fitJ <- brm(bformJ, data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
    
saveRDS(fitJ, "tmp/modelJ_FPB.rds")

bformA <- bf(mvbind(Fungi_A, Parasite_A, Bacteria_A)~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)))+ set_rescor(TRUE)  
fitA <- brm(bformA, data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(fitA, "tmp/modelA_FPB.rds")

bformB <- bf(mvbind(Fungi_B, Parasite_B, Bacteria_B)~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)))+ set_rescor(TRUE)  
fitB <- brm(bformB, data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(fitB, "tmp/modelB_FPB.rds")

###doing separate the fungi and bacteria model because I can't  plot the interaction otherwise

modelB_P<-brm(Parasite_B~ 1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(modelB_P, "tmp/modelB_P.rds")

modelB_F<-brm(Fungi_B~ 1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(modelB_F, "tmp/modelB_F.rds")

#### model
modelJ<-brm(MS_J~ 1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
    saveRDS(modelJ, "tmp/modelJ.rds")

modelA<-brm(MS_A~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
    saveRDS(modelA, "tmp/modelA.rds")

modelB<-brm(MS_B~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
saveRDS(modelB, "tmp/modelB.rds")
} else
    data.dyad <- readRDS("tmp/data.dyad.rds")


modelJ_FPB <- readRDS("tmp/modelJ_FPB.rds")
modelA_FPB <- readRDS("tmp/modelA_FPB.rds")
modelB_FPB <- readRDS("tmp/modelB_FPB.rds")
modelA <- readRDS("tmp/modelA.rds")
modelJ <- readRDS("tmp/modelJ.rds")
modelB <- readRDS("tmp/modelB.rds")


print(summary(modelB), digits=3)

print(summary(modelJ), digits=3)

print(summary(modelA), digits=3)

print(summary(modelB_FPB), digits=3) # same as individual models but with autcorr between dependent variables

print(summary(modelA_FPB), digits=3) # same as individual models but with autcorr between dependent variables

############################# plotting
### making dataframe for Jaccard first
para <- summary(modelB_FPB)$fixed

library(sjPlot)
tab_model(modelB_FPB, collapse.ci = TRUE, digits=3) # neat function to see the table results, but mixes up the titles of the models, for the multivariate model.
tab_model(modelB, collapse.ci = TRUE, digits=3)

res.df <- data.frame(Effect=rownames(para),
           Estimate=para$Estimate,
           lCI=para$'l-95% CI',
           uCI=para$'u-95% CI')

res.df$Domain <- gsub("B_.*", "", res.df$Effect)
res.df$Effect <-gsub(".*_", "", res.df$Effect)
res.df <- res.df[!res.df$Effect=="Intercept",] # removing intercepts

res2 <- summary(modelB)$fixed
res2$Domain <- "Overall"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]

names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

coul= c( "#ffbf00", "#293885", "#d63232", "#744700")

res.df$Domain <- factor(res.df$Domain, levels=c("Fungi", "Parasite", "Bacteria", "Overall"))

resdf_A <- res.df[res.df$Effect=="IgAP",]
IgA<-ggplot(resdf_A, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-IgA distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

resdf_M <- res.df[res.df$Effect=="MucinP",]
Mucin<-ggplot(resdf_M, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-Mucin distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

resdf_A <- res.df[res.df$Effect=="Age",]
Age<-ggplot(resdf_A, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

### Bray-Curtis distances (abundance-based)
para <- summary(modelB_FPB)$fixed

res.df <- data.frame(Effect=rownames(para),
           Estimate=para$Estimate,
           lCI=para$'l-95% CI',
           uCI=para$'u-95% CI')

res.df$Domain <- gsub("B_.*", "", res.df$Effect)
res.df$Effect <-gsub(".*_", "", res.df$Effect)
res.df <- res.df[!res.df$Effect=="Intercept",] # removing intercepts

res2 <- summary(modelB)$fixed
res2$Domain <- "Overall"
res2$Effect <- rownames(res2)
res2 <- res2[-1,c(9, 1, 3,4,8)]
names(res2) <- names(res.df)
res.df <- rbind(res.df, res2)

res.df$Domain <- factor(res.df$Domain, levels=c("Fungi", "Parasite", "Bacteria", "Overall"))

resdf_B <- res.df[res.df$Effect=="IgAP",]
IgAB<-ggplot(resdf_A, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-IgA distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

resdf_M <- res.df[res.df$Effect=="MucinP",]
MucinB<-ggplot(resdf_M, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-Mucin distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

resdf_A <- res.df[res.df$Effect=="Age",]
AgeB<-ggplot(resdf_A, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


resdf_AA <- res.df[res.df$Effect=="IgAP:Age",]
IgAAgeB<-ggplot(resdf_AA, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age dist:f-IgA dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

resdf_MA <- res.df[res.df$Effect=="MucinP:Age",]
MucAgeB<-ggplot(resdf_MA, aes(x=Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Domain),
                  linewidth=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age dist: f-mucin dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


Figure2 <- plot_grid(IgA, IgAB, Mucin, MucinB, Age, AgeB, labels="AUTO", ncol=2)
ggplot2::ggsave(file="fig/Figure2.pdf", Figure2, width = 190, height = 150, dpi = 300, units="mm")


AgeI <- cowplot::plot_grid(AgeB, IgAAgeB, MucAgeB, labels="AUTO", ncol=3)

###################################### plotting composition
## plotting compostion

## small adjustment here

Eukarya@tax_table[which(is.na(Eukarya@tax_table[,2])),2] <- "Unknown_phylum_in_Eukarya"
Eukarya@tax_table[which(Eukarya@tax_table[,2]=="Unknown_phylum"), 2] <- "Unknown_phylum_in_Eukarya"
Bacteria@tax_table[which(Bacteria@tax_table[,2]=="Unknown_phylum"), 2] <- "Unknown_phylum_in_Bacteria"


library(dplyr)

Euk.df <- data.frame(tax_table(Eukarya))
Euk.df %>% count(Phylum) -> Euk.df


Eukarya.t <- microbiome::transform(Eukarya, "compositional")
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
            "#202c6a", "#293885", "#535f9d"), c("#e5e5e5", "#cccccc", "#b2b2b2", "#999999", "#7f7f7f","#666666", "#bfaf7f", "#c5b78b", "#cbbf98", "#d2c7a5", "#d8cfb2", "#e68bc0", "#e998c7", "#eca5ce", "#eeb2d5", "#966fbb", "#a07dc1", "#997db4", "#ab8bc8", "#b59acf", "#c0a8d6", "#cab7dd"))

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
    guides(fill="none")

Bac.t <- microbiome::transform(Bacteria, "compositional")
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
bac.df %>% count(Phylum) -> bac.df
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
    guides(fill="none")

#library(ggpubr)
#lB <- get_legend(B)
#lE <- get_legend(E)
FigS2 <-  plot_grid(E1, E, B1, B, ncol=2, labels="AUTO")

DE <- plot_grid(IgAB, MucinB, labels=c("E", "F"), nrow=1)

Fig2 <- plot_grid(FigS2, DE, nrow=2, rel_heights=c(1,0.35))

ggplot2::ggsave(file="fig/FigureS_composition.pdf", Fig2, width = 180, height = 250, dpi = 300, units="mm")

# interaction
library(tidyr)
library(modelr)
library(tidybayes)

######################## Figure 3: interaction
newdata1 <- data.frame(IgAP=seq_range(0:1, n=51),
                       MucinP=rep(median(data.dyad$MucinP), n=51),
                       Age=rep(1, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
newdata0 <- data.frame(IgAP=seq_range(0:1, n=51),
                       MucinP=rep(median(data.dyad$MucinP), n=51),
                       Age=rep(0, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
pred.df0 <- add_epred_draws(newdata0, modelB)
pred.df1 <- add_epred_draws(newdata1, modelB)
pred.df <- rbind(pred.df0, pred.df1)
pred.df$Age <- as.factor(pred.df$Age)

IgA_B <- ggplot(pred.df, aes(x=IgAP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Overall microbiome similarity")+
    xlab("f-IgA distances")+
    theme_classic()

newdata1 <- data.frame(MucinP=seq_range(0:1, n=51),
                       IgAP=rep(median(data.dyad$IgAP), n=51),
                       Age=rep(1, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
newdata0 <- data.frame(MucinP=seq_range(0:1, n=51),
                       IgAP=rep(median(data.dyad$IgAP), n=51),
                       Age=rep(0, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
pred.df0 <- add_epred_draws(newdata0, modelB)
pred.df1 <- add_epred_draws(newdata1, modelB)

pred.df <- rbind(pred.df0, pred.df1)
pred.df$Age <- as.factor(pred.df$Age)

Mucin_B <-    ggplot(pred.df, aes(x=MucinP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Overall microbiome similarity")+
    xlab("f-mucin distances")+
    theme_classic()

AgeI2 <- cowplot::plot_grid(IgA_B, Mucin_B, labels=c("D", "E"))

Fig3 <- cowplot::plot_grid(AgeI, AgeI2, rel_heights=c(0.6, 1), ncol=1)

ggplot2::ggsave(file="fig/Figure3.pdf", Fig3, width = 190, height = 150, dpi = 300, units="mm")

#############################
# random forest regression
library(caret)
library(ranger)
library(pdp)
library(patchwork)

PMS.t <- transform(PMS, "compositional") # transform to make all taxa from 0 to 1
# this is to make it easier to compare among taxa
otu <- PMS.t@otu_table

colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS.t@tax_table[,6], sep="_")
df <- data.frame(otu, IgA=PMS.t@sam_data$IgA_inputed)
df$IgA <- log(df$IgA)

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

corr <-ggplot(pred_obs, aes(x=predicted, y=observed))+
        geom_point(size=3, color="black", alpha=0.4)+
#    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
    stat_poly_line(color="black", size=1, alpha=0.2) +
    stat_poly_eq() +
    theme_classic()

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

topImp <-
    ggplot(imp20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()

Fig4 <- cowplot::plot_grid(corr, topImp, labels="AUTO", rel_widths=c(0.6, 1))

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

fig4_2 <- wrap_plots(pd_plots, ncol=4)
fig4 <- cowplot::plot_grid(Fig4, fig4_2, nrow=2, rel_heights=c(0.5, 1))


# for mucin
df <- data.frame(otu, Mucin=PMS.t@sam_data$mucin_inputed)
df$Mucin <- log(df$Mucin)

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

corrM <- 
    ggplot(pred_obs, aes(x=predicted, y=observed))+
    geom_point(size=3, color="black", alpha=0.4)+
#    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
    stat_poly_line(color="black", size=1, alpha=0.2) +
        stat_poly_eq() +
#    stat_smooth(method="lm", formula=y~x, geom="smooth")+
    theme_classic()

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

topMImp <- ggplot(impM20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()

Fig5 <- cowplot::plot_grid(corrM, topMImp, labels=c("C", "D"), rel_widths=c(0.6, 1))

Figure4 <- cowplot::plot_grid(Fig4, Fig5, ncol=1)
ggplot2::ggsave(file="fig/Figure4.pdf", Figure4, width = 200, height = 200, dpi = 300, units="mm")

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

fig5_2 <- wrap_plots(pd_plots, ncol=4)

ggplot2::ggsave(file="fig/FigurePDP_IgA.pdf", fig4_2, width = 200, height = 180, dpi = 300, units="mm")

ggplot2::ggsave(file="fig/FigurePDP_Muc.pdf", fig5_2, width = 200, height = 180, dpi = 300, units="mm")


#####################################################################
###### network

colnames(otu) # this is defined above.

taxa_names(PMS) <- colnames(otu)

## prevalebce filtering of 10%
KeepTaxap <- microbiome::prevalence(PMS)>0.2
PS<- phyloseq::prune_taxa(KeepTaxap, PMS)
PS

#### spiec easi
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)

## mb
#se.net <- spiec.easi(PS, method="mb", pulsar.params=pargs)
#saveRDS(se.net, "tmp/se.fnet.rds")
se.net <- readRDS("tmp/se.fnet.rds")

se.net$select$stars$summary # lambda path
se.net$select

#coding bacteria/eukaryote nodes
group <- PS@tax_table[,1]

group[which(PS@tax_table[,6]%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))] <- "Parasite"

group[which(PS@tax_table[,2]%in%c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))] <- "Fungi"

all(rownames(group)==taxa_names(PS)) # sanity check

Genus <- as.vector(PS@tax_table[,6])

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
V(net)$id <- taxa_names(PS)

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

V(net)$size<-0.5
V(net)$size[which(V(net)$id%in%as.vector(imp20$taxa))] <- 1
V(net)$size[which(V(net)$id%in%as.vector(impM20$taxa))] <- 1

close <- closeness(net)
betw <- betweenness(net)

# node importance in the network
vip <- data.frame(degree=deg, closeness=close, betweness=betw, genus=V(net)$genus, name=rownames(vip))

vip$hub_score <- hub_score(net)$vector

rownames(vip) <- NULL

all(vip$name==V(net)$id)

vip$IgA <- "no"
vip$IgA[which(V(net)$id%in%as.vector(imp20$taxa))] <- "yes"

vip$muc <- "no"
vip$muc[which(V(net)$id%in%as.vector(impM20$taxa))] <- "yes"

vip

# testing importance of the predictors in the network
summary(aov(hub_score~muc, data=vip))
summary(aov(hub_score~IgA, data=vip))

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
    width =10, height = 10)
set.seed(100)
plot(net,
     layout=layout.fruchterman.reingold,
     vertex.shape=V(net)$domain,
     vertex.label=label,
#     vertex.label.dist=0.4,
#     vertex.label.degree=-pi/2,
#     vertex.size=as.integer(cut(hub_score(net)$vector, breaks=10))+1,
     vertex.size=1,
     vertex.color=adjustcolor(mc,0.8),
     edge.width=as.integer(cut(E(net)$weight, breaks=6))/3,
     margin=c(0,0,0,0))
legend(x=-2, y=1, legend=levels(V(net)$group), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)
dev.off()



modules =cluster_fast_greedy(net, weights=E(net)$weight)
modularity(modules)
sizes(modules)

V(net)$cluster=modules$membership

nodes <- V(net)$name

cluster_id <- V(net)$cluster

nodes<-as.data.frame(cbind(nodes, cluster_id, group=as.vector(group), Genus))


nodes[cluster_id=="1",]

names(nodes)

gen.cl <- as.data.frame(table(nodes$group, nodes$cluster_id))


Mod.m <- ggplot(gen.cl) +
    geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', width = 0.5, color="black") +
        labs(x = "Module",
             y = "Count") +
        guides(fill=guide_legend(title="Family/Antibiotic type")) +
#        scale_fill_manual(values =mycolor)+
    theme_bw()
#    guides(fill="none")


#### plotting figure 1
age.df <- PMS@sam_data %>%
    group_by(hyena_ID) %>%
    mutate(min_age=min(age_sampling),
           n= n(),
           max_age=max(age_sampling))

PMS@sam_data$TimeP <- 0
PMS@sam_data$TimeP[age.df$age_sampling==age.df$min_age] <- 1
PMS@sam_data$TimeP[age.df$age_sampling>age.df$min_age] <- 3
PMS@sam_data$TimeP[age.df$age_sampling==age.df$max_age] <- 2

repeatedS <- ggplot(data=sample_data(PMS),
                    aes(x=age_sampling/365,
                        y=reorder(hyena_ID, age_sampling),
                        group=hyena_ID,
                        fill=age_sampling_cat))+
    geom_point(colour="white",pch=21,size=3, alpha=0.7)+
    geom_line(size=0.5,alpha=0.5)+
    geom_vline(xintercept=2, colour="firebrick", size=2, alpha=0.5)+
    ylab("Hyena ID")+
    xlab("Age at sampling (years)")+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"))+
    scale_x_continuous(breaks=0:16)+
    theme_classic()+
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black"))+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())


IgA_hist <- ggplot(sample_data(PMS), aes(x=log(IgA_inputed)))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Faecal IgA (RU, log scale)", y="Sample count")+
    theme_classic()

mucin_hist <- ggplot(sample_data(PMS), aes(x=log(mucin_inputed)))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Faecal mucin (OE, log scale)", y="Sample count")+
    theme_classic()


Age_diff <- ggplot(data.dyad, aes(x=Age))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Age distance", y="Pair count")+
    theme_classic()

IgA_diff <- ggplot(data.dyad, aes(x=IgAP))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="f-IgA distance", y="Pair count")+
    theme_classic()

Muc_diff <- ggplot(data.dyad, aes(x=MucinP))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="f-mucin distance", y="Pair count")+
    theme_classic()

AgeSD <- plot_grid(repeatedS, Age_diff, nrow=1, rel_widths=c(1, 0.5), labels=c("a", "b"))
ImmSD <- plot_grid(IgA_hist, IgA_diff, mucin_hist, Muc_diff, labels=c("c", "d", "e", "f"))
FigureSD <- plot_grid(AgeSD, ImmSD, nrow=2)

ggplot2::ggsave(file="fig/FigureSD.pdf", FigureSD, width = 190, height = 190, dpi = 300, units="mm")


######################################################
### Individual repeatability
library(dplyr)
samplecount <- data.frame(sample_data(PMS)) %>%
    group_by(hyena_ID) %>%
    summarise(NoSamples=n())

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
Core <- core(PMS.r, detection = 0, prevalence = .3)


library("GUniFrac", lib="/usr/local/lib/R/site-library")


jac <- vegdist(Core@otu_table, method="bray")

BRAP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="bray"))
BRAP[is.na(BRAP)]<- 0 # defining those as 0 distances

BRAF <- vegan::vegdist(Fungi@otu_table,
                       method="bray")
BRAF[is.na(BRAF)]<- 0 # defining those as 0 distances

BRAB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="bray"))
BRAB[is.na(BRAB)]<- 0 # defining those as 0 distances

# mucin and IgA ICC

dICC(BRAF, strata=Fungi@sam_data$hyena_ID)
dICC(BRAP, strata=Parasite@sam_data$hyena_ID)
dICC(BRAB, strata=Bacteria@sam_data$hyena_ID)

dICC(jac, strata=Core@sam_data$hyena_ID)



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

