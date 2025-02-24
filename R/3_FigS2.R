## in this screen we plot the reads we are removing in the quality filtering
## before we were removing all reads with less than 0.005% of abundance and 1% prevalence
## Most removed reads were actually singletons which also had very low abundance
## We now changed the filtering to remove only singletons
library(phyloseq)
library(ggplot2)
library(dplyr)

PMS.l <- readRDS("tmp/PMS.l_decontan.rds")

## quality filtering of ASVs remove singletons
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # remove singletons (prevalnce filter at 0.5%)
    KeepTaxap <- microbiome::prevalence(ps)>0.005
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
    # subset samples based on total read count (keep samples with more than 100 reads)
    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}


## filtering MA by amplicon
fPMS.l <- list()
for (i in 1:length(PMS.l)) {
    try(fPMS.l[[i]] <- fil(PMS.l[[i]]), silent=TRUE)
}


#### let's normalise by amplicon
TPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(fPMS.l)) {
    try(TPMS.l[[i]] <- transform_sample_counts(fPMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}

# And now merge
TPMS <- TPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(TPMS.l)){
    TPMS <- try(merge_phyloseq(TPMS,TPMS.l[[i]]))
}

TPMS # 4706 out of 199 samples

### singletons
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # keep singletons (prevalnce filter at 0.5%)
    KeepTaxap <- microbiome::prevalence(ps)<0.005
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
    # subset samples based on total read count (keep samples with more than 100 reads)
#    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}

## filtering MA by amplicon
sPMS.l <- list()
for (i in 1:length(PMS.l)) {
    try(sPMS.l[[i]] <- fil(PMS.l[[i]]), silent=TRUE)
}


#### let's normalise by amplicon
SPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(sPMS.l)) {
    try(SPMS.l[[i]] <- transform_sample_counts(sPMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}

# And now merge
SPMS <- SPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(SPMS.l)){
    SPMS <- try(merge_phyloseq(SPMS,SPMS.l[[i]]))
}

SPMS # 16869 out of 211 samples

# then we see the low abundant taxa (that might overlap with singletons)
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) < 0.00005)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
    ps 
}

## filtering MA by amplicon
lPMS.l <- list()
for (i in 1:length(PMS.l)) {
    try(lPMS.l[[i]] <- fil(PMS.l[[i]]), silent=TRUE)
    lPMS.l
}

#### let's normalise by amplicon
LPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(lPMS.l)) {
    try(LPMS.l[[i]] <- transform_sample_counts(lPMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}

# And now merge
L.PMS <- LPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(LPMS.l)){
    L.PMS <- try(merge_phyloseq(L.PMS,LPMS.l[[i]]))
}

L.PMS

SPMS # singletons
L.PMS # low abundant

### let's see the overlap between singletons and low abundant taxa
summary(taxa_names(SPMS)%in%taxa_names(L.PMS)) # it's a massive overlap. only 62 ASVs that are low abundant are not singletons.

Keep <- taxa_names(SPMS)[-which(taxa_names(SPMS)%in%taxa_names(L.PMS))]

Low <- prune_taxa(Keep,SPMS)


# remove those ugly handlers for the singletons too
tax_table(SPMS)[,colnames(tax_table(SPMS))] <- gsub(tax_table(SPMS)[, colnames(tax_table(SPMS))], pattern="[a-z]__", replacement="")
all(sample_names(SPMS)==sample_names(SPMS))
tax <- as.data.frame(tax_table(SPMS))
tax$Kingdom[is.na(tax$Kingdom)] <- "Unknown_domain"
tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Archaea", "Unknown_domain")] <- "Eukarya"
tax$Phylum[is.na(tax$Phylum)] <- "Unknown_phylum"
tax$Class[is.na(tax$Class)] <- "Unknown_class"
tax$Order[is.na(tax$Order)] <- "Unknown_order"
SPMS@tax_table <-tax_table(as.matrix(tax))


###################################### plotting composition of singletons
## small adjustment here

Sin.df <- data.frame(tax_table(SPMS))
Sin.df %>% dplyr::count(Phylum) -> Sin.df

Sin1 <- ggplot(Sin.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
    coord_polar()+
#        scale_fill_manual(values=coul)+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")
Sin1

ggplot2::ggsave(file="fig/FigureS2_compositionSing.pdf", Sin1, width = 180, height = 180, dpi = 300, units="mm")
