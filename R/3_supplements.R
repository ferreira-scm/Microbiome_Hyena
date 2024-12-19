## in this screen we plot the reads we are removing in the quality filtering
## before we were removing all reads with less than 0.005% of abundance and 1% prevalence
## Most removed reads were actually singletons which also had very low abundance
## We now changed the filtering to remove only singletons
library(phyloseq)
library(ggplot2)
library(dplyr)

PMS.l <- readRDS("tmp/PMS.l_decontan.rds")

## quality filtering 1% prevalence and 0.005% relative abundance plus samples with less than 100 reads
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
#    summary(keepTaxa)
   ps = phyloseq::prune_taxa(keepTaxa, ps)
    # plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
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

# now I want to see the singletons which corresponds to prevalence less than 1%
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
#    summary(keepTaxa)
    ps.f = phyloseq::prune_taxa(keepTaxa, ps)
    ps.low <- phyloseq::prune_taxa(!keepTaxa, ps)
    # plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps.f)>0.01
    ps.f <- phyloseq::prune_taxa(KeepTaxap, ps.f)
    KeepTaxap_S <- microbiome::prevalence(ps)>0.01
    ps.S <- phyloseq::prune_taxa(!KeepTaxap_S, ps)
    # subset samples based on total read count (keep samples with more than 100 reads)
    ps.f <- phyloseq::prune_samples(sample_sums(ps.f)>100, ps.f)
 #   ps.f # filtered
    ps.S # all singletons
#    ps.low # all low abundant taxa
}

## filtering MA by amplicon
singPMS.l <- list()
for (i in 1:length(PMS.l)) {
    try(singPMS.l[[i]] <- fil(PMS.l[[i]]), silent=TRUE)
}
#### let's normalise by amplicon
TsingPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(singPMS.l)) {
    try(TsingPMS.l[[i]] <- transform_sample_counts(singPMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}
# And now merge
TsingPMS <- TsingPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(TsingPMS.l)){
    TsingPMS <- try(merge_phyloseq(TsingPMS,TsingPMS.l[[i]]))
}

# then we see the low abundant taxa (that might overlap with singletons)
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
#    summary(keepTaxa)
    ps.f = phyloseq::prune_taxa(keepTaxa, ps)
    ps.low <- phyloseq::prune_taxa(!keepTaxa, ps)
    # plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps.f)>0.01
    ps.f <- phyloseq::prune_taxa(KeepTaxap, ps.f)
    KeepTaxap_S <- microbiome::prevalence(ps)>0.01
    ps.S <- phyloseq::prune_taxa(!KeepTaxap_S, ps)
    # subset samples based on total read count (keep samples with more than 100 reads)
    ps.f <- phyloseq::prune_samples(sample_sums(ps.f)>100, ps.f)
 #   ps.f # filtered
    #ps.S # all singletons
    ps.low # all low abundant taxa
}
## filtering MA by amplicon
lowPMS.l <- list()
for (i in 1:length(PMS.l)) {
    try(lowPMS.l[[i]] <- fil(PMS.l[[i]]), silent=TRUE)
}
#### let's normalise by amplicon
TlowPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(lowPMS.l)) {
    try(TlowPMS.l[[i]] <- transform_sample_counts(lowPMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}
# And now merge
TlowPMS <- TlowPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(TlowPMS.l)){
    TlowPMS <- try(merge_phyloseq(TlowPMS,TlowPMS.l[[i]]))
}

### and now without filtering
#### let's normalise by amplicon
nfPMS.l <- list() # a list of normalised amplicons
for (i in 1:length(PMS.l)) {
    try(nfPMS.l[[i]] <- transform_sample_counts(PMS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}
# And now merge
nfPMS <- nfPMS.l[[1]] # one phyloseq objects with all the amplicons normalised individually
for (i in 2:length(nfPMS.l)){
    nfPMS <- try(merge_phyloseq(nfPMS,nfPMS.l[[i]]))
}


nfPMS
TPMS
TlowPMS
TsingPMS

### let's see the overlap between singletons and low abundant taxa
summary(taxa_names(TlowPMS)%in%taxa_names(TsingPMS)) # it's a massive overlap. only 62 ASVs that are low abundant are not singletons.

Keep <- taxa_names(TlowPMS)[-which(taxa_names(TlowPMS)%in%taxa_names(TsingPMS))]

Low <- prune_taxa(Keep, TlowPMS)


# remove those ugly handlers
tax_table(TlowPMS)[,colnames(tax_table(TlowPMS))] <- gsub(tax_table(TlowPMS)[, colnames(tax_table(TlowPMS))], pattern="[a-z]__", replacement="")
all(sample_names(TlowPMS)==sample_names(TlowPMS))
tax <- as.data.frame(tax_table(TlowPMS))
tax$Kingdom[is.na(tax$Kingdom)] <- "Unknown_domain"
tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Archaea", "Unknown_domain")] <- "Eukarya"
tax$Phylum[is.na(tax$Phylum)] <- "Unknown_phylum"
tax$Class[is.na(tax$Class)] <- "Unknown_class"
tax$Order[is.na(tax$Order)] <- "Unknown_order"
TlowPMS@tax_table <-tax_table(as.matrix(tax))

# remove those ugly handlers for the singletons too
tax_table(TsingPMS)[,colnames(tax_table(TsingPMS))] <- gsub(tax_table(TsingPMS)[, colnames(tax_table(TsingPMS))], pattern="[a-z]__", replacement="")
all(sample_names(TsingPMS)==sample_names(TsingPMS))
tax <- as.data.frame(tax_table(TsingPMS))
tax$Kingdom[is.na(tax$Kingdom)] <- "Unknown_domain"
tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Archaea", "Unknown_domain")] <- "Eukarya"
tax$Phylum[is.na(tax$Phylum)] <- "Unknown_phylum"
tax$Class[is.na(tax$Class)] <- "Unknown_class"
tax$Order[is.na(tax$Order)] <- "Unknown_order"
TsingPMS@tax_table <-tax_table(as.matrix(tax))


###################################### plotting composition
## small adjustment here

Sin.df <- data.frame(tax_table(TsingPMS))
Sin.df %>% count(Phylum) -> Sin.df

Sin1 <- ggplot(Sin.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
    coord_polar()+
#        scale_fill_manual(values=coul)+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")
Sin1

ggplot2::ggsave(file="fig/FigureS2_compositionSing.pdf", Sin1, width = 180, height = 180, dpi = 300, units="mm")
