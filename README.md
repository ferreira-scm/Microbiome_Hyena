# Microbiome_Hyena

Intestinal crosstalks: associations between mucosal immunity and microbiome in free-ranging spotted hyenas
Soares et al 2024

data/ --> data necessary to run these scripts; these are phyloseq objects. Each phyloseq object corresponds to OTU table, taxonomic table and metadata for each amplicon, after identification and quality screening of ASVs and taxonomic annotation.

R/ --> scripts


tmp/ --> temporary files created from running the scripts

fig/ --> figures

Briefly:
R/1_preprocess <- this script is where we preprocess the reads. Needs raw files to be reproducible
raw files are found in NCBI SRA project number XXXXX (too large to be stores here), paths need to be adjusted then
Merge ASVs into cASV that are likely from the same taxa: Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.
taxonomic annotation
removal of contaminants, from here it is reprodocible, temporary file is under tmp/PMS_decontan.rds
Removing ASVs with less than 0.005% abundance and removing samples with less than 100 reads
Removing taxonomic handlers from silva (e.g. g__; s__).
Total sum scaling per amplicon (relative abundances)
Merge into one phyloseq object.
Small adjustments to metadata.
Merge ASVs into cASV that are likely from the same taxa: Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.
final phyloseq object for further analysis is "tmp/fPMS.rds"
