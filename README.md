# Microbiome_Hyena

Mucosal immune responses and intestinal microbiome associations in a wild carnivore

## Project description
The mammalian mucosal immunity and gut microbiome are in constant cross-talk. The resulting interactions are influenced by host characteristics and by the host’s biotic and abiotic environments. Knowledge on symbiont-symbiont and symbiont-immune interactions is limited, especially outside the scope of human and model organism studies. Wild populations are exposed to and harbour a much more diverse array of macro- and microorganisms, live in heterogeneous environments and have diverse genetic backgrounds. Here, we explore the associations between two important and broad-acting measures of intestinal mucosal immunity, IgA and mucin, and the gut microbiome, while accounting for host characteristics, social rank, and environmental factors in Serengeti hyenas. 

## Repository structure
data/ --> data necessary to run these scripts; these are phyloseq objects. Each phyloseq object corresponds to OTU table, taxonomic table and metadata for each amplicon, after identification and quality screening of ASVs and taxonomic annotation.

R/ --> scripts

tmp/ --> temporary files created from running the scripts

fig/ --> figures

Briefly:
R/1_preprocess <- this script is where we preprocess the reads. Needs raw files to be reproducible
raw files are found in NCBI SRA project number XXXXX (too large to be stored here), paths need to be adjusted then
taxonomic annotation
removal of contaminants, from here it is reprodocible, temporary file is under tmp/PMS_decontan.rds
Removing ASVs with less than 0.005% abundance and removing samples with less than 100 reads
Removing taxonomic handlers from silva (e.g. g__; s__).
Total sum scaling per amplicon (relative abundances)
Merge into one phyloseq object.
Small adjustments to metadata.
Merge ASVs into cASV that are likely from the same taxa: Co-occurrence networks with only positive edges of ASV abundace per genus. Merge ASVs that cluster together.
final phyloseq object for further analysis is "tmp/fPMS.rds"

R/2_analysis.R -> this script contains all statistical analysis and visualisations for the manuscript. Starting point: "tmp/fPMS.rds"

R/3_FigS2.R -> plotting figure S2.
