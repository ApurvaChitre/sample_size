setwd("path/to/sample_size/qtl2") # Replace with the actual directory path

library(qtl2)
library(qtl2fst) # For memory-intensive calc_genoprob with many SNPs

CORES <- 6
GRID_SPACING <- 0.1 # Distance between pseudomarkers used for mapping, in centimorgans.

cross <- read_cross2("control.yaml")

cat("Calculating haplotype probabilities...\n")
gmap <- insert_pseudomarkers(cross$gmap, step = GRID_SPACING)
pr_dirname="path/to/sample_size/qtl2/pr" # Replace with the actual directory path
pr <- calc_genoprob_fst(cross, "pr", pr_dirname, map = gmap, error_prob = 0.01, quiet = FALSE, cores = CORES)
## Uncomment this and comment out the above line to load a previously calculated haplotype probability array:
# pr <- readRDS("qtl2/pr/pr_fstindex.rds")
saveRDS(pr, "path/to/sample_size/qtl2/pr/pr_fstindex.rds") # Replace with the actual directory path


cat("Thinning out dense regions of SNPs...\n")
grid <- calc_grid(cross$gmap, step = GRID_SPACING)
pr <- probs_to_grid(pr, grid)  # Subset genotype probabilities to grid
saveRDS(pr, "path/to/sample_size/qtl2/pr/pr_grid.rds") # Replace with the actual directory path
## Uncomment this and comment out the above lines to reload previously computed probabilities:
# pr <- readRDS("qtl2/pr/pr_grid.rds")

cat("Calculating LOCO kinship matrices...\n")
kinship <- calc_kinship(pr, "loco", quiet = FALSE, cores = CORES)
saveRDS(kinship, "path/to/sample_size/qtl2/kinship/kinship.rds") # Replace with the actual directory path
