#k represents each of the sub-sampled phenotype 

args <- commandArgs(trailingOnly = TRUE)
#print(args)
k <- args[1]


setwd("path/to/sample_size/qtl2") # Replace with the actual directory path



library(qtl2)
library(qtl2fst) # For memory-intensive calc_genoprob with many SNPs
library(yaml)

CORES <- 8
GRID_SPACING <- 0.1 # Distance between pseudomarkers used for mapping, in centimorgans.


control_filename=paste0("control_",k,".yaml")
cross <- read_cross2(control_filename)
yaml_file=yaml::yaml.load_file(input=control_filename)
pheno_name=gsub(".csv","",yaml_file$pheno)
pr <-readRDS("/path/to/sample_size/qtl2/pr/pr_grid.rds")
## Uncomment this and comment out the above lines to reload previously computed probabilities:
# pr <- readRDS("qtl2/pr/pr_grid.rds")


kinship<-readRDS("/path/to/sample_size/qtl2/kinship/kinship.rds")



cat("Running genome scan...\n")
out <- scan1(pr, cross$pheno[,"bmi_bodylength_w_tail"], kinship, quiet = FALSE, cores = CORES)

out_filename=paste0("/path/to/sample_size/qtl2/out/scan_",k,"_out.rds")
saveRDS(out, out_filename)

gmap <- insert_pseudomarkers(cross$gmap, step = GRID_SPACING)


#out<-readRDS(out_filename)

all_lod<-find_peaks(out, gmap, threshold=2, drop=1.5)

write.csv(all_lod,paste0("/path/to/sample_size/qtl2/peaks/",pheno_name,".csv"),row.names = F,quote=T)





