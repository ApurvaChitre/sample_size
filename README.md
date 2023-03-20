
# Prerequisites

The code was run on a shared high-performance computing cluster (HPC system) with a PBS Resource Manager. On this system, we had access to four 2 x 12-core Intel Haswell processors with 64GB of main memory.  

Each standalone script of the pipeline can be adapted for other systems with adequate adjustments such as inputting the paths for the locations of files/programs and configuring appropriate system-specific scheduler directives.  

**Input files for GWAS analysis that utilized SNPs**: The Genotype, Phenotype and the Genetic Relatedness Matrix (GRM) files can be accessed from the __UC San Diego Library Digital Collections repository__ at https://library.ucsd.edu/dc/object/bb9156620z  

**Input files required for linkage mapping using haplotypes** are present in this __Zenodo repository__.   



\vspace{0.5cm}

# Programs Used & Version Information

- **GEMMA**: 0.97.3
- **R**: 3.6.1
- **Python**: 3
- **R/QTL2**: 0.28
- **R/QTL2FST**: 0.26
- **PLINK**: 1.90

\vspace{0.5cm}

# Generating Random Subsamples 
We performed 100 random subsamples in which we retained 500, 1,000, 1,500, 2,000, or 2,500 individuals (for fasting glucose we could not include 2,000 and 2,500 because the total sample size was smaller than 2,000) using the `sample()` function in `R`.  
\vspace{0.1cm}





# Genetic Relatedness Matrix (GRM)
We performed GWAS for each sub-sampled dataset using GEMMA. 
Since we used the LOCO (leave-one-chromosome-out) method when we do GWAS, we use a relatedness matrix (-k) that has been formed using all the SNPs EXCLUDING those on the chromosome currently being tested. We need 20 dosage files [bimbam file format], each missing a chromosome of SNPs, in order to make the matrix. As well as 20 individual dosage files for when you do the testing chromosome by chromosome. I found this to be very easily accomplished with a simple `grep` command. If you use `grep -v "${chrom}"` it will give you the original file, minus the chromosome you want to leave out. The dosage files and the GRMs can be downloaded from https://library.ucsd.edu/dc/object/bb9156620z  




_Script:_

```{r GRM, engine = 'bash', eval = FALSE,echo=T}

#Constructing GRM
gemma \
-g allExcept.${chrom}.P50_round2_3473_unpruned.bimbam \
-p pheno.txt \
-gk 1 \
-o allExcept.${chrom}.round2_unpruned_3473.dosages



```

_Arguments:_  
-g : genotype file.  
-p : phenotype file _Each line is a number indicating the phenotype value for each individual in turn, in the same order as in the genotype file. Missing phenotype information is denoted as NA. The number of rows should be equal to the number of individuals in the genotype file._.  
-gk 1 : Calculates the centered relatedness matrix.  
-o : Specifies the output file prefix.  

\vspace{0.1cm}


# GWAS for each subsampled dataset

_Script:_
```{r GWAS, engine = 'bash', eval = FALSE,echo=T}

gemma \
-g ${chrom}.round2_impute2_3473.bimbam \
-p pheno.txt \
-n 1 \
-k allExcept.${chrom}.round2_impute2_LDpruned0.95_3473.dosages.cXX.txt  \
-a ${chrom}.round2_128447.snpinfo \
-lmm 4 \
-o ${chrom}.output


```



_Arguments:_   
-g : genotype file.  
-p : phenotype file and should just be 1 columns per phenotype in the same sample order as the genotype file, NO HEADER.  
-k : relatedness matrix.  
-lmm 4 tells it to run all the tests and give you a beta estimate for the SNP effect.  
-a is a SNPINFO file that looks like this: The first column is SNP id, the second column is its base-pair position, and the third column is its chromosome number.  
-o is the output prefix.  

# Calling QTLs for each GWAS 
- An automated pipeline implemented in R was used to record the number of significant QTLs in each subsampled dataset.
- A fixed threshold of -log10P = 5.6, derived by permutation in Chitre et al. 2020, was used for all traits that were quantile normalized.
- Each chromosome was scanned to determine if there was at least one SNP that exceeded the threshold of –log10(p) > 5.6.
- To avoid situations where only a single, presumably anomalous, SNP showed a significant association, it was required that at least one other SNP within 0.5 Mb have a p-value that was within 2 –log10(p) of the index SNP.
- If a second supporting SNP was found, the identification of a QTL for that dataset was recorded.
- All SNPs with r2 > 0.4 relative to the just identified index SNP were excluded to avoid counting the same locus twice.
- The chromosome was then rescanned to see if any additional SNPs on this chromosome exceeded the threshold of –log10(p) > 5.6. If they did and were supported by a second SNP within 0.5 Mb that had a p-value that was within 2 –log10(p) of the index SNP, an additional QTL for that dataset was recorded.
- These steps were repeated as often as needed until no further significant QTLs could be identified on a given chromosome.
- This process was then continued for all subsequent chromosomes.
- After scanning the last chromosome, the number of QTLs detected for that dataset was tabulated.

_Script name:_ calling_qtls.R

# Linkage mapping with haplotypes 
- We performed linkage mapping with haplotypes using `R/qtl2`. Input files required for linkage mapping using haplotypes are present in this __Zenodo repository__.   
- We estimated founder haplotypes using the `calc_genoprob_fst` function with the cohort and founder strain genotypes.    
- The kinship matrices were derived using the “leave one chromosome out” method with the `calc_kinship` function.  

_Script name:_ geno_prob_and_kinship.R  

- For each sub-sampled dataset for the trait _BMI with tail_, we performed a genome scan using a linear mixed model with the `scan1` function.  

- We estimated the threshold using the `scan1perm` function.  
_Script name:_ permutations.R  

- We used the function `find_peaks` to identify LOD peaks that exceeded the permutation derived threshold of **18.2**.  
_Script name:_ scanone_sub_samples.R


