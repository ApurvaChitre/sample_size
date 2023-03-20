library(data.table)



all_traits_chr<-list.files(path=".",full.names = F)


all_traits_chr<-gsub(".txt","",all_traits_chr)

column_names<-c('chr','snp','ps','beta','p_score')

#read gemma gwas assoc file
readdata<-function(x){
  
  data <- fread(paste0('cat *.',x,'.assoc.txt* | grep -v p_score | cut -f 1,2,3,8,15 '),header=F,stringsAsFactors =F)
  
  setnames(data, column_names)
  
}




all_data<-lapply(all_traits_chr,readdata)

#apply object names
names(all_data)<-all_traits_chr


for(i in seq_along(names(all_data))){
  all_data[[i]]$log10P <- -log10(all_data[[i]]$p_score)
  
  
}



####Following code will log which files have significant results

topSNP_part1<-data.frame()

for(i in seq_along(all_traits_chr)){
  
  
  if(length(which(all_data[[all_traits_chr[i]]]$log10P > 5.6  )) > 0){
    
    trait=all_traits_chr[i]
    above_threshold="yes"
    
    df <- data.frame(trait,above_threshold,stringsAsFactors = F)
    
  }else{
    trait=all_traits_chr[i]
    above_threshold="no"
    df <- data.frame(trait,above_threshold,stringsAsFactors = F)
  }
  
  topSNP_part1  <- rbind (topSNP_part1,df)
  
}

###find the chromosomes which contain significant results and store them in a list
signi_traits_n<-length(which(topSNP_part1$above_threshold=="yes"))
signi_traits<-topSNP_part1$trait[which(topSNP_part1$above_threshold=="yes")]

#this will contain chrs with signi results
signi_chrs <- setNames(replicate(signi_traits_n,data.frame()),signi_traits)


#5.6 is the genome wide significance threshold
for(i in seq_along(signi_traits)){
  
  chrs<-unlist(unique(all_data[[signi_traits[i]]][which(all_data[[signi_traits[i]]]$log10P > 5.6  ) ,"chr" ]))
  chrs<-unname(chrs)
  
  signi_chrs[[signi_traits[i]]]<-chrs
  
  
}



###subset assoc files to only signi chromosomes
t
master_list <- list()



sub_chr<-list()


for(i in seq_along(signi_traits)){
  chrs <-signi_chrs[[signi_traits[i]]]
  for(j in seq_along(chrs)){
    
    sub_chr[[chrs[j]]]<-all_data[[signi_traits[i]]][which(all_data[[signi_traits[i]]]$chr == chrs[j])]
    
  }
  master_list[[signi_traits[i]]]<-sub_chr
  
}

names(master_list)<-signi_traits


summary_QTL <- list()
topSNP_part2_list <- list ()



#define parameters that can be soft coded
#This is log10P - 1.5 (to see if there are any supporting SNPs)
#sup_P_term<-1.5
sup_P_term<-2
#This is 0.5 MB flanking dist- to see if there are any supporting snps for top snp or if it is a rogue snp
sup_dist<-500000
#LD command parameters for plink
ld_window_kb=11000
ld_r2=0.4
#QTL boundary distance 
qtl_dist<-1000000



for(i in seq_along(signi_traits)){
  topSNP_part2_list <- list ()  
  chrs <-signi_chrs[[signi_traits[i]]]
  for(j in seq_along(chrs)){
    
    #master_list[[signi_traits[1]]][1] 
    CHR <- all_data[[signi_traits[i]]][which(all_data[[signi_traits[i]]]$chr == chrs[j])]
    
    
    topSNP_part2<-data.frame()
    
    loop_run <- TRUE
    while(loop_run){
      #run until there are no topsnps
      #find top snp
      
      
      topsnp<-CHR$snp[which((CHR$log10P > 5.6) & (CHR$log10P==max(CHR$log10P)))]
      
      head(final_assoc)
      
      
      if(length(topsnp)>1){
        dups<-CHR[which(CHR$snp %in% topsnp),]
        dups<-dups[order(dups$ps),]
        topsnp<-dups$snp[which.max((abs(dups$beta)))]  
      }
      
      topsnp_log10P<-CHR$log10P[which(CHR$snp== topsnp)]
      sup_P<-(topsnp_log10P - sup_P_term)
      
      
      chr_topsnp<-gsub("chr","",strsplit(topsnp,split=":")[[1]][1])
      ps_topsnp<-strsplit(topsnp,split=":")[[1]][2]
      ps_topsnp<-as.numeric(ps_topsnp)
      
      start<- ps_topsnp -sup_dist
      stop<- ps_topsnp +sup_dist
      
      
      ##This evaluates if it is not a rogue SNP
      if(nrow(CHR[which((CHR$ps %in% seq(start,stop)) & (CHR$log10P > sup_P)),])>2){
        trait=signi_traits[i]  
        topsnp=topsnp
        QTL="yes"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)
        
      }else{
        trait=signi_traits[i]  
        topsnp=topsnp
        QTL="no"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)
      }
      
      topSNP_part2  <- rbind (topSNP_part2,df)
      
      system(paste0("plink --bfile P50_round2_3473_unpruned_sorted --chr ",chr_topsnp," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 1000000 -ld-window-kb ",ld_window_kb," --ld-window-r2 ",ld_r2," --out ","temp_qtl_n"),ignore.stdout = T,ignore.stderr = T,wait = T)
      
      

      #find next topsnp on the same chr which is not in the following list
      
      f<-fread("wc -l temp_qtl_n.ld")
      f<-f$V1
      
      
      
      if (f > 1){
        
        
        qtl_snps<-fread(paste0("cat temp_qtl_n.ld | grep -v R2 | awk '{ print $6 }'"),sep=" ",header = F,stringsAsFactors = F)
        
        #0.5 
        #1MB
        sub_start<-ps_topsnp-qtl_dist
        sub_stop<-ps_topsnp+qtl_dist
        
        nrow_next<-nrow(CHR[which((!CHR$snp %in% c(qtl_snps$V1,topsnp)) & (CHR$log10P>5.6) & (!CHR$ps %in% seq(sub_start,sub_stop))),])
        
        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$snp %in% c(qtl_snps$V1,topsnp)) & (!CHR$ps %in% seq(sub_start,sub_stop))),]
        }else{
          
          loop_run<-FALSE
        }
        
      }
      else{
        nrow_next<-nrow(CHR[which((!CHR$snp %in% topsnp) & (CHR$log10P>5.6) ),])
        
        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$snp %in% topsnp)),] 
        }else{
          loop_run<-FALSE 
        }
        
        
      }
      
      
    }
    
    topSNP_part2_list[[j]]<-topSNP_part2
    
  }
  
  summary_QTL[[i]] <- topSNP_part2_list
  
}





names(summary_QTL) <- signi_traits




temp <- unlist(summary_QTL, recursive = FALSE)
final <- do.call("rbind", temp)

write.csv(final, file.path(output_dir, paste0("QTL_vs_N_", traits_pattern, "_results.csv")), row.names = FALSE, quote = FALSE)