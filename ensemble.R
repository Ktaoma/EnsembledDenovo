library(msa)
library(Biostrings)
library(dplyr)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
print("# collect the result")
vec <- c()
for (j in c(20,18,16,14,12,10)) {

    for (k in c(127,111,99,77,55,33,11)) {
      
      fl_n <- paste0(args[1],"/assembly/",j,"_",k,"/QC_AA_final.fasta")
      print(fl_n)
      seq_ <- readAAStringSet(fl_n) 
      
      if (length(seq_) == 0) {
        print("")
      }else{
        names(seq_) <- paste0("Kmer_",j,"_Qual_",k)
        vec <- c(vec,seq_) 
      }
      
    }
}

print("# filter AA length")
cvec <- do.call(c,vec)
cvec_ <- cvec@ranges %>% as.data.frame() %>% 
        filter(width >= 3391) %>% select(names) %>% 
        unlist() %>% as.vector()

print("# check the similarity")
ProtDF <- cvec[cvec_]
AAcomb <- gtools::combinations(n=length(ProtDF), r=2, v=names(ProtDF))

score_c <- data.frame()
for (idx in c(1:nrow(AAcomb))) {
    dat_AA <- AAcomb[idx,]
    AA1 <- dat_AA[1]
    AA2 <- dat_AA[2]

    msa_AA <- c(ProtDF[AA1],ProtDF[AA2]) %>% msa(method="Muscle")
    msa_AA_df <- msa_AA %>% as.matrix() %>% t() %>% as.data.frame() 
    colnames(msa_AA_df) <- c("V1","V2")
    sim_score <- 100*(sum(ifelse(msa_AA_df$V1 == msa_AA_df$V2,1,0))/nrow(msa_AA_df))

    score <- c(AA1,AA2,sim_score) %>% 
            as.data.frame() %>% 
            t() %>% as.data.frame()

    score_c <- rbind(score_c,score)

}

print("# voting")
vote_ <- score_c %>% filter(V3 == 100)
candidate <- c(vote_$V1,vote_$V2) %>% unique()

print("# check depth")
dat <- data.frame()
for (x in candidate) {
    qual_ <- strsplit(x,"_")[[1]][2]
    kmer_ <- strsplit(x,"_")[[1]][4]
    cov_1 <- data.table::fread(paste0(args[1],"/assembly/",qual_,"_",kmer_,"/depth_virus.txt"))
    dat_ <- c(qual_,kmer_,mean(cov_1$V3)) %>% as.data.frame() %>% t() %>% as.data.frame() 
    dat <- rbind(dat,dat_)
}


print("# select best")
dat$V5 <- paste0("Quality_",dat$V1,"_Kmer_",dat$V2) 
dat %>% write.csv(paste0(args[1],"/assembly/final_candidate.csv"))
df <- data.table::fread(paste0(args[1],"/assembly/final_candidate.csv"))[,-1] %>% as.data.frame()
max_ <- df %>% as.data.frame() %>% group_by(V5) %>% slice_max(V5,n=1) %>% 
  as.data.frame() %>%
  select(V5) %>%
  unlist() %>% as.vector()
df$label <- ifelse(df$V5 %in% max_,df$V5,"")



df_sel <- df %>% filter(label != '') %>% group_by(V1) %>% top_n(1)
dat_final <- data.frame()
AA_vec <- c()
DNA_vec <- c()
for (z in df_sel$V5) {

  qual_ <- strsplit(z,"_")[[1]][2]
  kmer_ <- strsplit(z,"_")[[1]][4]
  
  
  #Amino Acid
  AA_ <- readAAStringSet(paste0(args[1],"/assembly/",qual_,"_",kmer_,"/QC_AA_final.fasta"))
  AA_vec <- c(AA_vec,AA_)

  
  #DNA 
  DNA_ <- readDNAStringSet(paste0(args[1],"/assembly/",qual_,"_",kmer_,"/QC_scaffold_final.fasta"))
  DNA_vec <- c(DNA_vec,DNA_)

  
  #Depth
  cov_1 <- data.table::fread(paste0(args[1],"/assembly/",qual_,"_",kmer_,"/depth_virus.txt"))
  dat_final <- rbind(dat_final,cov_1)
}


# write output
combineAA <- do.call(c,AA_vec)
combineDNA <- do.call(c,DNA_vec)
writeXStringSet(combineAA,paste0(args[1],"/assembly/AA_best.fa"))
writeXStringSet(combineDNA,paste0(args[1],"/assembly/DNA_best.fa"))
dat_final %>% write.csv(paste0(args[1],"/assembly/Depth_final.txt"))














