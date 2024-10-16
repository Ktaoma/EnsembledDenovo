library(msa)
library(Biostrings)
library(dplyr)
library(stringr)

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

args <- commandArgs(trailingOnly = TRUE)

print("# collect the result")
get_best <- function(output_name) {
  
  vec <- c()
  vec_DNA <- c()
  
  for (j in c(20,18,16,14,12,10)) {
    for (k in c(127,111,99,77,55,33,11)) {
      
      fl_n <- paste0(output_name,"/assembly/",j,"_",k,"/QC_AA_final.fasta")
      
      print(fl_n)
      
      seq_ <- readAAStringSet(fl_n) 
      fl_DNA <- paste0(output_name,"/assembly/",j,"_",k,"/QC_scaffold.fasta")
      seq_DNA <- readDNAStringSet(fl_DNA)
      
      
      
      print(names(seq_))
      if (length(seq_) == 0) {
        print("")
      }else{
        
        frame <- strsplit(strsplit(names(seq_),"_begin")[[1]][1],"frame=")[[1]][2] %>% as.numeric()
        names(seq_) <- paste0("Kmer_",j,"_Qual_",k,"_frame_",frame)
        names(seq_DNA) <- paste0("Kmer_",j,"_Qual_",k,"_frame_",frame)
        vec <- c(vec,seq_) 
        
        if (frame < 0) {
          vec_DNA <- c(vec_DNA,seq_DNA %>% reverseComplement())
        } else{
          vec_DNA <- c(vec_DNA,seq_DNA)
        }
        
      }
    }
  }

  AAvec <- do.call(c,vec)
  DNAvec <- do.call(c,vec_DNA)
  
  writeXStringSet(DNAvec,paste0(output_name,"/candidate_DNA.fa"))
  writeXStringSet(AAvec,paste0(output_name,"/candidate_AA.fa"))
  
  #######################################################################
  #AA check
  #######################################################################
  
  print("# filter AA length")  
  AAvec_ <- AAvec@ranges %>% as.data.frame() %>% 
    filter(width == 3391) %>% select(names) %>% 
    unlist() %>% as.vector()
  
  print("# check the similarity")
  ProtDF <- AAvec[AAvec_]
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
  
  #######################################################################
  #UTR check
  #######################################################################
  cand_DNAvec <- DNAvec[candidate]
  ref_den2 <- readDNAStringSet("ref_den2.fa")
  alnM <- c(cand_DNAvec,ref_den2)
  total_seq <- length(alnM)
  rep_seq <- total_seq-1
  
  tvec <- c()
  for (l in 1:rep_seq) {
    print(l)
    aln_pair <- alnM[c(l,total_seq)] %>% msa()
    alp_df <- aln_pair %>% as.matrix() %>% t() %>% as.data.frame()
    
    seq_01 <- table(alp_df[,1]) %>% as.data.frame() %>% filter(Var1 == "-")
    if (nrow(seq_01) == 0) {
      seq_01 <- 0
    } else{
      seq_01 <- seq_01$Freq
    }
    
    seq_02 <- table(alp_df[,2]) %>% as.data.frame() %>% filter(Var1 == "-") 
    if (nrow(seq_02) == 0) {
      seq_02 <- 0
    } else{
      seq_02 <- seq_02$Freq
    }
    
    Tot <- seq_01+seq_02
    
    tvec <- c(tvec,Tot)
    
  }
  
  cnt_dash <- tvec %>% as.data.frame()
  colnames(cnt_dash) <- "V1"
  cnt_dash$V2 <- rownames(cnt_dash)
  index_ <- cnt_dash %>% filter(V1 == min(V1))
  shortest_ <- index_$V2 %>% as.numeric()
  
  bestDNA <- alnM[shortest_]
  bestAA <- AAvec[names(bestDNA)]
  
  writeXStringSet(bestDNA,paste0(output_name,"/best_DNA.fa"))
  writeXStringSet(bestAA,paste0(output_name,"/best_AA.fa"))
  
  
  
}


get_best(args[1])


