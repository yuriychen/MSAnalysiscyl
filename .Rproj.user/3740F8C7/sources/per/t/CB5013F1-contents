data_clean <- function(prot_raw){
  temp <- (prot_raw$Potential.contaminant!='+') & (prot_raw$Reverse!='+') & (prot_raw$Only.identified.by.site!='+')
  prot_dat_clean <- prot_raw[temp,]
  return(prot_dat_clean)
}

data_annotation <- function(prot_dat_extract){
  prot_annotation <- prot_dat_extract[[2]]
  prot_annotation$ACCID <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'\\|'))[2]})
  prot_annotation$PROTEIN <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'\\|'))[3]})
  prot_annotation$PROTEIN <- apply(prot_annotation,1,function(x){unlist(strsplit(x[4],'_HUMAN'))[1]})
  prot_annotation$FULLNAME <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'_HUMAN'))[2]})
  prot_annotation$FULLNAME <- apply(prot_annotation,1,function(x){unlist(strsplit(x[5],'OS=Homo'))[1]})
  
  return(list(prot_dat_extract[[1]],prot_annotation))
}

data_extract <- function(prot_dat,meta){
  prot_dat_tmt <- prot_dat[,meta$channel]
  colnames(prot_dat_tmt) <- meta$sample
  
  return(list(prot_dat_tmt,prot_dat[,c('id','Fasta.headers')]))
}

data_calCV_draw <- function(prot_dat,condition_num,repeat_num,state){
  prot_dat_reorder <- data_rearrange(prot_dat,condition_num,repeat_num)
  
  i <- 1
  mean_list <- apply(prot_dat_reorder[,1 : repeat_num], 1, mean)
  sd_list <- apply(prot_dat_reorder[,1 : repeat_num], 1, sd)
  cv_list <- sd_list / mean_list * 100
  prot_dat_cv <- cbind(as.numeric(rownames(prot_dat)),cv_list)
  while (i < condition_num) {
    mean_list <- apply(prot_dat_reorder[,(i * repeat_num + 1) : (i * repeat_num + repeat_num)], 1, mean)
    sd_list <- apply(prot_dat_reorder[,(i * repeat_num + 1) : (i * repeat_num + repeat_num)], 1, sd)
    cv_list <- sd_list / mean_list * 100
    prot_dat_cv <- cbind(prot_dat_cv,cv_list)
    i <- i + 1
  }
  prot_dat_cv <- prot_dat_cv[,-1]
  cond_seq <- seq(1,condition_num,1)
  cond_seq_title <- paste('cond',cond_seq,sep = '')
  colnames(prot_dat_cv) <- cond_seq_title
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  boxplot(prot_dat_cv, col = col_vector[9:(condition_num+8)], main = state)
}

data_rearrange <- function(prot_dat,condition_num,repeat_num){
  i <- 1
  order_seq <- seq(1,(condition_num * repeat_num - condition_num + 1),condition_num)
  while(i < condition_num){
    i <- i + 1
    temp <- seq(i,(condition_num * repeat_num - condition_num + i),condition_num)
    order_seq <- c(order_seq,temp)
  }
  prot_dat_reorder <- prot_dat[,order_seq]
  return(prot_dat_reorder)
}

data_scale_mat <- function(prot){
  m <- apply(prot, 1, mean)
  s <- apply(prot, 1, sd)
  prot <- (prot - m) / s
  return(prot)
}

#' @title Zero Substitution in 3 Repeats.
#' @description Zero value substitution of 1 repeat in 3 repeats with average.
#' @details Input dataframe, and column start and end index of each repeat, then return a set of dataframe.
#' @param protdat A DataFrame from with only numeric value.
#' @param s1 A integer of start index of repeat 1.
#' @param e1 A integer of end index of repeat 1.
#' @param s2 A integer of start index of repeat 2.
#' @param e2 A integer of end index of repeat 2.
#' @param s3 A integer of start index of repeat 3 need substitution.
#' @param e3 A integer of end index of repeat 3 need substitution.
#' @param condition_num A integer of condition numbers.
#' @param least A integer for zero substitution for repeat 1 and 2.
#' @export
data_zero_substitution_3repeats <- function(protdat,s1,e1,s2,e2,s3,e3,condition_num,least=30){
  temp_1 <- protdat[,s1:e1]
  temp_1[temp_1 == 0] <- least
  temp_2 <- protdat[,s2:e2]
  temp_2[temp_2 == 0] <- least
  temp_3 <- protdat[,s3:e3]
  
  temp_12 <- cbind(temp_1,temp_2)
  temp_12_norm <- data_norm_irs(temp_12,condition_num,2)
  
  temp_1 <- temp_12_norm[,1:condition_num]
  temp_2 <- temp_12_norm[,(condition_num + 1):(condition_num * 2)]
  
  i <- 1
  while(i <= condition_num){
    temp_3[,i] <- (temp_1[,i] + temp_2[,i]) / 2
    i <- i + 1
  }
  
  protdat[,s1:e1] <- temp_1
  protdat[,s2:e2] <- temp_2
  protdat[,s3:e3] <- temp_3
  
  return(protdat)
}

data_norm_irs <- function(prot_dat,condition_num,repeat_num){
  i <- 1
  irs <- cbind(as.numeric(rownames(prot_dat)),rowSums(prot_dat[,1:condition_num]))
  while (i < repeat_num){
    start_site <- i * condition_num + 1
    stop_site <- i * condition_num + condition_num
    irs <- cbind(irs, rowSums(prot_dat[,start_site:stop_site]))
    i <- i + 1
  }
  irs <- irs[,-1]
  irs$average <- apply(irs,1,function(x){exp(mean(log((x))))})
  
  i <- 1
  prot_dat_irs <- cbind(rownames(prot_dat),(prot_dat[,1:condition_num] * irs$average / rowSums(prot_dat[,1:condition_num])))
  while (i < repeat_num){
    start_site <- i * condition_num + 1
    stop_site <- i * condition_num + condition_num
    prot_dat_irs <- cbind(prot_dat_irs,(prot_dat[,start_site:stop_site] * irs$average / rowSums(prot_dat[,start_site:stop_site])))
    i <- i + 1
  }
  
  rownames(prot_dat_irs) <- prot_dat_irs[,1]
  prot_dat_irs <- prot_dat_irs[,-1]
  
  return(prot_dat_irs)
}

#' @title Normalize MS TMT data with ERS.
#' @description Normalize MS TMT data with External Reference Signal.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with ERS.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @return A DataFrame.
#' @export
data_norm_ers <- function(prot_dat,condition_num,repeat_num,ref_col = c(1,7,13)){
  prot_ref <- prot_dat[,ref_col]
  temp <- prot_ref
  temp$rowsum <- rowSums(temp)
  prot_ref <- temp$rowsum / prot_ref
  
  i <- 1
  while(i <= ncol(prot_ref)){
    start_site <- (i - 1) * (condition_num + 1) + 1
    stop_site <- i * (condition_num + 1)
    prot_dat[,start_site : stop_site] <- prot_dat[,start_site : stop_site] * prot_ref[,i]
    i <- i + 1
  }
  
  return(prot_dat)
}