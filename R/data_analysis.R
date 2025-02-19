#' @title Calculate mean and sd value by conditions
#' @description Calculate mean and sd value by conditions.
#' @details Input dataframe from previous processing, then return a data set with mean and sd values by conditions.
#' @param prot_dat A DataFrame from previous processing.
#' @param meta A DataFrame of meta file.
#' @return A DataFrame of mean and sd by conditions.
#' @export
data_calMean_calSD <- function(prot_dat, meta){
  rownames(meta) <- meta$sample
  condition_list <- unique(meta[colnames(prot_dat),'condition'])

  prot_mean_sd <- as.numeric(rownames(prot_dat))
  for (c in condition_list) {
    sample_list <- meta$sample[meta$condition == c]
    temp <- prot_dat[,sample_list]
    mean_list <- apply(temp, 1, mean)
    sd_list <- apply(temp, 1, sd)
    mean_sd <- cbind(mean_list,sd_list)
    mean_coln <- paste('mean_',c)
    sd_coln <- paste('sd_',c)
    colnames(mean_sd) <- c(mean_coln,sd_coln)
    prot_mean_sd <- cbind(prot_mean_sd,mean_sd)
  }
  rownames(prot_mean_sd) <- prot_mean_sd[,1]
  prot_mean_sd <- prot_mean_sd[,-1]

  return(prot_mean_sd)
}

#' @title Draw volcano plot
#' @description Draw volcano plot.
#' @details Input dataframe from previous processing, then return a ggplot object of volcano plot.
#' @param prot A MSDataSet object.
#' @param condition1 A character of condition 1.
#' @param condition2 A character of condition 2.
#' @param df A integer of freedom dgree.
#' @param xleft A numeric for x-axis limit on left.
#' @param xright A numeric for x-axis limit on right.
#' @param ydown A numeric for x-axis limit on bottom.
#' @param ytop A numeric for x-axis limit on top.
#' @param confidence A numeric of confidence interval.
#' @param s0 A numeric of fudge factor s0.
#' @param smooth_curve A boolean of significant curve or not.
#' @return A list of dataframe and ggplot object.
#' @export
#' @import ggplot2
#' @import ggrepel
data_volcano <- function(prot,condition1,condition2, df, xleft=-5, xright=5, ydown=0, ytop=5, confidence=.95, s0=.1,smooth_curve=TRUE){
  prot_dat_1 <- prot$dataset[,meta$sample[meta$condition == condition1]]
  rownames(prot_dat_1) <- prot$dataset_annotation$PROTEIN
  prot_dat_2 <- prot$dataset[,meta$sample[meta$condition == condition2]]
  rownames(prot_dat_2) <- prot$dataset_annotation$PROTEIN
  meta <- prot$meta
  ta <- qt(confidence,df)

  prot_dat_com <- cbind(prot_dat_1,prot_dat_2)
  colnum1 <- ncol(prot_dat_1)
  colnum2 <- ncol(prot_dat_2)
  pvalue <- apply(prot_dat_com, 1, function(x){t.test(x[1:colnum1],x[(colnum1 + 1):(colnum2+colnum1)],paired = F)$p.value})
  lgpvalue <- -log10(pvalue)

  prot_mean_sd <- data_calMean_calSD(cbind(prot_dat_1,prot_dat_2),meta)
  fc <- prot_mean_sd[,3] / prot_mean_sd[,1]
  l2fc <- log2(fc)

  prot_fudge <- data.frame(l2fc,lgpvalue)
  colnames(prot_fudge) <- c('l2fc','lgpvalue')
  rownames(prot_fudge) <- rownames(prot_dat_1)
  prot_fudge$fudge <- apply(prot_fudge,1, function(x){data_calsmoothcurve(as.numeric(x[1]),ta,s0,df)})
  prot_fudge$sig <- '0'
  prot_fudge$sig <- apply(prot_fudge,1,function(x){ifelse(((as.numeric(x[2]) > as.numeric(x[3])) & ((as.numeric(x[1])) > ta*s0)),'1',x[4])})
  prot_fudge$sig <- apply(prot_fudge,1,function(x){ifelse(((as.numeric(x[2]) > as.numeric(x[3])) & ((as.numeric(x[1])) < -ta*s0)),'-1',x[4])})

  prot_fudge_sig1 <- subset(prot_fudge,sig=='1')
  prot_fudge_sig2 <- subset(prot_fudge,sig=='-1')

  p <- ggplot(prot_fudge,aes(x=l2fc,y=lgpvalue,color=sig))+geom_point()+theme_bw()+
    scale_color_manual(values = c('red','grey60','green'))+theme(legend.position = 'none',title = element_text(size=20),axis.text = element_text(size=16))+
    xlab('Log2(Fold Change)')+ylab('-Log10(P-value)')+
    geom_text_repel(prot_fudge_sig2,mapping=aes(label=rownames(prot_fudge_sig2)),vjust=1,size=5,color='red')+
    geom_text_repel(prot_fudge_sig1,mapping=aes(label=rownames(prot_fudge_sig1)),vjust=1,size=5,color='green')+
    xlim(xleft,xright)+ylim(ydown,ytop)

  if (smooth_curve == TRUE){
    p <- p + geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(.02,xright),linetype='dashed',n=10000,size=2,alpha=.3,color='darkred')+
      geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(xleft,-.02),size=2,alpha=.3,color='darkgreen',n=10000,linetype='dashed')
  }

  return(list(prot_fudge,p))
}

#' @title Draw volcano plot with Labels
#' @description Draw volcano plot with manual labels.
#' @details Input dataframe from previous processing, then return a ggplot object of volcano plot.
#' @param prot A MSDataSet object.
#' @param condition1 A character of condition 1.
#' @param condition2 A character of condition 2.
#' @param df A integer of freedom dgree.
#' @param xleft A numeric for x-axis limit on left.
#' @param xright A numeric for x-axis limit on right.
#' @param ydown A numeric for x-axis limit on bottom.
#' @param ytop A numeric for x-axis limit on top.
#' @param confidence A numeric of confidence interval.
#' @param s0 A numeric of fudge factor s0.
#' @param smooth_curve A boolean of significant curve or not.
#' @param label A boolean of label or not.
#' @param label_list A list of labels.
#' @return A list of dataframe and ggplot object.
#' @export
#' @import ggplot2
#' @import ggrepel
data_volcano_label <- function(prot,condition1,condition2, df, xleft=-5, xright=5, ydown=0, ytop=5, confidence=.95, s0=.1,smooth_curve=TRUE,label=TRUE,label_list = 'none'){
  prot_dat_1 <- prot$dataset[,meta$sample[meta$condition == condition1]]
  rownames(prot_dat_1) <- prot$dataset_annotation$PROTEIN
  prot_dat_2 <- prot$dataset[,meta$sample[meta$condition == condition2]]
  rownames(prot_dat_2) <- prot$dataset_annotation$PROTEIN
  meta <- prot$meta
  ta <- qt(confidence,df)

  prot_dat_com <- cbind(prot_dat_1,prot_dat_2)
  colnum1 <- ncol(prot_dat_1)
  colnum2 <- ncol(prot_dat_2)
  pvalue <- apply(prot_dat_com, 1, function(x){t.test(x[1:colnum1],x[(colnum1 + 1):(colnum2+colnum1)],paired = F)$p.value})
  lgpvalue <- -log10(pvalue)

  prot_mean_sd <- data_calMean_calSD(cbind(prot_dat_1,prot_dat_2),meta)
  fc <- prot_mean_sd[,3] / prot_mean_sd[,1]
  l2fc <- log2(fc)

  prot_fudge <- data.frame(l2fc,lgpvalue)
  colnames(prot_fudge) <- c('l2fc','lgpvalue')
  rownames(prot_fudge) <- rownames(prot_dat_1)
  prot_fudge$fudge <- apply(prot_fudge,1, function(x){data_calsmoothcurve(as.numeric(x[1]),ta,s0,df)})
  prot_fudge$sig <- '0'
  prot_fudge$sig <- apply(prot_fudge,1,function(x){ifelse(((as.numeric(x[2]) > as.numeric(x[3])) & ((as.numeric(x[1])) > ta*s0)),'1',x[4])})
  prot_fudge$sig <- apply(prot_fudge,1,function(x){ifelse(((as.numeric(x[2]) > as.numeric(x[3])) & ((as.numeric(x[1])) < -ta*s0)),'-1',x[4])})

  prot_fudge_sig <- subset(prot_fudge,sig=='1')

  prot_fudge$label <- ifelse(rownames(prot_fudge) %in% label_list,'10','0')
  prot_fudge_label <- subset(prot_fudge,label=='10')

  prot_fudge$symbol <- as.character(as.numeric(prot_fudge$sig) + as.numeric(prot_fudge$label))

  if (label==TRUE){
    p <- ggplot(prot_fudge,aes(x=l2fc,y=lgpvalue,color=symbol))+geom_point(size=3,shape=1,stroke=1.5)+theme_bw()+
      scale_color_manual(values = c('red','grey60','green','blue','blue'))+
      theme(legend.position = 'none',title = element_text(size=20),axis.text = element_text(size=16))+
      xlab('Log2(Fold Change)')+ylab('-Log10(P-value)')+
      geom_text_repel(prot_fudge_label,mapping=aes(label=rownames(prot_fudge_label)),vjust=-0.5,hjust=-0.5,size=5,color='blue')+
      xlim(xleft,xright)+ylim(ydown,ytop)

    if (smooth_curve == TRUE){
      p <- p + geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(.02,xright),linetype='dashed',n=10000,size=2,alpha=.3,color='darkred')+
        geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(xleft,-.02),size=2,alpha=.3,color='darkgreen',n=10000,linetype='dashed')
    }

  }
  else {
    p <- ggplot(prot_fudge,aes(x=l2fc,y=lgpvalue,color=symbol))+geom_point(size=3,shape=1,stroke=1.5)+theme_bw()+
      scale_color_manual(values = c('red','grey60','green'))+theme(legend.position = 'none',title = element_text(size=20),axis.text = element_text(size=16))+
      xlab('Log10(Fold Change)')+ylab('-Log10(P-value)')+
      xlim(xleft,xright)+ylim(ydown,ytop)

    if (smooth_curve == TRUE){
      p <- p + geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(.02,xright),linetype='dashed',n=10000,size=2,alpha=.3,color='darkred')+
        geom_function(fun=function(x){data_calsmoothcurve(x,ta,s0,df)},xlim=c(xleft,-.02),size=2,alpha=.3,color='darkgreen',n=10000,linetype='dashed')
    }
  }

  return(list(prot_fudge,p))
}

data_calsmoothcurve = function(x,ta,s0,df){
  if (x > 0) {
    y = ta * (1 + (s0 / ((x/ta)-s0)))
    y = -log10(2*(1-pt(y,df=df)))
    return(y)
  }
  else if (x < 0) {
    y = ta * (1 + (s0 / ((x/-ta)-s0)))
    y = -log10(2*(1-pt(y,df=df)))
    return(y)
  }
}

#' @title Draw Specific Protein Intensity plot
#' @description Draw Specific Protein Intensity plot.
#' @details Input dataframe from data_transform_melt, then return a ggplot object of intensity plot.
#' @param prot_dat_specific A DataFrame from data_transform_melt.
#' @param id A integer or string.
#' @return A ggplot object.
#' @export
#' @import ggplot2
data_draw_specific <- function(prot_dat_specific, id){
  p <- ggplot(prot_dat_specific[prot_dat_specific[,1]==id,],aes(x=condition,y=value))+geom_point(size=2)+
    theme_bw()+theme(axis.title = element_text(size=20),axis.text = element_text(size=18,angle = 60,hjust = 1),title = element_text(size=22))+
    xlab('')+ylab('Relative Intensity')+ggtitle(id)
  return(p)
}

#' @title Transform DataFrame by melt
#' @description Transform DataFrame by melt.
#' @details Input dataframe from previous processing, then return a DataFame of conditions melt.
#' @param prot_dat A DataFrame from previous processing.
#' @param meta A DataFrame of meta file.
#' @return A DataFrame of melt.
#' @export
#' @import reshape2
data_transform_melt <- function(prot_dat, meta){
  prot_dat <- cbind(rownames(prot_dat),prot_dat)
  prot_dat_specific <- melt(prot_dat)
  prot_dat_specific$condition <- apply(prot_dat_specific,1,function(x){meta$condition[(meta$sample == x[2])]})
  return(prot_dat_specific)
}
