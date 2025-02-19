#' A Reference Class to represent MS dataset.
#'
#' @field balance A length-one numeric vector.
#' @import limma
#' @import edgeR
#' @import RColorBrewer
#' @import psych
#' @import ggplot2
MSDataSet <- setRefClass('MSDataSet',
                         fields = list(
                           name='character',
                           analysis_type='character',
                           exclusion_status='character',
                           zero_substitution_status='character',
                           normalization_status='vector',
                           reference='logical',
                           meta='data.frame',
                           dataset='data.frame',
                           dataset_raw='data.frame',
                           dataset_annotation='data.frame',
                           dataset_annotation_raw='data.frame',
                           dataset_count='data.frame',
                           dataset_count_venn='data.frame',
                           dataset_count_excluded='data.frame',
                           samples='vector',
                           conditions='vector',
                           repeats='vector',
                           reference_samples='vector',
                           num_sample='numeric',
                           num_condition='numeric',
                           num_repeat='numeric',
                           mean_condition='data.frame',
                           mean_condition_scaled='data.frame',
                           sd_condition='data.frame'
                         ),
                         methods = list(
                           summary = function(){
                             cat('name: ',name,'\n')
                             cat('analysis_type: ',analysis_type,'\n')
                             cat('exclusion_status: ',exclusion_status,'\n')
                             cat('zero_substitution_status: ',zero_substitution_status,'\n')
                             cat('normalization_status: ',normalization_status,'\n')
                             cat('samples: ',samples,'\n')
                             cat('conditions: ',conditions,'\n')
                             cat('reference_samples: ',reference_samples,'\n')
                             cat('num_sample: ',num_sample,'\n')
                             cat('num_condition: ',num_condition,'\n')
                             cat('num_repeat: ',num_repeat,'\n')
                             cat('dataset: ',str(dataset),'\n')
                             cat('dataset annotation: ',str(dataset_annotation),'\n')
                           },
                           venn_repeat = function(zerotolinrep=0){
                             meta_noref <- meta[is.na(meta$reference),]
                             datacount <- dataset
                             datacount[datacount > 0] <- 1
                             i <- 1
                             for (r in repeats){
                               datacount[,as.character(r)] <- rowSums(datacount[,meta[meta$repeat. == r,'sample']])
                               datacount[(datacount[,as.character(r)] < (num_condition - zerotolinrep)),as.character(r)] <- 0
                               datacount[(datacount[,as.character(r)] >= (num_condition - zerotolinrep)),as.character(r)] <- i
                               i <- i * 10
                             }
                             datacount[,'venn'] <- rowSums(datacount[,as.character(repeats)])
                             dataset_count_venn <<- datacount
                           },
                           exclude = function(zerotolinrep=0,zerotolamongreps=0){
                             meta_noref <- meta[is.na(meta$reference),]
                             datacount <- dataset
                             datacount[datacount > 0] <- 1
                             for (r in repeats){
                               datacount[,as.character(r)] <- rowSums(datacount[,meta[meta$repeat. == r,'sample']])
                               datacount[(datacount[,as.character(r)] < (num_condition - zerotolinrep)),as.character(r)] <- 0
                               datacount[(datacount[,as.character(r)] >= (num_condition - zerotolinrep)),as.character(r)] <- 1
                             }
                             datacount[,'total'] <- rowSums(datacount[,as.character(repeats)])
                             datacount[(datacount[,'total'] < (num_repeat - zerotolamongreps)),'total'] <- 0
                             datacount[(datacount[,'total'] >= (num_repeat - zerotolamongreps)),'total'] <- 1
                             dataset <<- dataset[datacount[,'total'] == 1,]
                             dataset_annotation <<- dataset_annotation[datacount[,'total'] == 1,]
                             dataset_count <<- datacount
                             dataset_count_excluded <<- datacount[datacount[,'total'] == 1,]
                             exclusion_status <<- paste('Exluded by',zerotolinrep,'zero tol in rep and',zerotolamongreps,'zero tol among rep.')
                           },
                           normalization_sl = function(){
                             target <- mean(colSums(dataset))
                             norm_facs <- target / colSums(dataset)
                             prot_dat_sl <- sweep(dataset,2,norm_facs,FUN = '*')
                             dataset <<- prot_dat_sl
                             normalization_status <<- append(normalization_status,'SL')
                           },
                           normalization_tmm = function(){
                             #need edgeR
                             target <- calcNormFactors(dataset)
                             prot_dat_tmm <- sweep(dataset,2,target,FUN='/')
                             dataset <<- prot_dat_tmm
                             normalization_status <<- append(normalization_status,'TMM')
                           },
                           normalization_irs = function(){
                             i <- 1
                             irs <- cbind(as.numeric(rownames(dataset)),rowSums(dataset[,1:num_condition]))
                             while (i < num_repeat){
                               start_site <- i * num_condition + 1
                               stop_site <- i * num_condition + num_condition
                               irs <- cbind(irs, rowSums(dataset[,start_site:stop_site]))
                               i <- i + 1
                             }
                             irs <- irs[,-1]
                             irs$average <- apply(irs,1,function(x){exp(mean(log((x))))})

                             i <- 1
                             prot_dat_irs <- cbind(rownames(dataset),(dataset[,1:num_condition] * irs$average / rowSums(dataset[,1:num_condition])))
                             while (i < num_repeat){
                               start_site <- i * num_condition + 1
                               stop_site <- i * num_condition + num_condition
                               prot_dat_irs <- cbind(prot_dat_irs,(dataset[,start_site:stop_site] * irs$average / rowSums(dataset[,start_site:stop_site])))
                               i <- i + 1
                             }

                             rownames(prot_dat_irs) <- prot_dat_irs[,1]
                             prot_dat_irs <- prot_dat_irs[,-1]

                             dataset <<- prot_dat_irs
                             normalization_status <<- append(normalization_status,'IRS')
                           },
                           zero_substitution = function(least=30){
                             for (r in repeats) {
                               temp <- (dataset_count_excluded[,as.character(r)] == 1)
                               dataset[dataset[temp,] == 0] <<- least
                               zero_substitution_status <<- paste('zero values were substituted by',least)
                             }
                           },
                           graphic_draw = function(state=''){
                             #need RColorBrewer
                             qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
                             col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

                             par(mfrow = c(2, 2))
                             boxplot(log2(dataset),col = rep(col_vector[9:(num_repeat+8)],each=num_condition), main = state)
                             #need limma
                             plotDensities(log2(dataset),col = rep(col_vector[9:(num_repeat+8)],num_condition), main = state)
                             plotMDS(log2(dataset), col = col_vector[9:(num_condition+8)], main = state)
                             data_calCV_draw(dataset,num_condition,num_repeat, state)
                             par(mfrow = c(1, 1)) # reset to default
                           },
                           pearson_draw = function(state='',lm=TRUE){
                             prot_dat_reorder <- data_rearrange(dataset,num_condition,num_repeat)
                             i <- 1
                             #need psych
                             pairs.panels(log2(prot_dat_reorder[,1:num_repeat]), lm = lm,main=state)
                             while (i < num_condition) {
                               start_site <- i * num_repeat + 1
                               stop_site <- i * num_repeat + num_repeat
                               pairs.panels(log2(prot_dat_reorder[,start_site:stop_site]), lm = lm,main=state)
                               i <- i + 1
                             }
                           },
                           proteins_summary = function(){
                             temp <- dataset_raw[,meta$sample]
                             colnames(temp) <- meta$sample
                             temp_count <- temp
                             temp_count[temp_count > 0] <- 1
                             p_count <- cbind(colnames(temp_count),colSums(temp_count))
                             p_count <- data.frame(p_count)
                             p_count[,1] <- factor(p_count[,1],levels=p_count[,1])
                             p_count[,2] <- as.numeric(p_count[,2])
                             p_c <- ggplot(data=p_count,aes(x=p_count[,1],y=p_count[,2]))+geom_bar(stat='identity',aes(fill=p_count[,1])) +
                               theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust=1),legend.position = 'none',axis.text = element_text(size=12))+
                               xlab('')+ylab('')+ggtitle('Protein Count')

                             temp_inten <- temp
                             p_inten <- cbind(colnames(temp_inten),colSums(temp_inten))
                             p_inten <- data.frame(p_inten)
                             p_inten[,1] <- factor(p_inten[,1],levels=p_inten[,1])
                             p_inten[,2] <- as.numeric(p_inten[,2])
                             p_i <- ggplot(data=p_count,aes(x=p_inten[,1],y=p_inten[,2]))+geom_bar(stat='identity',aes(fill=p_inten[,1])) +
                               theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust=1),legend.position = 'none',axis.text = element_text(size=12))+
                               xlab('')+ylab('')+ggtitle('Protein Total Intensity')

                             return(list(p_c,p_i))
                           },
                           data_calMean_calSD = function(){
                             datameta <- meta
                             rownames(datameta) <- datameta$sample
                             condition_list <- unique(datameta[colnames(dataset),'condition'])

                             datamean <- as.numeric(rownames(dataset))
                             datasd <- as.numeric(rownames(dataset))
                             datamean <- data.frame(datamean)
                             datasd <- data.frame(datasd)
                             mean_coln <- c()
                             sd_coln <- c()
                             for (c in condition_list) {
                               sample_list <- datameta$sample[datameta$condition == c]
                               temp <- dataset[,sample_list]
                               mean_list <- apply(temp, 1, mean)
                               sd_list <- apply(temp, 1, sd)
                               mean_coln <- append(mean_coln,paste('mean_',c,sep=''))
                               sd_coln <- append(sd_coln,paste('sd_',c,sep=''))
                               datamean <- cbind(datamean,mean_list)
                               datasd <- cbind(datasd,sd_list)
                             }
                             rownames(datamean) <- datamean[,1]
                             datamean <- datamean[,-1]
                             rownames(datasd) <- datasd[,1]
                             datasd <- datasd[,-1]
                             colnames(datamean) <- mean_coln
                             colnames(datasd) <- sd_coln
                             mean_condition <<- datamean
                             sd_condition <<- datasd
                           },
                           data_mean_scaled = function(){
                             scaled_mean <- data_scale_mat(mean_condition)
                             colnames(scaled_mean) <- paste(colnames(scaled_mean),'_scaled',sep='')
                             mean_condition_scaled <<- scaled_mean
                           }
                         )
)

#' @title Read, clean and annotate MS TMT data.
#' @description Read, clean and annotate MS TMT data with Maxquant process data and meta file.
#' @details Input MaxQuang processed MS TMT data and meta file, then return a MSDataSet object.
#' @param prot_raw A DataFrame of proteinGroups from MaxQuant.
#' @param meta A DataFrame of metafile.
#' @param name A character of the name.
#' @param reference A boolean.
#' @param analysis_type A character.
#' @param exclusion_status A character.
#' @param normalization_status A vector.
#' @export
read_maxquant_prot <- function(prot_raw,meta,name='defalut',reference=FALSE,analysis_type='TMT',exclusion_status='none',normalization_status=c('none')){
  prot_dat_clean <- data_clean(prot_raw)
  prot_dat_extract <- data_extract(prot_dat_clean,meta)
  prot_annotation <- data_annotation(prot_dat_extract)

  msdataset <- MSDataSet$new(
    name=name,
    analysis_type=analysis_type,
    exclusion_status=exclusion_status,
    zero_substitution_status='none',
    normalization_status=normalization_status,
    reference=reference,
    meta=meta,
    dataset=prot_annotation[[1]],
    dataset_raw=prot_annotation[[1]],
    dataset_annotation=prot_annotation[[2]],
    dataset_annotation_raw=prot_annotation[[2]],
    dataset_count=data.frame(),
    dataset_count_venn=data.frame(),
    dataset_count_excluded=data.frame(),
    samples=meta$sample[is.na(meta$reference)],
    conditions=unique(meta$condition[is.na(meta$reference)]),
    repeats=unique(meta[,'repeat.']),
    reference_samples= ifelse(reference,meta$sample[meta$reference==1],'none'),
    num_sample=nrow(meta[is.na(meta$reference),]),
    num_condition=(nrow(meta[is.na(meta$reference),]) / max(meta[,'repeat.'])),
    num_repeat=length(unique(meta[,'repeat.'])),
    mean_condition=data.frame(),
    mean_condition_scaled=data.frame(),
    sd_condition=data.frame()
  )

  return(msdataset)
}

