#' @importFrom ggplot2 ggplot aes xlab ylab theme_bw coord_flip ggsave geom_boxplot
#' @importFrom rlang .data
#' @importFrom reshape2 melt
#' @importFrom pheatmap pheatmap
#' @importFrom stats pairwise.wilcox.test
#' @importFrom utils write.table

#' @title DistBoxplot
#' @param dm A squared distance matrix or a 'dist' object.
#' @param group A factor indicates a metadata variable corresponds to samples in the distance matrix.
#' @param dm_name An optional character string for the distance name.
#' @param group_name An optional character string for the metadata variable name.
#' @param IndividualID Another metadata variable that can make this function return only paired distance values within this group.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' require(vegan)
#' dm<-vegdist(x, method="bray")
#' dm<-data.matrix(dm)
#' DistBoxplot(dm, y, dm_name="bray")
#' @author Shi Huang
#' @export
DistBoxplot<-function(dm, group, dm_name='Dist', group_name='', IndividualID=NULL, outdir=''){
  # Data check
  if("dist" %in% class(dm))
    dm <- data.matrix(dm)
  if(ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE){
    stop('The distance matrix is not squared or contains NA values!')
  }
  if(length(unique(group))==1){
    stop('At least two levels in your metadata variable are required.')}
  if(length(group)!=nrow(dm))
    stop('The number of rows in metadata and distance matrix are not equal')
  group<-factor(group)
  if(nlevels(group)>length(group)*0.95){
    stop('The number of levels in a certain category can not exceed 95% of total number of samples')
  }

  names(group)<-rownames(dm)
  if(!is.null(IndividualID)){
    names(IndividualID)<-rownames(dm)
  }

  #--------------------------------
  if(is.null(IndividualID)){
    colnames(dm)<-rownames(dm)<-paste(rownames(dm), group, sep="____")
    }else{
      colnames(dm)<-rownames(dm)<-paste(rownames(dm), group, IndividualID, sep="____") }

  dm[lower.tri(dm)]<-NA
  melt_dm<-reshape2::melt(data.matrix(dm))
  melt_dm<-melt_dm[!is.na(melt_dm$value),]
  melt_dm<-melt_dm[which(melt_dm[,1]!=melt_dm[,2]),]

  Row_Info <- data.frame(do.call(rbind,strsplit(as.character(melt_dm[,1]),"____",fixed=TRUE)))
  Col_Info <- data.frame(do.call(rbind,strsplit(as.character(melt_dm[,2]),"____",fixed=TRUE)))
  VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
  dm_value <- data.frame(VS, Row_Info, Col_Info, d=melt_dm$value)

  if(is.null(IndividualID)){
    colnames(dm_value) <- c("GroupPair", "Sample_1", "Group_1", "Sample_2", "Group_2", "Dist")
    DistType <- as.factor(dm_value$Group_1==dm_value$Group_2)
    DistType <- factor(DistType, levels=levels(DistType), labels=c("AllBetween", "AllWithin"))
    dm_value <- data.frame(DistType,dm_value)
  }else{
    colnames(dm_value) <- c("GroupPair", "Sample_1", "Group_1", "Subject_1", "Sample_2", "Group_2", "Subject_2", "Dist")
    Ind_dm_value <- dm_value[which(as.character(dm_value$Subject_1)==as.character(dm_value$Subject_2)), ]
  }
  #--------------------------------Output distance table and boxplot
  if(!is.null(IndividualID)){
    group_name<-paste(group_name, "IndividualID", sep="4")
    filepath <- sprintf('%s%s%s%s%s', outdir,dm_name, '.', group_name,'.Ind_dm_values.xls')
    sink(filepath); write.table(Ind_dm_value,quote=FALSE,sep='\t', row.names = FALSE); sink()

    #-----Output distance boxplot
    plot <- ggplot(Ind_dm_value, aes(x=.data$GroupPair, y=.data$Dist)) + geom_boxplot() +
                  xlab("Group pair") +
                  ylab(paste(dm_name, '_Distance', sep='')) + coord_flip() + theme_bw()
    suppressMessages(ggsave(filename=paste(outdir,dm_name,'.',group_name,'.boxplot.ggplot.pdf',sep=''), plot=plot, height=ifelse(nlevels(group)>2,nlevels(group),2)))

    invisible(Ind_dm_value)

  }else{
    #-----Objects to return
    if(length(dm_value$Dist)>4){
      p_t <- with(dm_value,t.test(Dist~DistType))$p.value
      p_w <- with(dm_value,wilcox.test(Dist~DistType))$p.value
      }else{
        p_t <- NA
        p_w <- NA}

    #-----Output distance-values table
    filepath<-sprintf('%s%s%s%s%s', outdir, dm_name, '.', group_name, '.dm_values.xls')
    sink(filepath); write.table(dm_value, quote=FALSE, sep='\t', row.names=FALSE); sink()

    #-----Output distance boxplot
    if(nlevels(group)<30){
      plot <- ggplot(dm_value, aes(x=.data$GroupPair, y=.data$Dist)) + geom_boxplot() +
        xlab("Group pair") +
        ylab(paste(dm_name, '_Distance', sep='')) + coord_flip() + theme_bw()
      suppressMessages(ggsave(filename=paste(outdir,dm_name,'.',group_name,'.boxplot.ggplot.pdf',sep=''), plot=plot))
      #-----Output statistical results
      p <- stats::pairwise.wilcox.test(dm_value$Dist,factor(dm_value$GroupPair))$p.value
      sink(paste(outdir, dm_name, '.', group_name, '.dm_stats_P_values--Wilcoxon-test.xls', sep=''))
      cat('\t'); write.table(p, quote=FALSE, sep='\t'); sink()
      if(length(unique(as.vector(p)[!is.na(p)]))>1)
        pheatmap::pheatmap(p, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers=T, main=group_name,
                           filename=paste(outdir,dm_name, '.', group_name, '.dm_stats_P_values--Wilcoxon-test.pdf', sep=''))
    }else{
      cat("Warning: neither a distance boxplot or p-value heatmap
          will output as the number of levels in the group has exceeded 30.")
    }

    result <- list()
    result$dm_value <- dm_value
    result$p_t <- p_t
    result$p_w <- p_w

    invisible(result)
  }
}
