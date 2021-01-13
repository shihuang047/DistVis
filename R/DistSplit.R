#' @importFrom vegan vegdist

#' @title DistBoxplot
#' @param dm A squared distance matrix or a 'dist' object.
#' @param metadata A dataframe with > two columns corresponds to samples (rownames) in the distance matrix.
#' @param split_factor A metadata variable corresponds to samples in the distance matrix and used for split the distance matrix.
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
#' z<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' metadata<-data.frame(y, z)
#' dm<-dist(x)
#' split_distance_matrix_by_metadata(dm, metadata, split_factor="y")
#' @author Shi Huang
#' @export
split_distance_matrix_by_metadata<-function(dm, metadata, split_factor){
  #if(is.null(rownames(dm))#
  if(!is.element(split_factor, colnames(metadata)))
    stop("The split_factor you specified should be one of the column names of input metadata.")
  if(class(dm)=="dist") dm <- data.matrix(dm)
  f<-metadata[, split_factor]
  sub_dm_name_list<-split(1:nrow(dm), f, drop=TRUE)
  sub_dm_list<-
    lapply(sub_dm_name_list, function(x){
    list(dm[x, x], metadata[x, ])
  })

  sub_dm_list

}
