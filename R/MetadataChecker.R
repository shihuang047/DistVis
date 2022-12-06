#' @importFrom vegan vegdist


#' @title check_metadata
#' @param metadata A dataframe with > two columns corresponds to samples (rownames) in the biological data.
#' @param more_missing_values A optional string(s) can be added to define the missing values.
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
#' a<-factor(c(rep("A", 29), NA, rep("B", 29), NA))
#' b<-factor(c(rep("A", 27), NA, "Not applicable", "Missing:not collected", rep("B", 28), NA, NA))
#' c<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' d<-c(sample(1:18), NA, "Not applicable", sample(5:44))
#' e<-c(rnorm(40), NA, "Not applicable", rnorm(18, 4))
#' e0<-c(NA, "Not applicable", rnorm(58, 4))
#' f<-rep("C", 60)
#' g<-rep(4, 60)
#' h<-c(rep("C", 59), "B")
#' i<-rep("Missing:not collected", 60)
#' metadata<-data.frame(a, b, c, d, e, e0, f, g, h, i)
#' metadata_summ<-check_metadata(metadata)
#' metadata_summ
#' @author Shi Huang
#' @export
check_metadata <-function(metadata, more_missing_values=NULL){
    missing_values<-c("not applicable", "Not applicable", "Missing:not collected",
                      "Not provided", "missing: not provided",
                      "not provided", "not collected", "NA", NA, "")
    if(!is.null(more_missing_values)){
      missing_values<-c(miss_values, more_missing_values)
      missing_values
    }

    check_integers <- function(vector, all=TRUE){
      checker <- grepl("^[0-9]+$", as.character(vector), perl = T)
      total_len <- length(vector)
      n_integers <- sum(checker)
      n_non_integers <- total_len - n_integers
      #cat("Number of integers: ", n_integers, "\n")
      #cat("Number of non-integers: ", n_non_integers, "\n")
      if(all) if_all_integers <- all(checker) ## if any TRUE in the checker
      out <- list(if_all_integers=if_all_integers,
                  n_integers=n_integers,
                  n_non_integers=n_non_integers)
      return(out)
    }

    check_characters <- function(vector){
      checker <- grepl("[A-Za-z-$_]+", as.character(vector), perl = T) # [A-Za-z[:punct:]]+
      total_len <- length(vector)
      n_characters <- sum(checker)
      n_non_characters <- total_len - n_characters
      #cat("Number of characters: ", n_characters, "\n")
      #cat("Number of non-characters: ", n_non_characters, "\n")
      if_character_existed <- any(checker) ## if any TRUE in the checker
      out <- list(if_character_existed=if_character_existed,
                  n_characters=n_characters,
                  n_non_characters=n_non_characters)
      return(out)
    }

    check_numeric <- function(vector){
      checker <- suppressWarnings(as.numeric(as.character(vector)) %% 1 != 0)
      total_len <- length(vector)
      n_numeric <- length(checker[!is.na(checker)]) # if any NA values it means a string appear among the numeric values.
      n_non_numeric <- total_len - n_numeric
      #cat("Number of numeric: ", n_numeric, "\n")
      #cat("Number of non-numeric: ", n_non_numeric, "\n")
      if_numeric_existed <- any(!is.na(checker)) ## if any TRUE in the checker
      out <- list(if_numeric_existed=if_numeric_existed,
                  n_numeric=n_numeric,
                  n_non_numeric=n_non_numeric)
      return(out)
    }

    n_unique_values <- sapply(metadata, function(x) nlevels(factor(x)))
    n_missing_values <- sapply(metadata, function(x) sum(x %in% missing_values))
    n_real_values <- sapply(metadata, function(x) sum(!x %in% missing_values))
    #n_numeric_values <- suppressWarnings(sapply(metadata, function(x) sum(!is.na(as.numeric(x)))))
    #n_unique_real_values <- sapply(metadata, function(x) nlevels(as.factor(x[!x %in% missing_values])))
    if_character_existed <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[1]])
    n_characters <- sapply(metadata, function(x) check_characters(x[!x %in% missing_values])[[2]])
    if_all_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[1]])
    n_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[2]])
    if_numeric_existed <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[1]])
    n_numeric <- sapply(metadata, function(x) check_numeric(x[!x %in% missing_values])[[2]])
    all_values_identical <- apply(metadata, 2, function(x) length(unique(x))==1)
    numeric_var <- if_numeric_existed | (if_all_integers & n_unique_values >= 30)
    categorical_var <- (if_character_existed & !if_numeric_existed) | (if_all_integers & n_unique_values < 30)
    metadata_summ<-data.frame(metadata=colnames(metadata),
                              total_n=nrow(metadata),
                              n_missing_values,
                              n_real_values,
                              n_unique_values,
                              if_character_existed,
                              n_characters,
                              if_all_integers,
                              n_integers,
                              if_numeric_existed,
                              n_numeric,
                              numeric_var,
                              categorical_var,
                              all_values_identical)
    metadata_summ
    }

filter_samples_by_sample_ids_in_metadata <- function(data, dm=FALSE, metadata, ids_col=NA){
      if(dm & class(data)=="dist") data<-data.matrix(data)
      if(is.na(ids_col)){
        shared_ids<-intersect(rownames(data), rownames(metadata))
        metadata_idx<-which(rownames(metadata) %in% shared_ids)
      }else{
        shared_ids<-intersect(rownames(data), metadata[, ids_col])
        metadata_idx<-which(metadata[, ids_col] %in% shared_ids)
      }
      if(dm==TRUE){
        data_matched<-data[shared_ids, shared_ids]
        data_matched<-data_matched[order(rownames(data_matched)),order(colnames(data_matched))]
        cat("The number of samples in distance matrix (after filtering out samples with no metadata): ",
            nrow(data_matched) ,"\n")
      }else{
        data_matched<-data[shared_ids, ]
        data_matched<-data_matched[order(rownames(data_matched)),]
        cat("The number of samples in feature table (after filtering out samples with no metadata): ",
            nrow(data_matched) ,"\n")
      }

      metadata_matched<-metadata[metadata_idx, ]

      cat("The number of samples metadata (after filtering out samples with no metadata): ",
          nrow(metadata_matched) ,"\n")

      if(is.na(ids_col)){
        metadata_matched<-metadata_matched[order(rownames(metadata_matched)),]
        cat("The sample IDs are idenical in feature table and metadata: ",
            identical(rownames(data_matched), rownames(metadata_matched)), "\n")
      }else{
        metadata_matched<-metadata_matched[order(as.character(metadata_matched[, ids_col])),]
        cat("The sample IDs are idenical in feature table and metadata: ",
            identical(rownames(data_matched), as.character(metadata_matched[, ids_col])), "\n")
      }

      result<-list()
      result$data<-data_matched
      result$metadata<-metadata_matched
      return(result)
}

filter_dm_by_NA_in_target_field_of_metadata <- function(dm, metadata, target_field, ids_col=NA){
  if(is.na(ids_col)){
    SampleIDs<-rownames(metadata)
  }else{
    SampleIDs<-metadata[, ids_col]
  }
  if(!identical(rownames(dm), SampleIDs)) stop("The sample IDs should be idenical in feature table and metadata!")
  NAN_values<-c("not applicable", "Not applicable", "Missing:not collected",
                "Not provided", "missing: not provided",
                "not provided", "not collected", "NA", NA, "")
  idx<-which(!metadata[, target_field] %in% NAN_values)
  metadata_k<-metadata[idx, ]
  # check if the target field is numeric data
  temp <- as.numeric(metadata_k[, target_field])
  if(all(!is.na(temp))) metadata_k[, target_field] <- temp
  dm_k<-dm[idx, idx]
  cat("The number of kept samples (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$dm<-dm_k
  result$metadata<-metadata_k
  result
}

#' @title effect_size_eval
#' @param metadata A dataframe with > two columns corresponds to samples (rownames) in the biological data.
#' @param x A dataframe containing the biolgical data
#' @author Shi Huang
#' @example
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
#' a<-factor(c(rep("A", 29), NA, rep("B", 29), NA))
#' b<-factor(c(rep("A", 27), NA, "Not applicable", "Missing:not collected", rep("B", 28), NA, NA))
#' c<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' d<-c(sample(1:18), NA, "Not applicable", sample(5:44))
#' e<-c(rnorm(40), NA, "Not applicable", rnorm(18, 4))
#' e0<-c(NA, "Not applicable", rnorm(58, 4))
#' f<-rep("C", 60)
#' g<-rep(4, 60)
#' h<-c(rep("C", 59), "B")
#' i<-rep("Missing:not collected", 60)
#' metadata<-data.frame(a, b, c, d, e, e0, f, g, h, i)
#' metadata_summ<-check_metadata(metadata)
#' metadata_summ
#' out<-effect_size_eval(x, metadata)
#' out
#' @export
effect_size_eval <- function(x, metadata){
  dm<-vegan::vegdist(x)
  dm_list<-filter_samples_by_sample_ids_in_metadata(dm, dm=TRUE, metadata)
  metadata_summ<-check_metadata(metadata)
  all_group <- as.character(metadata_summ[which(!metadata_summ$all_values_identical), "metadata"])
  all_group_f <- as.character(metadata_summ[which(!metadata_summ$all_values_identical & metadata_summ$categorical_var), "metadata"])
  all_group_n <- as.character(metadata_summ[which(!metadata_summ$all_values_identical & metadata_summ$numeric_var), "metadata"])

  #--------------------------------
  # Statistical test: Adonis and Anosim
  #--------------------------------
  stat_summ <- data.frame(matrix(NA, nrow=length(all_group), ncol=9))
  rownames(stat_summ) <- all_group
  colnames(stat_summ) <- c("raw_sample_size", "filtered_sample_size", "num_class", "class_distribution",
                           "Adonis.F", "Adonis.R2", "Adonis.P","Anosim.R","Anosim.P")
  dm <- dm_list[[1]]
  metadata <- dm_list[[2]]
  #--------------------------------
  #suppressWarnings(
  for(group in all_group){
    stat_summ[group, 1] <- dim(dm)[1]
    # filter samples with NA in the metadata
    filtered_dm_list <- filter_dm_by_NA_in_target_field_of_metadata(dm, metadata, group)
    if(all(is.na(filtered_dm_list))){
      stat_summ[group, 2] <- 0
      stat_summ[group, 3:9] <- NA
    }else{
      dm_f <- filtered_dm_list$dm
      metadata_f <- filtered_dm_list$metadata
      y <- metadata_f[, group]
      stat_summ[group, 2] <- length(y)
      #--------------------------------
      if(nlevels(factor(y))==1){
        stat_summ[group, 5:9] <-NA #next
        cat("All values are identical in ", group,"!\n")
      }else{
        if(is.element(group, all_group_f)){
          y <- factor(y)
          stat_summ[group, 3] <- nlevels(y)
          stat_summ[group, 4] <- paste0(levels(y), collapse="|")
          if(all(table(y)!=1)){
            ano <- vegan::anosim(dm_f, y)
            stat_summ[group, 8] <- ano.R <- ano$statistic
            stat_summ[group, 9] <- ano.P <- ano$signif
            cat("ANOSIM (", group, "): \n")
            cat("--------------------------------")
            print(ano)
            ado <- vegan::adonis(dm_f ~ y)
            stat_summ[group, 5] <- ado.F <- ado$aov.tab$F.Model[1]
            stat_summ[group, 6] <- ado.F <- ado$aov.tab$R2[1]
            stat_summ[group, 7] <- ado.P <- ado$aov.tab$P[1]
            cat("ADONIS/PERMANOVA (",group,"): \n")
            cat("--------------------------------\n")
            print(ado$aov.tab)
            cat("--------------------------------\n\n")
          }else{
            stat_summ[group, 5:9] <- NA
          }
        }else{
          y<-as.numeric(as.character(y))
          stat_summ[group, 3]  <-  NA
          stat_summ[group, 4]  <-  NA
          ado <- vegan::adonis(dm_f ~ y)
          stat_summ[group, 5] <- ado.F <- ado$aov.tab$F.Model[1]
          stat_summ[group, 6] <- ado.F <- ado$aov.tab$R2[1]
          stat_summ[group, 7] <- ado.P <- ado$aov.tab$P[1]
          cat("ADONIS/PERMANOVA (",group,"): \n")
          cat("--------------------------------\n")
          print(ado$aov.tab)
          cat("--------------------------------\n\n")
        }
      }
    }
  }
  sink(paste(outpath, "Beta_diversity_summ.xls",sep=""));
  cat("\t");
  write.table(stat_summ, quote=FALSE,sep='\t',row.names=TRUE);
  sink()
  stat_summ
}

