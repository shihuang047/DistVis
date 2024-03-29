% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetadataChecker.R
\name{check_metadata}
\alias{check_metadata}
\title{check_metadata}
\usage{
check_metadata(metadata, more_missing_values = NULL)
}
\arguments{
\item{metadata}{A dataframe with > two columns corresponds to samples (rownames) in the biological data.}

\item{more_missing_values}{A optional string(s) can be added to define the missing values.}
}
\description{
check_metadata
}
\examples{
set.seed(123)
x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
a<-factor(c(rep("A", 29), NA, rep("B", 29), NA))
b<-factor(c(rep("A", 27), NA, "Not applicable", "Missing:not collected", rep("B", 28), NA, NA))
c<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
d<-c(sample(1:18), NA, "Not applicable", sample(5:44))
e<-c(rnorm(40), NA, "Not applicable", rnorm(18, 4))
e0<-c(NA, "Not applicable", rnorm(58, 4))
f<-rep("C", 60)
g<-rep(4, 60)
h<-c(rep("C", 59), "B")
i<-rep("Missing:not collected", 60)
metadata<-data.frame(a, b, c, d, e, e0, f, g, h, i)
metadata_summ<-check_metadata(metadata)
metadata_summ
}
\author{
Shi Huang
}
