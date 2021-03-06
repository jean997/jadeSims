% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot_roc_curves.R
\name{plot_roc_curves}
\alias{plot_roc_curves}
\title{Plot ROC curves for simulations}
\usage{
plot_roc_curves(agg.obj, which.stats, cols, ltys, lwd = 1.5,
  points.pch = 16, direction = "vertical", n.bars = 11, point.cex = 2.5,
  main = "", make.legend = FALSE)
}
\arguments{
\item{agg.obj}{An object produced by aggregate_sims}

\item{which.stats}{Names of statistics to plot}

\item{cols}{Colors}

\item{direction}{Which direction should the error bars go}

\item{main}{Plot title}

\item{make.legend}{Make a legend?}
}
\value{
Nothing
}
\description{
Plot ROC curves for JADE and any other methods
}

