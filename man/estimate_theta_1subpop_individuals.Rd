% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{estimate_theta_1subpop_individuals}
\alias{estimate_theta_1subpop_individuals}
\title{Estimate theta from individuals}
\usage{
estimate_theta_1subpop_individuals(individuals)
}
\arguments{
\item{individuals}{Individuals to get haplotypes for.}
}
\value{
List:
\itemize{
\item estimate: Vector of length 1 containing estimate of theta or NA if it could not be estimated
\item error: true if an error happened, false otherwise
\item details: contains description if an error happened
}
}
\description{
Estimate theta for one subpopulation given a sample of genotypes.
}
