% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{search_target}
\alias{search_target}
\title{Retrieve Herbs, Molecules, and Targets Based on Gene List}
\usage{
search_target(gene_list)
}
\arguments{
\item{gene_list}{A character vector containing gene symbols which can be determined freely by users.}
}
\value{
A \code{data.frame} with three columns: \code{herb}, \code{molecule}, and \code{target},
containing information related to the specified genes. Rows with \code{NA} values are excluded.
}
\description{
This function retrieves herb, molecule, and target information from the internal dataset, \code{tcm_data},
based on a provided list of target genes.
}
\examples{
\dontrun{
# Example usage with a list of gene symbols
genes <- c("TP53", "EGFR", "BRCA1")
herbs_targets <- search_target(genes)
print(herbs_targets)
}

}
