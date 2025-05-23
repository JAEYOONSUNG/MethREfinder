% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REBASE_sequence_generator.R
\name{filter_meth_by_mod_type}
\alias{filter_meth_by_mod_type}
\title{Filter methylation sensitivity data by modification type}
\usage{
filter_meth_by_mod_type(methylation_sensitivity, mod_type)
}
\arguments{
\item{methylation_sensitivity}{A data frame from REBASE containing enzyme information.}

\item{mod_type}{A character string specifying the modification type (e.g., "4mC").}
}
\value{
A filtered data frame containing only enzymes with the specified modification type.
}
\description{
Filters methylation sensitivity data based on the specified modification type.
}
\examples{
df_filtered <- filter_meth_by_mod_type(methylation_sensitivity, "4mC")

}
