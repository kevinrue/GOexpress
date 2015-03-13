\name{AlvMac_allgenes}

\alias{AlvMac_allgenes}

\docType{data}

\title{
Example of custom gene feature annotations.
}

\description{
An example \code{data.frame} providing the identifier, name and description
corresponding to Ensembl gene identifiers present in the \code{AlvMac}
example \code{ExpressionSet}.
}

\usage{data("AlvMac_allgenes")}

\value{
    A data frame detailing information about the 100 gene features present in
    the \code{AlvMac} example \code{ExpressionSet}:
    \itemize{
    \item{\code{gene_id}}{ are Ensembl gene identifiers.}
    \item{\code{external_gene_name}}{ contains the corresponding gene name.}
    \item{\code{description}}{ contains the corresponding description.}
    }
}

\details{
This data-frame includes only Ensembl gene identifiers present in the
\code{AlvMac} example \code{ExpressionSet}.

See the help page of the \code{GO_analyse} function for an example usage.
}

\source{
These annotations were obtained from the Ensembl BioMart server using the
\code{biomaRt} package to access the Ensembl release 75 
\url{http://feb2014.archive.ensembl.org}.
}

\examples{
# Load the data
data(AlvMac_allgenes)

# Structure of the data
str(AlvMac_allgenes)

# First few rows
head(AlvMac_allgenes)
}

\keyword{datasets}
\keyword{GOexpress}
