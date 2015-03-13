\name{AlvMac_allGO}

\alias{AlvMac_allGO}

\docType{data}

\title{
Example of custom gene ontology annotations.
}

\description{
An example \code{data.frame} providing the identifier, name and namespace
corresponding
to gene ontology identifiers, compatible with the \code{AlvMac} example
\code{ExpressionSet}.
}

\usage{data("AlvMac_allGO")}

\value{
    A data frame detailing information about 13,302 gene ontologies:
    \itemize{
    \item{\code{go_id}}{ are gene ontology identifiers.}
    \item{\code{name_1006}}{ contains the corresponding gene ontology name.}
    \item{\code{namespace_1003}}{ contains the corresponding gene ontology
    namespace (i.e. "biological_process", "molecular_function", or
    "cellular_component")}
    }
}

\details{
This data-frame includes all gene ontologies present in the
\code{btaurus_gene_ensembl} dataset of the Ensembl BioMart server, including
those associated with no gene identifier in the \code{AlvMac} example
\code{ExpressionSet}.

See the help page of the \code{GO_analyse} function for an example usage.
}

\source{
These annotations were obtained from the Ensembl BioMart server using the
\code{biomaRt} package to access the Ensembl release 75 
\url{http://feb2014.archive.ensembl.org}.
}

\examples{
# Load the data
data(AlvMac_allGO)

# Structure of the data
str(AlvMac_allGO)

# First few rows
head(AlvMac_allGO)

}

\keyword{datasets}
\keyword{GOexpress}
