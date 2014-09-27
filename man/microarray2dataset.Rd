\name{microarray2dataset}

\alias{microarray2dataset}

\docType{data}

\title{
Table mapping probeset identifier prefixes to datasets in the Ensembl BioMart.
}

\description{
The species corresponding to an probeset identifier can often be identified
from the prefix of the identifier (e.g. Bt.457.1.S1_at corresponds to
\emph{Bos taurus}). This table maps some known unique prefix to the
corresponding species.
}

\usage{data(microarray2dataset)}

\value{
    A data frame with 105 observations on the following 5 variables.
    \describe{
    \item{\code{dataset}}{BioMart species dataset.}
    \item{\code{microarray}}{Microarray identifier in the Ensembl BioMart
    dataset.}
    \item{\code{sample}}{An example probeset of the microarray.}
    \item{\code{prefix}}{The prefix pattern.}
    \item{\code{unique}}{Whether the prefix pattern is unique to the
    microarray.}
    }
}

\details{
All Agilent microarray share the same prefix pattern, making it very difficult
to differentiate. Many affymetrix microarrays also share the same prefix
pattern for several probesets.
}

\source{
A separate script was used to query the Ensembl BioMart server and build this
table.
}

\examples{
data(microarray2dataset)
microarray2dataset
}

\keyword{datasets}
\keyword{GOexpress}
