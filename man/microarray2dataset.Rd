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
    A data frame with 105 observations of the 5 variables described below:
    \itemize{
    \item{\code{dataset} contains species-specific biomaRt dataset names.}
    \item{\code{microarray} contains microarray identifiers in the Ensembl
    BioMart dataset described above.}
    \item{\code{sample} contains an example probeset sampled from the
    corresponding Ensembl BioMart dataset and microarray.}
    \item{\code{prefix} a manually curated prefix pattern for the
    corresponding Ensembl BioMart dataset and microarray. This pattern may or
    may not be unique to the microarray in the BioMart database. Therefore
    we encourage users to use the 'microarray' .....}
    \item{\code{unique} contains a boolean value stating whether the prefix
    pattern was found curated as unique to the microarray or not.}
    }
}

\details{
All Agilent microarray share the same prefix pattern, making it very difficult
to differentiate. Many Affymetrix microarrays also share the same prefix
pattern for several probesets. 

The \code{dataset} and \code{microarray} arguments of the \code{GO_analyse}
method are the best way to specify which BioMart information should be
used to annotate the features in your expression dataset.
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
