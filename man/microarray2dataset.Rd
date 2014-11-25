\name{microarray2dataset}

\alias{microarray2dataset}

\docType{data}

\title{
Table mapping probeset identifier prefixes to datasets in the Ensembl BioMart.
}

\description{
The species corresponding to an probeset identifier can often be identified
from the prefix of the identifier (e.g. Bt.457.1.S1_at corresponds to
\emph{Bos taurus} which is only associated with a single microarray platform
in the Ensembl BioMart). This table maps some manually curated unique patterns
to the corresponding species and microarray platform.
}

\usage{data(microarray2dataset)}

\value{
    A data frame with 158 rows and the 6 columns. Rows are sorted
    alphabetically by species name and each refer to a unique combination of
    dataset and microarray identifier in the Ensembl BioMart. The columns are
    described below:
    \itemize{
    \item{\code{dataset} contains species-specific biomaRt dataset names.}
    \item{\code{microarray} contains the microarray identifier in the Ensembl
    BioMart dataset above.}
    \item{\code{sample} contains a sample probeset from this combination of
    Ensembl BioMart dataset and microarray.}
    \item{\code{species} contains an example probeset sampled from the
    corresponding Ensembl BioMart dataset and microarray.}
    \item{\code{pattern} contains the corresponding manually curated pattern
    representative of probesets for this microarray platform. This pattern
    may or may not be unique to the microarray. Therefore
    we encourage users to use the 'microarray' argument of the
    \code{GO_analyse} method to specify the Ensembl BioMart dataset to use.}
    \item{\code{unique} contains a boolean value stating whether the
    pattern was found unique to the microarray or not.}
    }
}

\details{
All Agilent microarray share the same prefix pattern, making it very difficult
to differentiate. Many Affymetrix microarrays also share the same prefix
"AFFX" for several probesets. 

The \code{dataset} and \code{microarray} arguments of the \code{GO_analyse}
method are the best way to specify which BioMart information should be
used to annotate the features in your expression dataset.
}

\source{
The \code{microarray2dataset.build} method stored in the \code{toolkit.R}
script was used to query the Ensembl BioMart server and build this table.
}

\examples{
data(microarray2dataset)
microarray2dataset
}

\keyword{datasets}
\keyword{GOexpress}
