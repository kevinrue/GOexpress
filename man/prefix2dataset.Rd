\name{prefix2dataset}

\alias{prefix2dataset}

\docType{data}

\title{
Table mapping Ensembl gene identifier prefixes to BioMart datasets.
}

\description{
The species corresponding to an Ensembl gene identifier can typically be
identified from the prefix of the identifier (e.g. ENSBTAG corresponds to
\emph{Bos taurus}).
This table maps each known unique prefix to the corresponding species.
}

\usage{data(prefix2dataset)}

\value{
    A data frame with 66 rows.
    \describe{
        \item{\code{prefix}}{Unique prefix.}
        \item{\code{dataset}}{BioMart dataset.}
    }
}

\details{
\emph{C. elegans}, \emph{D. melanogaster}, and \emph{S. cerevisiae} have
atypical identifier pattern and prefixes. Prefixes for those are included in
the table, but checked independently of the other species.
}

\source{
A separate script was used to query the Ensembl BioMart server and build this
table.
}

\examples{
data(prefix2dataset)
str(prefix2dataset)
prefix2dataset
}

\keyword{datasets}
\keyword{GOexpress}
