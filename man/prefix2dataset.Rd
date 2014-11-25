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
    A data frame with 69 rows and 4 columns. Each row refers to one dataset
    in the Ensembl BioMart. The columns are described below:
    \itemize{
    \item{\code{dataset} contains the biomaRt dataset name.}
    \item{\code{species} contains the corresponding species name.}
    \item{\code{prefix} contains the corresponding unique prefix.}
    \item{\code{sample} contains a sample Ensembl gene identifier from this
    dataset.}    
    }
}

\details{
\emph{C. elegans}, \emph{D. melanogaster}, and \emph{S. cerevisiae} have
atypical identifier pattern and prefixes in their Ensembl gene identifiers.
However, the automatically
extracted prefix for \emph{C. elegans} and \emph{D. melanogaster} --
respectively "WBgene" and "FBgn" -- can be used as such to identify datasets
from those species. On the oher hand, prefixes used for the
\emph{S. cerevisiae} include "YHR", "YAL", and many others. Consequently,
expression data from \emph{S. cerevisiae} species is identified without
referring to the "prefix2dataset"" table; instead, such datasets are
identified if the first gene identifier in the dataset starts with "Y".
}

\source{
The \code{prefix2dataset.build} method stored in the \code{toolkit.R}
script was used to query the Ensembl BioMart server and build this table.
}

\examples{
data(prefix2dataset)
str(prefix2dataset)
prefix2dataset
}

\keyword{datasets}
\keyword{GOexpress}
