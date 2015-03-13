\name{subEset}

\alias{subEset}

\title{
    Subset an ExpressionSet to values of phenotypic data columns.
}
\description{
    Given a list of column names and corresponding values present in the
    phenoData slot of an ExpressionSet object, this method returns the
    subset of the ExpressionSet restricted to samples associated with the
    given values in the given columns.
}

\usage{
    subEset(eSet, subset = list())
}

\arguments{
    \item{eSet}{
    \code{ExpressionSet} of the \code{Biobase} package including a
    gene-by-sample expression matrix in the \code{AssayData} slot, and a
    phenotypic information data-frame in the \code{phenodate} slot. In the
    expression matrix, row names are Ensembl gene identifiers or probeset
    identifiers, and column names are sample identifiers. In the phentypic
    data-frame, row names are sample idenfifiers, column names are grouping
    factors and phenotypic traits usable for the one-way ANOVA.
    }
    \item{subset}{
    A named list to subset \code{eSet}. Names must be column names existing
    in colnames(pData(eSet)). Values must be vectors of values existing in
    the corresponding column of pData(eSet).
    }
}

\value{
    A subset of the given ExpressionSet restricted to samples associated with
    the given values in the given columns
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Class \code{\link[Biobase:ExpressionSet]{ExpressionSet}}.
}

\examples{
# Load example data set
data(AlvMac)

# Subset it to only samples of "CN" and "MB" treatments, and also only "2H",
# "6H", and "24H" time-points
sub.AlvMac <- subEset(
    eSet=AlvMac,
    subset=list(
        Treatment=c("CN","MB"),
        Time=c("2H","6H")
        )
    )
}

\keyword{ GOexpress }
\keyword{ ExpressionSet }% __ONLY ONE__ keyword per line
