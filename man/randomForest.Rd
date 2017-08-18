\name{randomForest}
\alias{randomForest}
\alias{randomForest.ExpressionSet}
\alias{randomForest.SummarizedExperiment}

\title{
Adaptation of \code{randomForest} for common \emph{Bioconductor} classes
}
\description{
Wrappers around the
\code{\link[randomForest:randomForest]{randomForest::randomForest}} method
for common \emph{Bioconductor} classes including
\code{\linkS4class{ExpressionSet}}
and \code{\linkS4class{RangedSummarizedExperiment}}
(and classes that inherits from them).

Notably, these wrappers require additional arguments to access the
appropriate data slots. Refer to the individual methods for details.
}

\usage{
\method{randomForest}{ExpressionSet}(
    x, pheno, assay = "exprs", ..., do.trace = 100
)
\method{randomForest}{SummarizedExperiment}(
    x, pheno, assay = "exprs", ..., do.trace = 100
)
}

\arguments{
    \item{x}{
An object that inherits from
\code{\linkS4class{ExpressionSet}}
or \code{\linkS4class{RangedSummarizedExperiment}}.
}
    \item{pheno}{
A valid phenotype column name in
\code{pData(x)}
or \code{colData(x)}.
}
    \item{assay}{
A valid assay name in
\code{names(assayData(x))} or
\code{names(assays(x))}.
}
    \item{\dots}{
Optional parameters to be passed to the original
\code{\link[randomForest:randomForest]{randomForest::randomForest}}.
}
    \item{do.trace}{
Running output is printed for every \code{do.trace} trees.
See \code{\link[randomForest:randomForest]{randomForest::randomForest}}
for details.
}
}

\value{
The direct output of the original
\code{\link[randomForest:randomForest]{randomForest::randomForest}}
method
(compatible with all original downstream methods; \emph{e.g.}, \code{plot}).
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    \code{\link[randomForest:randomForest]{randomForest::randomForest}},
    \code{\link[randomForest:importance]{importance}}.
}

\examples{
# Example data ----

nrows <- 100; ncols <- 10
featureIds <- sprintf("F\%03d", 1:100)
sampleIds <- sprintf("Sample\%02d", 1:10)
mat <- matrix(
    runif(nrows * ncols), nrow=nrows, ncol=ncols,
    dimnames = list(featureIds, sampleIds)
)
pdata <- data.frame(
    group = rep(LETTERS[1:2], each = ncols/2),
    row.names = sampleIds
)
rRanges <- GRanges(
    rep(c("chr1", "chr2"), c(75, 25)),
    IRanges(seq(1:100)*100-99, width=100),
    strand=sample(c("+", "-"), 100, TRUE),
    feature_id=featureIds)

eset <- ExpressionSet(
    assayData=mat,
    phenoData = AnnotatedDataFrame(pdata)
)

rse <- SummarizedExperiment(
    assays=SimpleList(exprs=mat),
    rowRanges=rRanges, colData=DataFrame(pdata)
)

# Available phenoData columns
names(pData(eset))

# Examples ----

## ExpressionSet
rf_eset <- randomForest(eset, "group")
## SummarizedExperiment
rf_rse <- randomForest(rse, "group")

plot(rf_eset)
plot(rf_rse)
}

\keyword{ GOexpress }
\keyword{ randomForest }
