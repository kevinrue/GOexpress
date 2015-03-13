\name{AlvMac_results.pVal}

\alias{AlvMac_results.pVal}

\docType{data}

\title{
Sample output from the \code{pValue_GO()} function on an RNAseq experiment.
}

\description{
This variable may be used to test the filtering and visualisation methods
implemented in the package. It contains the output of 
successively applying the commands
\code{AlvMac_results = GO_analyse(eSet=AlvMac, f="Treatment")} and
\code{AlvMac_results.pVal = pValue_GO(AlvMac_results, N=100,
ranked.by=result$rank.by, rank.by='P')} to the toy input data \code{AlvMac}.
}

\usage{data("AlvMac_results.pVal")}


\value{
    A list of 9 slots summarising the input and results of the analysis:
    \itemize{
    \item{\code{GO} contains a table ranking all GO terms related to genes in
    the expression dataset based on the average ability of their related
    genes to cluster the samples according to the predefined grouping
    factor.}
    \item{\code{mapping} contains the table mapping genes present in the
    dataset to GO terms.}
    \item{\code{genes} contains a table ranking all genes present in the
    expression dataset based on their ability to cluster the samples
    according to the predefined grouping factor (see 'factor' below).}
    \item{\code{factor} contains the grouping factor analysed.}
    \item{\code{method} contains the statistical framework used.}
    \item{\code{subset} contains the  filters used to select a subset of
    samples from the original \code{ExpressionSet} for analysis.}
    \item{\code{rank.by} contains the metric used to rank the scoring tables.}
    \item{\code{ntree} contains number of trees built during the randomForest
    analysis.}
    \item{\code{mtry} contains the number of features randomly sampled as
    candidates at each split in each tree built during the randomForest
    analysis.}
    \item{\code{p.iterations} contains the number of permutations performed
    to compute the P-value in the \code{GO} slot.}
    }
}

\section{Warning}{
Running the above command again, you might obtain slightly different scores
and ranks due to the stochastic process of sampling used by the random forest
algorithm. However, the ranking metric was found to be robust and stable
across run, given adequate number of trees and predictor variables sampled.

To produce reproducible results, use the \code{set.seed()} function prior
to running any randomising or sampling function.
}


\source{
    Source data are part of a publication in review.
}

\examples{
data(AlvMac_results.pVal)
str(AlvMac_results.pVal)
head(AlvMac_results.pVal, n=20)
}

\keyword{datasets}
\keyword{GOexpress}
