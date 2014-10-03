\name{AlvMac_results}

\alias{AlvMac_results}

\docType{data}

\title{
Sample output from the \code{GO_analyse()} function on an RNAseq experiment.
}

\description{
This variable may be used to test the filtering and visualisation methods
implemented in the package. It contains the output of the command
\code{AlvMac_results = GO_analyse(eSet=AlvMac, f="Treatment")}
applied to the sample input data \code{data(AlvMac)}.
}

\usage{data(AlvMac_results)}

\value{
    A list of 7 slots summarising the input and results of the analysis:
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
    \item{\code{ntree} contains number of trees built during the randomForest
    analysis.}
    \item{\code{mtry} contains the number of features randomly sampled as
    candidates at each split in each tree built during the randomForest
    analysis.}
    }
}

\section{Warning}{
Running the above command again, you might obtain slightly different scores
and ranks due to the stochastic process of sampling used by the random forest
algorithm. However, the ranking metric was found to be robust and stable
across run, given adequate number of trees and predictor variables sampled.
}

\source{
    Source data are part of a publication in review.
}

\examples{
data(AlvMac_results)
str(AlvMac_results)
head(AlvMac_results$GO, n=20)
head(AlvMac_results$GO$genes, n=20)
}

\keyword{datasets}
