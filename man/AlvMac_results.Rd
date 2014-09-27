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
    \describe{
    \item{ GO }{A table ranking all GO terms related to genes in the
    expression dataset based on the average ability of their related genes to
    cluster the samples according to the predefined grouping factor.}
    \item{ mapping }{The table mapping genes present in the dataset to GO
    terms.}
    \item{ genes }{A table ranking all genes present in the expression dataset
    based on their ability to cluster the samples according to the predefined
    grouping factor.}
    \item{ factor }{"Treatment", the predefined grouping factor.}
    \item{ method }{"randomForest", the statistical framework used.}
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
