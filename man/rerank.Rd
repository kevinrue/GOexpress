\name{rerank}

\alias{rerank}

\title{
Reorder the result variable by alternative metrics.
}

\description{
Reorder the ranked tables of GO terms and genes either by increasing
(average) rank or decreasing (average) score.
}

\usage{
rerank(result, rank.by = "rank")
}

\arguments{
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.
    }
    \item{rank.by}{
    Either of "rank" or "score"; the metric to rank the GO terms and genes by.
    }
}

\value{
    A list formatted identically to the results of the analysis, but ordered
    by the chosen metric.
}

\author{
    Kevin Rue-Albrecht
}

\note{
The name reorder() was not used to avoid conflict with package \code{stats}.
}

\seealso{
    Method \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# Re-rank the GO terms and genes based on the actual score instead of the rank
reranked.byScore <- rerank(result=AlvMac_results, rank.by="score")
}

\keyword{ GOexpress }
