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
rerank(result, rank.by = 'rank')
}

\arguments{
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.
    }
    \item{rank.by}{
    Either of 'rank', 'score' or 'p.val'; the metric to rank the GO terms and
    genes. Note that 'pval' is only applicable on the output of the
    \code{pValue_GO()} function. See details for breaking ties.
    }
}

\value{
    A list formatted identically to the results of the analysis, but ordered
    by the chosen metric.
}

\author{
    Kevin Rue-Albrecht
}

\details{
Taking an example, to rank GO terms by P-value and break ties by average
rank, rerank first by 'rank', and then rerank the resulting object by 'p.val'.
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

# load the sample output data with p.values computed
data(AlvMac_results.pVal)

# To rank by P-value, while breaking the ties by ave_rank,
# rank first by ave_rank
reranked.byRank <- rerank(result=AlvMac_results.pVal, rank.by="rank")
# rank the result by P-value
reranked.pVal_rank <- rerank(result=reranked.byRank, rank.by="p.val")

}

\keyword{ GOexpress }
