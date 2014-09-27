\name{quantiles_scores}

\alias{quantiles_scores}

\title{
Returns the quantiles of scores following an GOexpress analysis.
}

\description{
Returns a set of quantiles in indicating the scores reached by given
proportions of the GO terms.

This function can also be used on the output of \code{subset_scores()}
function as it returns a value formatted identically to the output of the
\code{GO_analyse()} function.
}

\usage{
quantiles_scores(result, probs=c(0.9, 0.95, 0.99, 0.999, 0.9999),
    quartiles=FALSE)
}

\arguments{
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.
    }
    \item{probs}{
    Numeric vector of probabilities with values in [0,1]. (Values up
    to 2e-14 outside that range are accepted and moved to the nearby
    endpoint.) See \code{quantile}
    }
    \item{quartiles}{
    A numeric vector of the percentiles for which the scores are desired. 
    }
}

\value{
    A named vector of percentiles and corresponding scores.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[stats:quantile]{quantile}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# filter for Biological Processes associated with 10+ genes
filtered_results <- subset_scores(result=AlvMac_results, total_count=5,
    namespace="BP")

# Quantiles of scores
quantiles_scores(result=filtered_results)
}

\keyword{ GOexpress }