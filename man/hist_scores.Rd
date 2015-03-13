\name{hist_scores}

\alias{hist_scores}

\title{
Plots the distribution of scores following an GOexpress analysis.
}

\description{
Plots the an histogram representing the frequencies of scores in the output
variable of the {GO_analyse()} function.

This function can also be used on the output of \code{subset_scores()}
function as it returns a value formatted identically to the output of the
\code{GO_analyse()} function.
}

\usage{
hist_scores(result,
            main=paste("Distribution of average scores in",
                        deparse(substitute(result))),
            xlab="Average score", ...)
}

\arguments{
    \item{result}{
    The output of the \code{GO_analyse()} function.
    }
    \item{main, xlab}{
    These arguments to title have useful defaults here.
    }
    \item{\dots}{
    Additional arguments passed on to \code{hist()}.
    }
}

\value{
    Returns the output of the \code{hist()} function.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[graphics:hist]{hist}},
    and \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}

\examples{
# load the sample output data with p.values computed
data(AlvMac_results.pVal)

# Histogram of scores (labelled with counts)
hist_scores(result=AlvMac_results, breaks=20, labels=TRUE)

# filter for Biological Processes associated with 5+ genes and <=0.05 P-value
filtered_results <- subset_scores(
    result=AlvMac_results.pVal, total_count=5, p.val=0.05,
    namespace="BP")

# Histogram of scores (labelled with counts)
hist_scores(result=filtered_results, breaks=20, labels=TRUE)
}

\keyword{ GOexpress }