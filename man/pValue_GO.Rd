\name{pValue_GO}
\alias{pValue_GO}

\title{
Compute p-values for ontologies by randomising gene labels.
}
\description{
Permutes the rank of genes processed by \code{GO_analyse()} to compute
a P-value for each ontology
}
\usage{
pValue_GO(result, N=1000, ranked.by=result$rank.by, rank.by='P')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.}
    \item{N}{
    The number of permutation desired. Default to 1000.
    }
    \item{ranked.by}{
    Either of 'rank' or 'score'. The metric used to compare whether GO terms
    in the randomised data arelative to their original value. Default to the
    current ordering method of the \code{result} object.
    }
    \item{rank.by}{
    Either of 'rank', 'score', or 'p.val'. The metric used to order the GO
    terms after computing of the P-value. Default to 'p.val'.
    }
}
\details{
This function is relatively lengthy. Its procedure can be divided in several
steps: (1) assemble a copy of all genes used in the original scoring step,
then for each of the N permutations, (2) permute the gene labels, (3)
aggregate the rank -- or score -- of all genes belonging to each GO term,
(4) decide whether each GO term is ranked -- or scored --better or worst than
originally, (5) return a p-value based on the number of permutations where
each GO term was better ranked -- or scored -- than originally.
}
\value{
    A list formatted identically to the results of the analysis, with an
    added 'p.val' column in the GO slot, and re-ordered by the chosen metric.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}
\examples{
\dontrun{
# Time-consuming examples:
# Limited here to N=1 for time constraint
# Recommended N=1000 (or at least 100).

# Run the analysis on factor "Treatment" including all samples
AlvMac.Pvals <- pValue_GO(result, N=1, ranked.by=result$rank.by, rank.by='P')
}
}

\keyword{ GOexpress }

