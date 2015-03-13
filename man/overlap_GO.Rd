\name{overlap_GO}

\alias{overlap_GO}

\title{
Shared genes between a list of GO terms.
}

\description{
Given a list of two to five GO terms, overlap_GO() will produce a Venn diagram
showing the counts of overlapping genes associated with those GO terms.
}

\usage{
overlap_GO(go_ids, result, filename=NULL, mar=rep(0.1, 4), ...)
}

\arguments{
    \item{go_ids}{
    A character vector of GO term identifiers to compare. For instance,
    \code{head(result$scores$go_id, n=5)} to compare the first 5 top-
    ranked GO terms in the\code{result} variable.
    }
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.
    }
    \item{filename}{
    The output filename where the Venn diagram will be saved. Default is NULL,
    displaying the Venn diagram on screen.
    }
    \item{mar}{
    The margins around the Venn diagram. Some Venn diagrams place the GO term
    labels outside the visible frame.
    }
    \item{\dots}{
    Further parameters forwarded to \code{the venn.diagram()} function.
    }
}

\value{
    Returns the output of the \code{venn.diagram()} function.
}

\author{
    Kevin Rue-Albrecht
}

\section{Warning }{
    An error is returned if the list of GO term identifiers contains less than
    2 elements or more than 5, as the underlying \code{venn.diagram()} method
    does not support values outside that range.
}

\seealso{
    Method \code{\link[VennDiagram:venn.diagram]{venn.diagram}}.
}

\examples{
# load the sample output data with p.values computed
data(AlvMac_results.pVal)

# filter for Biological Processes associated with 5+ genes and <=0.05 P-value
filtered_results <- subset_scores(
    result=AlvMac_results.pVal, total_count=5, p.val=0.05,
    namespace="BP")

# Venn diagram of overlapping genes between top 5 GO terms
overlap_GO(
    go_ids=head(filtered_results$GO$go_id, n=5),
    result=filtered_results, filename="VennDiagram.tiff")
}

\keyword{ GOexpress }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ Venn }
\keyword{ ontology }
