\name{AlvMac_GOgenes}

\alias{AlvMac_GOgenes}

\docType{data}

\title{
Example of custom mapping between gene ontology identifiers and gene features.
}

\description{
An example \code{data.frame} associating Ensembl gene identifiers to gene
ontology identifiers, compatible with the \code{AlvMac} example
\code{ExpressionSet}.
}

\usage{data("AlvMac_GOgenes")}

\value{
  A data frame with 188,237 associations between the following 2 variables.
  \itemize{
    \item{\code{gene_id}}{ are Ensembl gene identifiers of the
    \emph{Bos taurus} species.}
    \item{\code{go_id}}{ are gene ontology identifiers.}
  }
}

\details{
This data-frame includes all annotations between Ensembl gene identifiers
and gene ontology identifiers
present in the \code{btaurus_gene_ensembl} dataset of the Ensembl BioMart
server, including gene identifiers absent from the \code{AlvMac} example
\code{ExpressionSet}.

Importantly, gene identifiers present in these annotations but absent from
the \code{ExpressionSet} will be given a score of 0 (minimum valid
score; indicates no power to discriminates the predefined groups of
sample) and a rank equal to the number of genes present in the entire
dataset plus one (worst rank, while preserving discrete continuity of
the ranking). This is helpful where features considered uninformative
were filtered out of the \code{ExpressionSet}.

See the help page of the \code{GO_analyse} function for an example usage.
}

\source{
These annotations were obtained from the Ensembl BioMart server using the
\code{biomaRt} package on 10-Dec-2014.
}

\examples{
# Load the data
data(AlvMac_GOgenes)

# Structure of the data
str(AlvMac_GOgenes)

# First few rows
head(AlvMac_GOgenes)

}

\keyword{datasets}
\keyword{GOexpress}
