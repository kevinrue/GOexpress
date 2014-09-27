\name{table_genes}

\alias{table_genes}

\title{
Returns a table listing the genes associated with a given Gene Ontology
}

\description{
Given a Gene Ontology (GO) identifier represented in the dataset and the
output variable of a \code{GO_analyse()} function, \code{table_genes()}
returns a table listing the genes associated with that go_id, their
associated gene name, and description.
}

\usage{
table_genes(go_id, result, data.only=FALSE)
}

\arguments{
    \item{go_id}{
    A Gene Ontology (GO) identifier.
    }
    \item{result}{
    The output of the \code{GO_analyse()} function.
    }
    \item{data.only}{
    Whether to return only the feature identifiers present in the given
    dataset or alternatively returns all feature identifiers associated with
    the GO term in the Ensembl BioMart.
    }
}

\value{
A data frame listing the statistics and annotations for the genes present in
the expression dataset and associated with the GO term.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# Table of result for genes associated with the GO term
# "toll-like receptor 4 signaling pathway"
table_genes(result=AlvMac_results, go_id="GO:0034142")
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ ontology }
