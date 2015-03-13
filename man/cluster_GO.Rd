\name{cluster_GO}

\alias{cluster_GO}

\title{
Generates a hierarchical clustering of the samples
}

\description{
Clusters the samples using only the expression levels of genes associated with
a given go_id.
}

\usage{
cluster_GO(
    go_id, result, eSet, f=result$factor, subset=NULL,
    method_dist="euclidean", method_hclust="average", cex=0.8,
    main=paste(go_id, result$GO[result$GO$go_id == go_id, "name_1006"]),
    xlab="Distance", cex.main=1, main.Lsplit=NULL, ...)

}

\arguments{
    \item{go_id}{
    A Gene Ontology (GO) identifier.
    }
    \item{result}{
    The output of \code{GO_analyse()} or a subset of it obtained from
    \code{subset_scores()}.
    }
    \item{eSet}{
    \code{ExpressionSet} of the \code{Biobase} package including a
    gene-by-sample expression matrix in the \code{AssayData} slot, and a
    phenotypic information data-frame in the \code{phenodata} slot. In the
    expression matrix, row names are Ensembl gene identifiers or probeset
    identifiers, and column names are sample identifiers. In the phentypic
    data-frame, row names are sample idenfifiers, column names are grouping
    factors and phenotypic traits usable for the one-way ANOVA. 
    }
    \item{f}{
    The grouping factor in \code{phenodata} to label the samples by.
    }
    \item{subset}{
    A named list to subset \code{eSet}. Names must be column names existing
    in colnames(pData(eSet)). Values must be vectors of values existing in
    the corresponding column of pData(eSet).
    }
    \item{method_dist}{
    The method used to calculate distance between samples. See the 
    \code{dist()} method from package \code{stats}.
    }
    \item{method_hclust}{
    The method used to cluster samples. See the 
    \code{hclust()} method from package \code{stats}.
    }
    \item{cex}{
    A numeric value defining the character expansion of text in the plot.
    }
    \item{main}{
    A character string for the main title of the plot.
    }
    \item{xlab}{
    A label for the x axis, defaults to "Distance".
    }
    \item{cex.main}{
    Scaling factor of the main title font size. Default is 1. We suggest to 
    use it in combination with the argument \code{main.Lsplit} for GO terms
    with long names.
    }
    \item{main.Lsplit}{
    Number of characters after which a new-line character will be inserted in
    the main title. If this would occur within a word, the new-line character
    will be inserted before this word. Default is NULL, leaving the title on a
    single line.
    }
    \item{\dots}{
    Additional parameters passed on to \code{dist()}, \code{hclust()}
    and \code{plot()}.
    }
}

\value{
    Returns the output of the \code{plot()} function.
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

# Hierarchical clustering of samples based on the same GO term
cluster_GO(
    go_id="GO:0034142", result=AlvMac_results, eSet=AlvMac, cex=0.7
    )

# Re-label sample by another factor
cluster_GO(
    go_id="GO:0034142", result=AlvMac_results, eSet=AlvMac, cex=0.7, f="Group"
    )

}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
