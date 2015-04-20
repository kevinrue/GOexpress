\name{heatmap_GO}

\alias{heatmap_GO}

\title{
Generates a heatmap and hierarchical clustering of the samples and the genes
}

\description{
Clusters the samples and the genes associated with a GO term using the
expression levels of genes related to a given ontology. Represents expression
levels of those genes in a heatmap.
}

\usage{
heatmap_GO(
    go_id, result, eSet, f=result$factor, subset=NULL, gene_names=TRUE,
    NA.names=FALSE, margins=c(7 ,5),
    scale="none", cexCol=1.2, cexRow=0.5, 
    labRow=NULL,
    cex.main=1, trace="none", expr.col=bluered(75), 
    row.col.palette="Accent",
    row.col=c(),
    main=paste(
        go_id, result$GO[result$GO$go_id == go_id,"name_1006"]
        ),
    main.Lsplit=NULL,
    ...)
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
    gene-by-sample expression matrix in the \code{assayData} slot, and a
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
    \item{gene_names}{
    A boolean value. Default is TRUE, to label genes by their associated gene
    name. If FALSE, labels the genes by their feature identifier in the
    expression dataset (i.e. Ensembl gene identifier or microarray probeset).
    }
    \item{NA.names}{
    A boolean value. Default is TRUE. If labelling genes by their associated
    gene name (see argument \code{gene_names}), whether to display the
    gene feature identifier for gene features without associated gene name.
    }
    \item{margins}{
    A numeric vector of length 2 specifying the number of lines of margins to
    apply for the bottom and right margins, respectively.
    Defaults are 7 (bottom) and 5 (right). For Ensembl gene identifiers, we
    suggest setting manually a bottom margin of 13. See \code{heatmap.2()}.
    }
    \item{scale}{
    Character indicating if the values should be centered and scaled in either
    the row direction or the column direction, or none. Default is "none".
    See \code{heatmap.2()}.
    }
    \item{cexCol, cexRow}{
    Positive numbers, used as cex.axis in for the row or column axis labeling.
    Defaults are 1.2 and 1, respectively. See \code{heatmap.2()}.
    }
    \item{labRow}{
    A character vector of names to re-label the rows (samples).
    If vector of length 1, it is assumed to be a column name in the
    \code{phenoData} slot. Default are the values of the factor \code{f}.
    See \code{heatmap.2()}. 
    }
    \item{cex.main}{
    Scaling factor of the main title font size. Default is 1. We suggest to 
    use it in combination with the argument \code{main.Lsplit} for GO terms
    with long names.
    }
    \item{trace}{
    Character string indicating whether a solid "trace" line should be drawn
    across each 'row' or down each 'column', both' or 'none'. The distance of
    the line from the center of each color-cell is proportional to the size of
    the measurement. Defaults to 'none'.
    }
    \item{expr.col}{
    Character vector indicating the colors to represent the different levels
    of gene expression. Defaults to a colormap of 75 shades ranging from blue
    (low) to red (high) and centered around white. If using differential
    expression data, you should probably use \code{greenred(75)} instead.
    }
    \item{row.col.palette}{
    A valid \code{RColorBrewer} palette name to fetch the colormap from, to
    color-code the groups of samples.
    }
    \item{row.col}{
    A vector of color names or codes. The number of colors provided must match
    the number of levels of the grouping factor. Default to an palette of up
    to 9 colors marking the different levels of the predefined grouping factor
    on the left side of the heatmap.
    }
    \item{main}{
    Main title of the figure. Default is paste(go_id, go_name).
    }
    \item{main.Lsplit}{
    Number of characters after which a new-line character will be inserted in
    the main title. If this would occur within a word, the new-line character
    will be inserted before this word. Default is NULL, leaving the title on a
    single line.
    }
    \item{\dots}{
    Additional arguments passed on to \code{heatmap.2()}.
    }
}

\value{
    Returns the output of the \code{heatmap.2()} function.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[gplots:heatmap.2]{heatmap.2}},
    \code{\link[GOexpress:GO_analyse]{GO_analyse}},
    and \code{\link[RColorBrewer:brewer.pal.info]{brewer.pal.info}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# Heatmap the top-ranked GO term (toll-like receptor 4 signaling pathway) as
# example
heatmap_GO(go_id="GO:0034142", result=AlvMac_results, eSet=AlvMac)

# Same with larger sample labels on the right hand side.
heatmap_GO(go_id="GO:0034142", result=AlvMac_results, eSet=AlvMac, cexRow=1)

# Change the color-coding to green-black-red gradient (more appropriate for
# differential expression values)
library(gplots)
heatmap_GO(
    go_id="GO:0034142", result=AlvMac_results, eSet=AlvMac,
    expr.col=greenred(75)
    )
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }