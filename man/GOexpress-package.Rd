\name{GOexpress-package}

\alias{GOexpress-package}

\alias{GOexpress}

\docType{package}

\title{
Visualise microarray and RNAseq data with gene ontology annotations.
}

\description{
Integrates gene expression data with gene ontology annotations to score and
visualise genes and gene ontologies best clustering groups of experimental
samples. Supports custom annotations, or alternatively provides an interface
to the Ensembl annotations using the \code{biomaRt} package. The default
scoring approach is based on the random forest framework, while a one-way
ANOVA is available as an alteranative. GO term are scored and ranked
according to the average rank (alternatively, average power) of all
associated genes. P-values can be generated to assess the significance of
GO term ranking. The ranked list of GO terms is returned, with tools
allowing to visualise the statistics on a gene- and ontology-basis.
}

\details{
    \tabular{ll}{
    Package: \tab GOexpress\cr
    Type: \tab Package\cr
    Version: \tab 1.1.9\cr
    Date: \tab 2014-03-19\cr
    License: \tab  GPL (>= 3)\cr
    }
    This package requires only two input variables
    \enumerate{
    \item An \code{ExpressionSet} containing assayData and phenoData.
    The former
    should be a gene-by-sample matrix providing gene expression values for
    each gene in each sample. The latter should be an 
    \code{AnnotatedDataFrame} from the \code{Biobase} package providing
    phenotypic information and grouping factors
    with two or more levels.
    \item The name of the grouping factor to investigate, which must be a
    valid column name in the \code{phenoData} slot of the above
    \code{ExpressionSet}.
    }
    Following analysis, visualisation methods include:
    \itemize{
    \item Histogram and quantiles representations of the scores of GO terms
    \item Permutation-based P-values to assess the significance of GO term
    ranking
    \item Filtering of results on various criteria (e.g. number of genes
    annotated to GO term)
    \item Re-ordering of GO terms and gene result tables based on score or
    rank metric
    \item Table of statistics for genes annotated to a given GO term
    \item Hierarchical clustering of samples based on the expression level of
    genes annotated to a given GO term
    \item Heatmap of samples and genes based on the expression level of genes
    annotated to a given GO term
    \item Expression profile of a gene against one given factor (e.g. Time)
    while grouping samples on another given factor (e.g. Treatment)
    \item Univariate analysis of the expression level of a gene in the
    different groups of each experimental factor.
    \item Venn diagram of the counts of genes shared between a list of GO
    terms.
    }
}

\author{
    Maintainer: Kevin Rue-Albrecht <\code{kevin.rue@ucdconnect.ie}>
}

\examples{
# Sample input data available with package:
data(AlvMac)

# Sample output data available with package:
data(AlvMac_results)

# Supported species and microarrays:
data(microarray2dataset)
data(prefix2dataset)
}

\keyword{ package }
\keyword{ GOexpress }
\keyword{ randomForest }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }

\seealso{
    Main method for an example usage:
    \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
    
    Packages
    \code{\link[Biobase]{Biobase}},
    \code{\link[ggplot2]{ggplot2}},
    \code{\link[randomForest]{randomForest}},
    \code{\link[RColorBrewer]{RColorBrewer}},
    \code{\link[VennDiagram]{VennDiagram}}.
    
    Methods
    \code{\link[biomaRt:getBM]{getBM}},
    \code{\link[gplots:heatmap.2]{heatmap.2}},
    \code{\link[gplots:bluered]{bluered}},
    \code{\link[gplots:greenred]{greenred}},
    \code{\link[grid:grid.newpage]{grid.newpage}},
    \code{\link[grid:grid.layout]{grid.layout}},
    \code{\link[stringr:str_extract]{str_extract}}.
}
