\name{expression_plot_symbol}

\alias{expression_plot_symbol}

\title{
Plots the expression profile of a gene by levels of a factor
}

\description{
This function will plot the expression profile of a gene across a valid
X-axis variable from the \code{AnnotatedDataFrame} while representing the mean
and confidence interval of groups of samples defined by levels of a valid
grouping factor from the AnnotatedDataFrame.

In the case of a gene name represented by multiple gene or probeset
identifiers in the dataset, a lattice of plots will be produced. Each of the
plots in this lattice can subsequently be plotted separately using its
associated index.
}

\usage{
expression_plot_symbol(
    gene_symbol, result, eSet, x_var, f=result$factor, subset=NULL,
    index=0, xlab=x_var, ylab="log2(cpm)", ylim=range(exprs(eSet)), 
    col.palette="Accent",
    col=brewer.pal(n=length(levels(pData(eSet)[,f])), name=col.palette),
    level=0.95, titles=c(), title.size=2, axis.title.size=20,
    axis.text.size=15, axis.text.angle=0,
    legend.title.size=20, legend.text.size=20, legend.key.size=30)
}

\arguments{
    \item{gene_symbol}{
    A gene name present in \code{rownames(AlvMac_results$genes)}.
    }
    \item{result}{
    An output of the \code{GO_analyse()} or \code{subset_scores()} function.
    }
    \item{eSet}{
    \code{ExpressionSet} of the \code{Biobase} package including a
    gene-by-sample expression matrix in the \code{AssayData} slot, and a
    phenotypic information data-frame in the \code{phenodate} slot. In the
    expression matrix, row names are Ensembl gene identifiers or probeset
    identifiers, and column names are sample identifiers. In the phentypic
    data-frame, row names are sample idenfifiers, column names are grouping
    factors and phenotypic traits usable for the one-way ANOVA. 
    }
    \item{x_var}{
    A column name in \code{phenodata} to plot on the X-axis. If representing
    time on the X-axis, users should store the time-points as numeric values
    in the \code{AnnotatedDataFrame} to adequately space the time-points.
    }
    \item{f}{
    A column name in \code{phenodata} to group the samples when representing
    mean and confidence interval. The factor specified in the initial
    \code{GO_analyse()} call is used by default. Unexpected grouping factors
    of samples can reveal interesting trends (e.g. "Animal", "Tissue",
    "CellType", ...).
    }
    \item{subset}{
    A named list to subset \code{eSet}. Names must be column names existing
    in colnames(pData(eSet)). Values must be vectors of values existing in
    the corresponding column of pData(eSet).
    }
    \item{index}{
    In the case where multiple gene or probeset identifiers are associated
    with the gene name, \code{index=2} will plot the expression profile of the
    second feature identifier alone, for instance.
    }
    \item{xlab}{
    Title of the X-axis. Default is tha value of \code{x_var}.
    }
    \item{ylab}{
    Title of the Y-axis. Default is "log2(cpm)".
    }
    \item{ylim}{
    Numeric vector of length 2 specifying the lower and upper bounds of the
    Y axis. Default is scaled to the full range of expression values in the 
    expression dataset, to ease comparison of different genes. If set to
    NULL, the axis will be scaled to fit the plotted data only.
    }
    \item{col.palette}{
    A valid \code{RColorBrewer} palette name to fetch the colormap from.
    }
    \item{col}{
    A vector of color names or codes. The number of colors provided must match
    the number of levels of the grouping factor. Default to a palette with an
    adequate set of colors. If specified, overrides argument
    \code{col.palette}.
    }
    \item{level}{
    The confidence interval level to visualise around the mean of each group.
    Default is 0.95.
    }
    \item{titles}{
    Character vector providing as many titles as there are plots to replace
    the default titles. Default is a combination of the gene id and the
    associated gene.
    }
    \item{title.size}{
    Changes the font size of the title. Default is 2.
    }
    \item{axis.title.size}{
    Changes the font size of the axes title. Defalt is 20.
    }
    \item{axis.text.size}{
    Changes the font size of the axes text labels. Default is 15.
    }
    \item{axis.text.angle}{
    Changes the angle of the X axis text labels. Default is 0 (horizontal).
    }
    \item{legend.title.size}{
    Changes the font size of the legend title. Default is 20.
    }
    \item{legend.text.size}{
    Changes the font size of the legend text labels. Default is 15.
    }
    \item{legend.key.size}{
    Changes the size of the legend keys (in points). Default is 30.
    }
}

\value{
    The ggplot object, or the vector of feature identifiers if multiple
    features are annotated to the gene symbol.
}

\references{
    \itemize{
    \item \href{http://cran.r-project.org/web/packages/ggplot2/index.html}{
    ggplot2} package.
    }
}

\author{
Kevin Rue-Albrecht
}

\section{Warning}{
    Common issues:
    \itemize{
    \item It may not be possible to produce plots where the combination of
    X-axis variable and grouping factor leaves too few replicates to compute
    a confidence interval for each X value. This is a limitation imposed by
    the ggplot2 package to produce proper statistics and confidence intervals.
    In such cases, it may be preferrable to use the
    expression_profiles_symbol() method.
    }
}

\seealso{
    Packages \code{\link[Biobase]{Biobase}} and \code{\link[ggplot2]{ggplot2}}
    , methods \code{\link[GOexpress:expression_plot]{expression_plot}} and
    \code{\link[GOexpress:GO_analyse]{GO_analyse}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# Expression by gene identifier (TNIP3)
expression_plot_symbol(gene_symbol="PIK3AP1", eSet=AlvMac, x_var="Timepoint",
                        result=AlvMac_results)

# Same gene, plotted by animal and grouped by treatment (merging time points)
expression_plot_symbol(gene_symbol="PIK3AP1", eSet=AlvMac, x_var="Animal",
                        result=AlvMac_results,
                        f="Treatment")

# Same gene, plotted by animal and grouped by time-point (merging treatments)
expression_plot_symbol(gene_symbol="PIK3AP1", eSet=AlvMac, x_var="Animal",
                        result=AlvMac_results,
                        f="Time")
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }