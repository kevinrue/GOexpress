\name{expression_profiles_symbol}

\alias{expression_profiles_symbol}

\title{
Plots the individual expression profile of a gene in samples series
}

\description{
This function will plot the expression profile of a gene in individual samples
series across a valid X-axis variable in the \code{phenodata}, while colouring
sample groups according to another variable in the \code{phenodata}, and 
using different line types according to yet another (or the same) variable
in the \code{phenodata}.
}

\usage{
expression_profiles_symbol(
    gene_symbol, result, eSet, x_var, seriesF, subset=NULL,
    colourF=result$factor, linetypeF=colourF, line.size=1.5,
    index=0, ylab="log2(cpm)", col.palette="Accent",
    col=brewer.pal(n=length(levels(pData(eSet)[,colourF])),
                    name=col.palette),
    lty=1:length(levels(pData(eSet)[,linetypeF])),
    titles=c(), title.size=2, axis.title.size=20,
    axis.text.size=15, axis.text.angle=0,
    legend.title.size=20, legend.text.size=15,
    legend.key.size=30)
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
    \item{seriesF}{
    A column name in \code{phenodata} to track samples of an identical series
    across the X-axis. The combination of \code{seriesF} and \code{x_var}
    should uniquely identify each sample in \code{eSet}.
    }
    \item{subset}{
    A named list to subset \code{eSet}. Names must be column names existing
    in colnames(pData(eSet)). Values must be vectors of values existing in
    the corresponding column of pData(eSet).
    }
    \item{colourF}{
    A column name in \code{phenodata} to colour the corresponding groups of
    samples. The factor specified in the initial \code{GO_analyse()} call is
    used by default.
    }
    \item{linetypeF}{
    A column name in \code{phenodata} to assign a different line type to the
    corresponding groups of samples. Default mirrors \code{colourF}.
    }
    \item{line.size}{
    The width of the plotted lines. Default is 1.5.
    }
    \item{index}{
    In the case where multiple gene or probeset identifiers are associated
    with the gene name, \code{index=2} will plot the expression profile of
    the second feature identifier alone, for instance.
    }
    \item{ylab}{
    Title of the Y-axis. Default is "log2(cpm)".
    }
    \item{col.palette}{
    A valid \code{RColorBrewer} palette name to fetch the colormap from.
    Default is palette "Accent".
    }
    \item{col}{
    A vector of color names or codes. The number of colors provided must match
    the number of levels of the factor \code{colourF}. If specified, overrides
    argument \code{col.palette}. Default colours are obtained from the 
    \code{col.palette}.
    }
    \item{lty}{
    A vector of numeric values corresponding to line types. The number of
    line types provided must match the number of levels of the factor
    \code{linetypeF}. Default line types are
    \code{1:length(levels(pData(eSet)[, linetypeF]))}.
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

\details{
In order to track and visualise individual sample series, each sample in the
\code{ExpressionSet} should be associated with a unique combination of
\code{seriesF} and \code{x_var}. This may require the generation of a new
factor in the \code{phenodata}, combining all experimental factors except
that plotted on the X-axis. See below for an example on the training dataset.
}

\value{
    The ggplot object.
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

\seealso{
    Packages \code{\link[Biobase]{Biobase}} and \code{\link[ggplot2]{ggplot2}}
    , methods
    \code{\link[GOexpress:expression_plot_symbol]{
    expression_plot_symbol}}
    and
    \code{\link[GOexpress:expression_plot_symbol]{
    expression_plot_symbol}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

AlvMac$Series = paste(AlvMac$Animal, AlvMac$Treatment, sep="_")

expression_profiles_symbol(gene_symbol="TNIP3", result = AlvMac_results,
        linetypeF="Animal", line.size=1.5,
        eSet=AlvMac, seriesF="Series", x_var = "Timepoint",
        title.size=1.5, legend.title.size=10, legend.text.size=10,
        legend.key.size=30)

expression_profiles_symbol(gene_symbol="TNIP3", result = AlvMac_results,
        linetypeF="Animal", line.size=1.5, lty=rep(1,10),
        eSet=AlvMac, seriesF="Series", x_var = "Timepoint",
        title.size=1.5, legend.title.size=10, legend.text.size=10,
        legend.key.size=30)
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
