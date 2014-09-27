\name{plot_design}

\alias{plot_design}

\title{
Plot Univariate Effects for genes associated with a Gene Ontology
}

\description{
Successively plots univariate effects of one or more
\code{\link[base:factor]{factors}}, typically for a designed experiment
as analyzed by \code{\link[stats:aov]{aov}}().
}

\usage{
plot_design(
    go_id, result, eSet,
    factors=colnames(pData(eSet)), main="", main.Lsplit=NULL, ...)
}

\arguments{
    \item{go_id}{
    A Gene Ontology (GO) identifier represented by at least one gene
    in the dataset.
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
    \item{factors}{
    A set of column names from \code{phenodata}. Each of these values will be
    represented on the X-axis to investigate its effect on the average
    expression of a given genes for each level of that factor.
    }
    \item{main}{
    Changes the main title of the plots.
    }
    \item{main.Lsplit}{
    Number of characters after which a new-line character will be inserted in
    the main title. If this would occur within a word, the new-line character
    will be inserted before this word. Default is NULL, leaving the title on a
    single line.
    }
    \item{\dots}{
    Additional arguments which will be passed on to the \code{plot.design}()
    function.
    }
}

\value{
    The output of the \code{plot.design}() function.
}

\author{
    Kevin Rue-Albrecht
}

\seealso{
    Method \code{\link[graphics:plot.design]{plot.design}}.
}

\examples{
# load the sample output data
data(AlvMac_results)

# Univariate plot
plot_design(go_id="GO:0034142", eSet=AlvMac, result=AlvMac_results)
}

\keyword{ GOexpress }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
