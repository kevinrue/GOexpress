\name{GO_analyse}

\alias{GO_analyse}

\title{
Identifies gene ontologies clustering samples according to predefined factor.
}

\description{
Combines local gene expression data with online Gene Ontology (GO) resources
to visualise and rank ontologies enriched for genes best clustering predefined
groups of samples based on gene expression levels. 

Using the \code{biomaRt} package, this method semi-automatically retrieves
gene ontology annotations from Ensembl for the species corresponding to the
expression dataset  prior to performing the analysis. A random forest analysis
is performed to evaluate the ability of each gene to cluster samples according
to a predefined grouping factor (one-way ANOVA available as an alteranative).
Each GO term is scored and ranked according to the average rank
(alternatively, average power) of all associated genes to cluster the
samples according to the factor. The ranked list of GO terms is returned,
with tools allowing to visualise the statistics on a gene- and ontology-basis.
}

\usage{
GO_analyse(
    eSet, f, subset=NULL, biomart_dataset="", microarray="",
    method="randomForest", rank.by="rank", do.trace=100, ntree=1000,
    mtry=ceiling(2*sqrt(nrow(eSet))), 
    ...)
}

\arguments{
    \item{eSet}{
    \code{ExpressionSet} of the \code{Biobase} package including a
    gene-by-sample expression matrix in the \code{AssayData} slot, and a
    phenotypic information data-frame in the \code{phenodate} slot. In the
    expression matrix, row names are Ensembl gene identifiers or probeset
    identifiers, and column names are sample identifiers. In the phentypic
    data-frame, row names are sample idenfifiers, column names are grouping
    factors and phenotypic traits usable for the one-way ANOVA. 
    }
    \item{f}{
    A column name in \code{phenodata} used as the grouping factor for the
    analysis.
    }
    \item{subset}{
    A named list to subset \code{eSet} for the analysis. Names must be
    column names existing in colnames(pData(eSet)). Values must be vectors of
    values existing in the corresponding column of pData(eSet). The original
    ExpressionSet will be left unchanged.
    }
    \item{biomart_dataset}{
    The BioMart Ensembl dataset corresponding to the species studied.
    Automatically detected if not specified. Use data(prefix2dataset) for
    valid values.
    }
    \item{microarray}{
    The identifier in the Ensembl BioMart corresponding to the microarray
    platform used. If not specified and Ensembl gene identifiers are not
    detected in the expression dataset, the method will attempt to
    automatically identify the platform used from the first probeset in the
    dataset. See \code{data(microarray2dataset)} for valid values.
    }
    \item{method}{
    Either "randomForest" or "rf" to use the random forest framework, or
    alternatively either of "anova" or "a" to use the one-way ANOVA model.
    Default is "randomForest".
    }
    \item{rank.by}{
    Either of "rank" or "score" to chose the metric used to order the gene and
    GO term result tables.
    }
    \item{do.trace}{
    Only used if method="randomForest". If set to TRUE, give a more verbose
    output as randomForest is run. If set to some integer, then running output
    is printed for every do.trace trees. Default is 100.
    }
    \item{ntree}{
    Only used if method="randomForest". Number of trees to grow. This should
    not be set to too small a number, to ensure that every input row gets
    predicted at least a few times
    }
    \item{mtry}{
    Only used if method="randomForest". Number of features randomly sampled as
    candidates at each split. Default value is 2*sqrt(gene_count) which is
    approximately 220 genes for a dataset of 12,000 genes.
    }
    \item{\dots}{
    Additional arguments passed on to the randomForest() method.
    }
}

\value{
A list containing the results of the analysis. Some elements are specific to
each analysis method.

Core elements:
    \item{GO}{A table ranking all GO terms related to genes in the expression
    dataset based on the average ability of their related genes to cluster the
    samples according to the predefined grouping factor.}
    \item{mapping}{The table mapping genes present in the dataset to GO
    terms.}
    \item{genes}{A table ranking all genes present in the expression dataset
    based on their ability to cluster the samples according to the predefined
    grouping factor.}
    \item{factor}{The predefined grouping factor.}
    \item{method}{The statistical framework used.}

Random Forest additional elements:
    \item{ntree}{Number of trees grown.}
    \item{mtry}{Number of variables randomly sampled as candidates at each
    split.}
    
One-way ANOVA does not have additional arguments.
}

\details{
    The current scoring functions strongly favor GO terms associated with
    fewer genes at the top of the ranking. This bias may actually be seen as
    a valuable feature which enables the user to browse through GO terms of
    increasing "granularity", i.e. associated with increasingly large sets of
    genes, although consequently being increasingly vague and blurry (e.g.
    "protein binding" molecular function associated with over 6,000 genes).
    
    It is suggested to use the \code{subset_scores()} function to
    filter out GO terms with fewer than 5+ genes associated with it. Indeed,
    those GO terms are more sensitive to outlier genes as they were scored on
    the average of a handful of genes.
}

\author{
    Kevin Rue-Albrecht
}

\section{Warning}{
    Make sure that the factor \code{f} is an actual factor in the R language
    meaning. This is important for the underlying statistical framework to
    identify the groups of samples defined by their level of this factor.
    
    If the column defining the factor (e.g. "Treatment") in code{phenodata} is
    not an R factor, use 
    \code{pData(targets)$Treatment = factor(pData(targets)$Treatment)}
    to convert the character values into an actual R factor with appropriate
    levels.
}

\seealso{
    Methods
    \code{\link[GOexpress:subset_scores]{subset_scores}},
    \code{\link[biomaRt:getBM]{getBM}},
    \code{\link[randomForest:randomForest]{randomForest}}, 
    and \code{\link[stats:oneway.test]{oneway.test}}.
}

\examples{
# Load example data subset
data(AlvMac)

# Run the analysis on factor "Treatment",
# considering only treatments "MB" and "TB" at time-point "48H"
# Factor 'f' must be one of the column names in pData(targets)
AlvMac_results <- GO_analyse(eSet=AlvMac, f="Treatment",
    subset=list(Time=c("48H"), Treatment=c("MB", "TB")))

\dontrun{
# Syntax examples without actual data:

# To force the use of the Ensembl BioMart for the human species, use:
GO_analyse(eSet, f, biomart_dataset = "hsapiens_gene_ensembl")

# To force use of the bovine affy_bovine microarray annotations use:
GO_analyse(eSet, f, microarray = "affy_bovine")

# Valid Ensembl BioMart datasets are listed in the following variable
data(prefix2dataset)

# Valid microarray= values are listed in the following variable
data(microarray2dataset)


# Additional time-consuming example usages, can be run manually
# with test data:

# Run the analysis on factor "Treatment" including all samples
# Factor must be one of the column names in pData(targets)
AlvMac_results <- GO_analyse(eSet=AlvMac, f="Treatment")

# Run the analysis on factor "Treatment" using ANOVA method
AlvMac_results <- GO_analyse(eSet=AlvMac, f="Treatment", method="a")
}
}

\keyword{ GOexpress }
\keyword{ randomForest }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }
