\name{GO_analyse}

\alias{GO_analyse}

\title{
Identifies gene ontologies clustering samples according to predefined factor.
}

\description{
Combines gene expression data with Gene Ontology (GO) annotations to rank and
visualise genes and GO terms enriched for genes best clustering predefined
groups of samples based on gene expression levels. 

This methods semi-automatically retrieves the latest information from Ensembl
using the \code{biomaRt} package, except if custom GO annotations are
provided. Custom GO annotations have two main benefits: firstly
they allow the analysis of species not supported in the Ensembl BioMart
server, and secondly they save time skipping calls to the Ensembl BioMart
server for species that are supported.
The latter also presents the possiblity of using an older release of the
Ensembl annotations, for example.

Using default settings, a random forest analysis
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
    mtry=ceiling(2*sqrt(nrow(eSet))), GO_genes=NULL, all_GO=NULL,
    all_genes=NULL, ...)
}

\arguments{
    \item{eSet}{
    \code{ExpressionSet} of the \code{Biobase} package including a
    gene-by-sample expression matrix in the \code{AssayData} slot, and a
    phenotypic information data-frame in the \code{phenodate} slot. In the
    expression matrix, row names are identifiers of expressed features, and
    column names are identifiers of the individual samples.
    In the phenotypic data-frame, row names are sample idenfifiers, column
    names are grouping factors and phenotypic traits usable for
    statistical tests and visualisation methods. 
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
    The Ensembl BioMart dataset identifier corresponding to the species
    studied. 
    If not specified and no custom annotations were provided, the method will
    attempt to automatically identify the
    adequate dataset from the first feature identifier in the dataset.
    Use data(prefix2dataset) to access a table listing valid choices.    
    }
    \item{microarray}{
    The identifier in the Ensembl BioMart corresponding to the microarray
    platform used. If not specified and no custom annotations were provided,
    the method will attempt to
    automatically identify the platform used from the first feature identifier
    in the dataset. 
    Use \code{data(microarray2dataset)} to access a table listing valid
    choices.
    }
    \item{method}{
    The statistical framework to score genes and gene ontologies.
    Either "randomForest" or "rf" to use the random forest algorithm, or
    alternatively either of "anova" or "a" to use the one-way ANOVA model.
    Default is "randomForest".
    }
    \item{rank.by}{
    Either of "rank" or "score" to chose the metric used to order the gene and
    GO term result tables.
    }
    \item{do.trace}{
    Only used if method="randomForest". If set to TRUE, gives a more verbose
    output as randomForest is run. If set to some integer, then running output
    is printed for every do.trace trees. Default is 100.
    }
    \item{ntree}{
    Only used if method="randomForest". Number of trees to grow. This should
    be set to a number large enough to ensure that every input row gets
    predicted at least a few times
    }
    \item{mtry}{
    Only used if method="randomForest". Number of features randomly sampled as
    candidates at each split. Default value is 2*sqrt(gene_count) which is
    approximately 220 genes for a dataset of 12,000 genes.
    }
    \item{GO_genes}{
    Custom annotations associating features present in the expression dataset
    to gene ontology identifiers. This must be provided as a data-frame of
    two columns, named \code{gene_id} and \code{go_id}. If provided, no call
    to the Ensembl BioMart server will be done, and arguments \code{all_GO}
    and \code{all_genes} should be provided as well, to enable all downstream
    features of \code{GOexpress}.
    An example is provided in \code{AlvMac_GOgenes}.
    }
    \item{all_GO}{
    Custom annotations used to annotate each GO identifier present in
    \code{GO_genes} with the ontology name (e.g. "apoptotic process") and
    namespace
    (i.e. "biological_process", "molecular_function", or
    "cellular_component").
    This must be provided as a data-frame containing
    at least one column named \code{go_id}, and preferably two more columns
    named \code{name_1006} and \code{namespace_1003} for consistency with
    the Ensembl BioMart. Supported alternative column headers are
    \code{name} and \code{namespace}.
    Respectively, \code{name} should be used to provide
    the description of the GO term, and \code{namespace} should contain
    one of "biological_process", "molecular_function" and
    "cellular_component". \code{name} is used to generate the title of
    ontology-based figured, and \code{namespace} is important to enable
    subsequent filtering of results by their corresponding value.
    An example is provided in \code{data(AlvMac_allGO)}.
    }
    \item{all_genes}{
    Custom annotations used to annotate each feature identifier in the
    expression dataset with the gene name or symbol (e.g. "TNF"), and an
    optional description. This must be provided as a data-frame containing at
    least a column named \code{gene_id} and preferably two more columns named
    \code{external_gene_name} and \code{description} for consistency with
    the Ensembl BioMart. A supported alternative header is \code{name}.
    While \code{external_gene_name} is important to enable
    subsequent visualisation of results by gene symbol, \code{description}
    is only displayed for readability of result tables.
    An example is provided in \code{data(AlvMac_allgenes)}.
    }
    \item{\dots}{
    Additional arguments passed on to the randomForest() method, if
    applicable.
    }
}

\value{
A list containing the results of the analysis. Some elements are specific to
the output of each analysis method.

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
    \item{subset}{The filters used to run the analysis only on a subet of the
    samples. NULL if no filter was applied.}

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
    
    It is suggested to use the \code{subset_scores()} function to subsequently
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
    
    If the column defining the factor (e.g. "Treatment") in \code{phenodata}
    is not an R factor, use 
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
# Load a local copy of annotations obtained from the Ensembl Biomart server
data(AlvMac_GOgenes)
data(AlvMac_allgenes)
data(AlvMac_allGO)

# Run the analysis on factor "Treatment",
# considering only treatments "MB" and "TB" at time-point "48H"
# using a local copy of annotations obtained from the Ensembl BioMart server
AlvMac_results <- GO_analyse(
    eSet=AlvMac, f="Treatment",
    subset=list(Time=c("48H"), Treatment=c("MB", "TB")),
    GO_genes=AlvMac_GOgenes, all_genes=AlvMac_allgenes, all_GO=AlvMac_allGO
    )

# Valid Ensembl BioMart datasets are listed in the following variable
data(prefix2dataset)

# Valid microarray= values are listed in the following variable
data(microarray2dataset)

\dontrun{
# Other valid but time-consuming examples:

# Run the analysis on factor "Treatment" including all samples
AlvMac_results <- GO_analyse(eSet=AlvMac, f="Treatment")

# Run the analysis on factor "Treatment" using ANOVA method
AlvMac_results <- GO_analyse(eSet=AlvMac, f="Treatment", method="anova")


# Syntax examples without actual data:

# To force the use of the Ensembl BioMart for the human species, use:
GO_analyse(eSet, f, biomart_dataset = "hsapiens_gene_ensembl")

# To force use of the bovine affy_bovine microarray annotations use:
GO_analyse(eSet, f, microarray = "affy_bovine")
}
}

\keyword{ GOexpress }
\keyword{ randomForest }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ clustering }
\keyword{ ontology }
