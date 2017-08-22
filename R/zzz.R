.onLoad <- function(libname, pkgname){
    packageStartupMessage(
"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ==== GOexpress is undergoing a progressive transition to S4! ====
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This will fundamentally alter the sequence of commands to perform an analysis,
yet facilitate the support and use of common classes implemented in other
Bioconductor packages (e.g., `DESeq2`, `edgeR`, `SummarizedExperiment`).

In particular, the main function `GO_analyse` is being split into separate
sub-methods to better control the individual steps of the workflow
(e.g., `randomForest`).

In addition, the semi-automated support for annotation through the Ensembl
`biomaRt` platform will be discontinued. Bioconductor annotation packages and
classes will be favoured (e.g., `org.Hs.eg.db`, `GSEABase::GOCollection`).

Please be mindful of messages and warnings emitted by the various functions.
Current functions will only be marked as 'Defunct' when appropriate
replacement methods are available.
"
    )
}
