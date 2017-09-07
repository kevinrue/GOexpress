
randomForest.ExpressionSet <- function(
    x, pheno, assay="exprs", ..., do.trace=100, verbose=FALSE
){
    stopifnot(pheno %in% colnames(pData(x)))
    stopifnot(assay %in% names(assayData(x)))
    rf <- .randomForest(
        t(assayData(x)[[assay]]),
        pData(x)[,pheno],
        ...,
        do.trace=do.trace,
        verbose=verbose
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

randomForest.SummarizedExperiment <- function(
    x, pheno, assay="exprs", ..., do.trace=100, verbose=FALSE
){
    stopifnot(pheno %in% colnames(colData(x)))
    stopifnot(assay %in% names(assays(x)))
    rf <- .randomForest(
        t(assay(x, assay)),
        colData(x)[,pheno],
        ...,
        do.trace=do.trace,
        verbose=verbose
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

randomForest.DESeqDataSet <- function(
    x, pheno, normalized=TRUE, ..., do.trace=100, verbose=FALSE
){
    stopifnot(requireNamespace("DESeq2"))
    stopifnot(pheno %in% colnames(colData(x)))
    rf <- .randomForest(
        t(DESeq2::counts(x, normalized=normalized)),
        colData(x)[,pheno],
        ...,
        do.trace=do.trace,
        verbose=verbose
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

randomForest.DESeqTransform <- function(
    x, pheno, ..., do.trace=100, verbose=FALSE
){
    stopifnot(pheno %in% colnames(colData(x)))
    rf <- .randomForest(
        t(assays(x)[[1]]),
        colData(x)[,pheno],
        ...,
        do.trace=do.trace,
        verbose=verbose
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

randomForest.DGEList <- function(
    x, pheno, normalized.lib.sizes=TRUE, ..., do.trace=100, verbose=FALSE
){
    stopifnot(requireNamespace("edgeR"))
    stopifnot(pheno %in% colnames(x[["samples"]]))
    rf <- .randomForest(
        t(edgeR::cpm(x, normalized.lib.sizes=normalized.lib.sizes)),
        x[["samples"]][,pheno],
        ...,
        do.trace=do.trace,
        verbose=verbose
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

# Default method
# x: matrix (predictor = features as columns)
# pdata: factor (length(pdata) == nrow(x))
.randomForest <- function(x, pheno, ..., do.trace=100, verbose=FALSE){
    stopifnot(nrow(x) >= 4)

    # Untested: pheno might be a numeric value (for regression)
    if (is.factor(pheno)){
        # Ensure all phenotype levels are populated
        tablePheno <- table(pheno)
        stopifnot(all(tablePheno != 0))
        if (verbose){
            levelsPheno <- names(tablePheno)
            message(length(levelsPheno), " levels")
            for (i in seq_along(levelsPheno)){
                message("  - ",levelsPheno[i]," : ",tablePheno[i]," samples")
            }
            message("")
        }
    }

    return(randomForest(
        x, y=pheno, importance=TRUE, do.trace=do.trace, ...
    ))
}
