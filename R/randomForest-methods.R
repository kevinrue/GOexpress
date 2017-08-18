
randomForest.ExpressionSet <- function(
    x, pheno, assay = "exprs", ..., do.trace = 100
){
    stopifnot(pheno %in% colnames(pData(x)))
    stopifnot(assay %in% names(assayData(x)))
    rf <- .randomForest(
        t(assayData(x)[[assay]]),
        pData(x)[,pheno],
        ...,
        do.trace = do.trace
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

randomForest.SummarizedExperiment <- function(
    x, pheno, assay = "exprs", ..., do.trace = 100
){
    stopifnot(pheno %in% colnames(colData(x)))
    stopifnot(assay %in% names(assays(x)))
    rf <- .randomForest(
        t(assay(x, assay)),
        colData(x)[,pheno],
        ...,
        do.trace = do.trace
    )
    rf$call$x <- substitute(x)
    rf$call$y <- substitute(pheno)
    return(rf)
}

# Default method
# x: matrix (predictor = features as columns)
# pdata: factor (length(pdata) == nrow(x))
.randomForest <- function(x, pheno, ..., do.trace = 100){
    stopifnot(nrow(x) >= 4)
    stopifnot(is.factor(pheno))

    tablePheno <- table(pheno)
    stopifnot(all(tablePheno != 0))

    levelsPheno <- names(tablePheno)
    message(length(levelsPheno), " levels")
    for (i in seq_along(levelsPheno)){
        message("  - ", levelsPheno[i], " : ", tablePheno[i], " samples")
    }
    message("")

    return(randomForest(
        x, y = pheno, importance = TRUE, do.trace = do.trace, ...
    ))
}
