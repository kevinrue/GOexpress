
# Test data ----

nrows <- 100; ncols <- 10
featureIds <- sprintf("F%03d", 1:100)
sampleIds <- sprintf("Sample%02d", 1:10)
mat <- matrix(
    runif(nrows * ncols), nrow=nrows, ncol=ncols,
    dimnames = list(featureIds, sampleIds)
)
pdata <- data.frame(
    group = rep(LETTERS[1:2], each = ncols/2),
    row.names = sampleIds
)
rRanges <- GRanges(
    rep(c("chr1", "chr2"), c(75, 25)),
    IRanges(seq(1:100)*100-99, width=100),
    strand=sample(c("+", "-"), 100, TRUE),
    feature_id=featureIds)

eset <- ExpressionSet(
    assayData=mat,
    phenoData = AnnotatedDataFrame(pdata)
)

rse <- SummarizedExperiment(
    assays=SimpleList(exprs=mat),
    rowRanges=rRanges, colData=DataFrame(pdata)
)

# Tests ----

## ExpressionSet
rf_eset <- randomForest(eset, "group")
## SummarizedExperiment
rf_rse <- randomForest(rse, "group")

test_that("randomForest works on ExpressionSet",{

    expect_s3_class(rf_eset, "randomForest")

})

test_that("randomForest works on (Ranged)SummarizedExperiment",{

    expect_s3_class(rf_rse, "randomForest")

})
