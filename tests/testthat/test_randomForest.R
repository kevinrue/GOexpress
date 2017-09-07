
stopifnot(requireNamespace("DESeq2"))

# Test data ----

nrows <- 100; ncols <- 10
featureIds <- sprintf("F%03d", 1:100)
sampleIds <- sprintf("Sample%02d", 1:10)
mat <- matrix(
    runif(nrows * ncols), nrow=nrows, ncol=ncols,
    dimnames=list(featureIds, sampleIds)
)
pdata <- data.frame(
    group=rep(LETTERS[1:2], each=ncols/2),
    row.names=sampleIds
)
rRanges <- GRanges(
    seqnames=rep(c("chr1", "chr2"), c(75, 25)),
    ranges=IRanges(seq(1:100)*100-99, width=100),
    strand=sample(c("+", "-"), 100, TRUE),
    feature_id=featureIds)

## ExpressionSet
eset <- ExpressionSet(
    assayData=mat,
    phenoData = AnnotatedDataFrame(pdata)
)
## SummarizedExperiment
rse <- SummarizedExperiment(
    assays=SimpleList(exprs=mat),
    rowRanges=rRanges, colData=DataFrame(pdata)
)
## DESeqDataSet
dds <- DESeq2::makeExampleDESeqDataSet(m=10, betaSD=1)
dds <- DESeq2::estimateSizeFactors(dds)
## DGEList
dgel <- edgeR::DGEList(counts=dds@assays[["counts"]], group=pdata$group)
rownames(dgel) <- rownames(dds)
dgel <- edgeR::calcNormFactors(dgel)

# Tests ----

## ExpressionSet
rf_eset <- randomForest(eset, "group")
## SummarizedExperiment

## DESeqDataSet

## DESeqTransform

test_that("randomForest works on ExpressionSet",{

    expect_s3_class(rf_eset, "randomForest")

})

test_that("randomForest works on (Ranged)SummarizedExperiment",{

    rf_rse <- randomForest(rse, "group", verbose=TRUE)

    expect_s3_class(rf_rse, "randomForest")

})

test_that("randomForest works on DESeqDataSet",{

    rf_dds <- randomForest(dds, "condition")

    expect_s3_class(rf_dds, "randomForest")

})

test_that("randomForest works on DESeqTransform",{

    rf_drlog <- randomForest(
        DESeq2::rlog(dds), "condition")
    rf_dvst <- randomForest(
        DESeq2::varianceStabilizingTransformation(dds), "condition")

    expect_s3_class(rf_drlog, "randomForest")
    expect_s3_class(rf_dvst, "randomForest")

})

test_that("randomForest works on DGEList",{

    rf_dgel <- randomForest(dgel, "group")

    expect_s3_class(rf_dgel, "randomForest")

})
