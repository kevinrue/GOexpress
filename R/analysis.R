
GO_analyse <- function(
    eSet, f, subset=NULL, biomart_dataset="", microarray="",
    method="randomForest", rank.by="rank", do.trace=100, ntree=1000,
    mtry=ceiling(2*sqrt(nrow(eSet))), GO_genes=NULL, all_GO=NULL,
    all_genes=NULL, FUN.GO=mean, ...){
    if (class(eSet) != "ExpressionSet"){
        stop("eSet must be an ExpressionSet of the Biobase package.")
    }
    # if less than 4 genes in data will cause mtry larger than number of genes
    # which is then impossible.
    # However, who uses a transcriptomics dataset of 4 genes?
    if (nrow(eSet) < 4 && method %in% c("randomForest","rf")){
        stop("Too few genes in dataset: ", nrow(eSet))
    }
    # If the user gave an invalid column name for the factor
    if (!f %in% colnames(pData(eSet))){
        stop("Invalid column name (absent from phenoData): ", f)
    }
    # The random forest requires the factor (f) to be an actual R factor
    if (!any(class(pData(eSet)[,f]) == "factor")){
        stop("pData(eSet)[,f] must be an actual R factor.
    See help(factor)")
    }
    # The random forest requires all levels of the factor to have at least one
    # sample, otherwise it crashes
    if (any(table(pData(eSet)[,f]) == 0)){
        stop("One level of pData(eSet)[,f] has no correspondig sample.
    Please remove that level.")
    }
    # Subset the data to the given values of the given factors, if existing
    if (!is.null(subset)){
        eSet <- subEset(eSet=eSet, subset=subset)
    }
    # If the user gave an invalid method name (allow specific abbreviations)
    if (!method %in% c("randomForest", "anova", "rf", "a")){
        stop("Invalid method: ", method)
    }
    # If the user gave an invalid ranking method for GO terms and genes
    if (!rank.by %in% c("rank","score")){
        stop("Invalid ranking method: ", rank.by)
    }
    # Remove rows with NA is the microarray2dataset table. Those rows
    # correspond to datasets (species) without any microarray
    microarray2dataset.clean <- microarray2dataset[
        !is.na(microarray2dataset$unique),]
    # If the user did not provide custom annotations
    if (is.null(GO_genes)){
        # if the user did not give a dataset name
        if (biomart_dataset == ""){
            # if the user did not give a microarray value
            if (microarray == ""){
                # automatically detect both
                # fetch the first gene id in the given expression dataset
                sample_gene <- rownames(eSet)[1]
                cat(
                    "First feature identifier in dataset:", sample_gene,
                    fill=TRUE
                    )
                # Try to find an appropriate biomaRt Ensembl dataset from the
                # gene prefix
                mart <- mart_from_ensembl(sample_gene)
                # if the gene id has not an identifiable Ensembl id prefix
                if (!class(mart) == "Mart"){
                    # Try to find an appropriate biomaRt microarray dataset
                    # from the gene prefix
                    microarray_match <- microarray_from_probeset(
                        sample_gene,
                        microarray2dataset.clean)
                    # if the gene id has an identifiable microarray id prefix
                    if (!is.null(nrow(microarray_match))){
                        # connect to biomart and set the microarray variable
                        cat("Looks like microarray data.", fill=TRUE)
                        cat(
                            "Loading detected dataset",
                            microarray_match$dataset,
                            "for detected microarray",
                            microarray_match$microarray, "...", fill=TRUE
                            )
                        microarray <- microarray_match$microarray
                        biomart_dataset <- microarray_match$dataset
                        mart <- useMart(
                            biomart="ensembl",
                            dataset=biomart_dataset
                            )
                    }
                    # if the gene id does not have an identifiable microarray
                    # gene id prefix
                    else{
                        # stop the program and throw an error
                        stop("Cannot guess origin of dataset:
        Please use \"biomart_dataset=\" and/or \"microarray=\" arguments.")
                    }
                }
                # if the gene id has an identifiable Ensembl gene id prefix
                # then, the connection to the mart is already established
                # leave microarray to "" to imply that we don't work with
                # microarray data
            }
            # if the user gave a microarray name
            else{
                # check if it exists
                if (!microarray %in% microarray2dataset.clean$microarray){
                    stop(
                        "Invalid microarray value. ",
                        "See data(microarray2dataset)"
                        )
                }
                # if it is unique to a dataset (some microarray have the same
                # column name
                if(sum(
                    microarray2dataset.clean$microarray == microarray
                    ) == 1){
                    biomart_dataset <- microarray2dataset.clean[
                        microarray2dataset.clean$microarray == microarray,
                        "dataset"]
                    cat(
                        "Loading requested microarray", microarray,
                        "from detected dataset", biomart_dataset, "...",
                        fill=TRUE
                        )
                    mart <- useMart(
                        biomart="ensembl",
                        dataset=biomart_dataset)
                    # Leave microarray to the current valid value
                }
                # if the microarray does not exist in the dataset
                else if(
                    sum(
                        microarray2dataset.clean$microarray == microarray
                        ) == 0){
                    stop(
                        "Microarray name not recognised. ",
                        "See data(microarray2dataset)."
                        )
                }
                # if the microarray name exists in multiple datasets
                else{
                    cat("Multiple datasets possible:", fill=TRUE)
                    print(
                        microarray2dataset.clean[
                        microarray2dataset.clean$microarray == microarray,
                        c("dataset", "microarray")]
                        )
                    stop(
                        "Cannot guess dataset. ",
                        "Please use \"biomart_dataset=\" argument."
                        )
                }
                }
        }
        # if the user gave a biomart_dataset value
        else{
            # Check that it exists
            if (!biomart_dataset %in% prefix2dataset$dataset){
                stop(
                    "Invalid biomart_dataset value. ",
                    "See data(prefix2dataset)"
                    )
            }
            cat("Using biomart dataset", biomart_dataset, fill=TRUE)
            # if the user did not give a microarray name
            if (microarray == ""){
                # Check if looks like microarray
                # fetch the first gene id in the given expression dataset
                sample_gene <- rownames(eSet)[1]
                cat(
                    "First feature identifier in dataset:",
                    sample_gene,
                    fill=TRUE
                    )
                microarray_match <- microarray_from_probeset(
                    sample_gene,
                    microarray2dataset.clean
                    )
                # if the data matches a known microarray pattern
                if (!is.null(nrow(microarray_match))){
                    # connect to biomart and set the microarray variable
                    cat("Looks like microarray data.", fill=TRUE)
                    # if the dataset/microarray pair exists
                    if (microarray_match$dataset == biomart_dataset){
                        cat(
                            "Loading annotations for detected microarray",
                            microarray_match$microarray,
                            "for requested dataset", biomart_dataset, "...",
                            fill=TRUE
                            )
                        mart <- useMart(
                            biomart="ensembl",
                            dataset=biomart_dataset
                            )
                        microarray <- microarray_match$microarray
                    }
                    # if the dataset/microarray pair does no exist
                    else{
                        # The dataset exists, the data matches a microarray
                        # but not a microarray of the dataset
                        cat(
                            "Detected microarray",
                            microarray_match$microarray,
                            "inexisting in requested dataset",
                            biomart_dataset,
                            ". Possible datasets are:", fill=TRUE)
                        return(microarray_match)
                    }
                }
                # if the data does not match a microarray pattern
                else{
                    cat(
                        sample_gene,
                        "feature identifier in expression data cannot",
                        "be resolved to a microarray. Assuming Ensembl gene",
                        "identifiers.",
                        fill=TRUE)
                }
                # If it does not look like microarray
                # assume it is Ensembl annotations
                # therefore do nothing more
                # in both cases load the requested mart dataset
                cat(
                    "Loading requested dataset", biomart_dataset, "...",
                    fill=TRUE
                    )
                mart <- useMart(biomart="ensembl", dataset=biomart_dataset)
                }
            # if the user gave a microarray name
            else{
                # Check that the pair dataset/microarray exists
                if (!biomart_dataset %in% microarray2dataset.clean[
                    microarray2dataset.clean$microarray == microarray,
                    "dataset"]){
                    stop(
                        "There is no microarray ", microarray,
                        " in dataset ", biomart_dataset
                        )
                }
                cat(
                    "Loading requested microarray", microarray,
                    "from requested biomart dataset", biomart_dataset,
                    fill=TRUE
                    )
                mart <- useMart(biomart="ensembl", dataset=biomart_dataset)
            }
        }
        print(mart)
    }
    else {
        if (biomart_dataset != "" || microarray != ""){
            warning(
                "Non-NULL GO_genes argument: Ignoring 'biomart_dataset' ",
                "and 'microarray' arguments."
                )
            biomart_dataset <- ""
            microarray <- ""
        }
        mart <- NULL
    }
    # If working with custom gene identifiers
    if (!is.null(GO_genes)){
        if (!all(colnames(GO_genes) == c("gene_id", "go_id"))){
            stop(
                "Column names of custom GO_genes must be ",
                "c(\"gene_id\", \"go_id\")"
                )
        }
        cat("Using custom GO_genes mapping ...", fill=TRUE)
    }
    # If working with BioMart-compatible feature identifiers
    else{
        # if working with Ensembl gene identifiers
        if (microarray == ""){
            # Prepare a mapping table between gene identifiers and GO terms
            cat(
                "Fetching ensembl_gene_id/go_id mappings from BioMart ...",
                fill=TRUE)
            GO_genes <- getBM(
                attributes=c("ensembl_gene_id", "go_id"),
                mart=mart)
        }
        # if working with microarray probesets
        else{
            # Prepare a mapping table between feature identifiers and GO terms
            cat(
                "Fetching probeset/go_id mappings from BioMart ...",
                fill=TRUE
                )
            GO_genes <- getBM(attributes=c(microarray, "go_id"), mart=mart)
        }
    }
    # Rename the first column which could be ensembl_id or probeset_id
    colnames(GO_genes)[1] <- "gene_id"
    # Remove over 1,000 rows where the go_id is ""
    GO_genes <- GO_genes[GO_genes$go_id != "",]
    # Remove rows where the gene_id is "" (happens)
    GO_genes <- GO_genes[GO_genes$gene_id != "",]    
    # Print how many features are present in the map
    cat(
        length(intersect(rownames(eSet), unique(GO_genes$gene_id))),
        "features from ExpressionSet found in the mapping table.",
        fill=TRUE
        )
    # If the user did not provide custom GO description
    # Note that users can provide custom annotations even if BioMart was used
    # to fetch gene-GO mappings
    if (is.null(all_GO)){
        if (is.null(mart)){
            warning(
                "No custom GO terms descriptions provided. ",
                "The \"namespace\" filter of the subset_scores method will ",
                "not be available."
                )
            all_GO <- data.frame(
                go_id=unique(GO_genes[,"go_id"]),
                name_1006=NA,
                namespace_1003=NA
                )
        } else{
            # Prepare a table of all the GO terms in BioMart (even if no gene
            # is annotated to it)
            cat("Fetching GO_terms description from BioMart ...", fill=TRUE)
            all_GO <- getBM(
                attributes=c("go_id", "name_1006", "namespace_1003"),
                mart=mart
            )
        }
    }
    # If working with custom GO descriptions
    else{
        if (! "go_id" %in% colnames(all_GO)){
            stop("Mandatory column \"go_id\" not found in custom all_GO")
        }
        if (! "name_1006" %in% colnames(all_GO)){
            # Allow the header "name" but internally convert it to name_1006
            if ("name" %in% colnames(all_GO)){
                colnames(all_GO)[colnames(all_GO) == "name"] <- "name_1006"
            }
            # else if could allow more headers
            else {
                warning(
                    "We encourage the use of a \"name\" column describing",
                    "the GO term in all_GO."
                    )
            }
        }
        if (! "namespace_1003" %in% colnames(all_GO)){
            # Allow the header "namespace" but internally convert it
            if ("namespace" %in% colnames(all_GO)){
                colnames(all_GO)[
                    colnames(all_GO) == "namespace"
                    ] <- "namespace_1006"
            }
            # else if could allow more headers
            else {
                warning(
                    "We encourage the use of a 'namespace' column",
                    "describing the ontology in all_GO",
                    "(e.g. 'biological_process')."
                )
            }
        }
        cat("Using custom GO terms description ...", fill=TRUE)
    }
    # Remove the GO terms which is ""
    all_GO <- all_GO[all_GO$go_id != "",]
    # Run the analysis with the desired method
    cat("Analysis using method", method ,"on factor", f,"for", nrow(eSet),
        "genes. This may take a few minutes ...", fill=TRUE)
    if (method %in% c("randomForest", "rf")){
        ## Similarly to the previous anova procedure (see below)
        # Run the randomForest algorithm
        rf <- randomForest(
            x=t(exprs(eSet)), y=pData(eSet)[,f], importance=TRUE,
            do.trace=do.trace, ntree=ntree, mtry=mtry,
            ...)
        # Save the importance value used as score for each gene in a
        # data.frame
        res <- data.frame("Score" = importance(rf)[,"MeanDecreaseGini"])
        # In the output variable, write the full method name
        method <- "randomForest"
    }
    else if (method %in% c("anova", "a")){
        # A vectorised calculation the F-value of an ANOVA used as score for
        # each gene
        res <- data.frame(
            "Score" = apply(
                X=exprs(eSet),
                MARGIN=1,
                FUN=function(x){
                    oneway.test(
                        formula=expr~group,
                        data = cbind(
                            expr=x, group=pData(eSet)[,f]
                            )
                        )$statistic
                    }
                )
            )
        # In the output variable, write the full method name
        method <- "anova"
    }
    # Calculate the rank of each gene based on their score
    res$Rank <- rank(-res$Score, ties.method="min")
    # Summary statistics by GO term
    ## Merge the table mapping GOterm to genes with the score of each gene,
    ## twice:
    # - First merge the tables while keeping all gene/GO mappings, even for
    # genes absent of the dataset. This will allow average F values to be
    # calculated on the basis of all ensembl genes annotated to the GO term,
    # even if not in the dataset (genes absent will be given score of 0 and
    # rank of (number of genes in eSet + 1)
    # Merge GO/gene mapping with randomForest results
    GO_gene_score_all <- merge(
        x=GO_genes, y=res, by.x="gene_id", by.y="row.names", all.x=TRUE
        )
    # Replace NAs (genes absent from dataset but present in biomaRt) by 0
    # (minimal valid value)
    GO_gene_score_all[is.na(GO_gene_score_all$Score), "Score"] <- 0
    # In addition, all genes absent from the dataset are assigned a value =
    # 1 + (the maximum rank of the genes in the dataset)
    GO_gene_score_all[is.na(GO_gene_score_all$Rank), "Rank"] <-
        nrow(eSet) + 1
    # - Second, merge the tables keeping only the genes present in the dataset
    # This will be used to count how many genes are present in dataset for
    # each GO term
    GO_gene_score_data <- merge(
        x=GO_genes, y=res, by.x="gene_id", by.y="row.names")
    # Results can now be summarised by aggregating rows with same GOterm
    # Appends gene annotations to rows of res
    # If the user did not provide gene annotations
    # Note that users can provide custom annotations even if BioMart was used
    # to fetch gene-GO mappings
    if (is.null(all_genes)){
        if (is.null(mart)){
            warning(
                "No custom gene descriptions provided. ",
                "Visualisation methods using gene symbols will not be ",
                "available."
            )
            all_genes <- data.frame(
                gene_id=unique(rownames(eSet)),
                external_gene_name=NA,
                description=NA
            )
        }
        else{
            cat("Fetching gene description from BioMart ...", fill=TRUE)
            # if working with Ensembl gene identifiers
            if (microarray == ""){
                all_genes <- getBM(
                    attributes=c(
                        "ensembl_gene_id",
                        "external_gene_name", # since Ensembl release 76
                        "description"
                    ),
                    filters="ensembl_gene_id",
                    values=rownames(res), mart=mart)
            }
            # if working with microarray probesets
            else{
                # Prepare a mapping table between probesets and annotations
                all_genes <- getBM(
                    attributes=c(
                        microarray,
                        "external_gene_name", # since Ensembl release 76
                        "description"
                    ),
                    filters=microarray,
                    values=rownames(res),
                    mart=mart)
            }
            # Rename the first column which could be ensembl_id or probeset_id
            colnames(all_genes)[1] <- "gene_id"
        }
    }
    # If working with custom gene identifiers
    else{
        if (! "gene_id" %in% colnames(all_genes)){
            stop("Mandatory column \"gene_id\" not found in custom all_genes")
        }
        if (! "external_gene_name" %in% colnames(all_genes)){
            # Allow the header "name" but internally convert it
            if ("name" %in% colnames(all_genes)){
                colnames(all_genes)[
                    colnames(all_genes) == "name"
                    ] <- "external_gene_name"
            }
            else if ("external_gene_id" %in% colnames(all_genes)){
                colnames(all_genes)[
                    colnames(all_genes) == "external_gene_id"
                    ] <- "external_gene_name"
            }
            # "else if" could allow more synonym headers
            else {
                warning(
                    "We encourage the use of a \"name\" column describing",
                    "the gene in all_genes"
                )
            }
        }
        cat("Using custom gene descriptions ...", fill=TRUE)
    }
    genes_score <- merge(
        x=res,
        all.x=TRUE,
        y=all_genes,
        by.x="row.names",
        by.y="gene_id"
        )
    # In the case of microarray, probesets can be annotated to multiple
    # gene symbols. For each probeset, keep only the first row (we're
    # losing one possible annotation here). But otherwise, we will not be
    # able to make unique row names for each feature measured.
    # Moreover, keeping the same probeset twice (because of two annotated
    # symbols) would affect the averaging of scores for GO terms.
    # Each probeset should only be there once anyway
    genes_score <- genes_score[ !duplicated(genes_score$Row.names), ]
    # Put the feature identifier back as the row name
    rownames(genes_score) <- genes_score$Row.names
    genes_score$Row.names <- NULL
    cat("Merging score into result table ...", fill=TRUE)
    # Total number of genes in the dataset annotated to each GO term
    GO_scores <- merge(
        x=aggregate(
            gene_id~go_id,
            data=GO_gene_score_data,
            FUN=length
            ),
        y=all_GO,
        by="go_id",
        all.y=TRUE
        )
    colnames(GO_scores)[2] <- "data_count"
    GO_scores[is.na(GO_scores$data_count), "data_count"] <- 0
    # Total number of genes annotated to each GO term in BioMart (not
    # necessarily in dataset)
    GO_scores <- merge(
        x=aggregate(
            Score~go_id,
            data=GO_gene_score_all,
            FUN=length
            ),
        y=GO_scores,
        by="go_id",
        all.y=TRUE)
    colnames(GO_scores)[2] <- "total_count"
    GO_scores[is.na(GO_scores$total_count), "total_count"] <- 0
    # Average score (denominator being the total of genes by GO term in
    # BioMart) being tested
    GO_scores <- merge(
        x=aggregate(
            Score~go_id,
            data=GO_gene_score_all,
            FUN=FUN.GO
            ),
        y=GO_scores,
        by="go_id",
        all.y=TRUE
        )
    colnames(GO_scores)[2] <- "ave_score"
    GO_scores[is.na(GO_scores$ave_score), "ave_score"] <- 0
    ## Average rank (denominator being the total of genes by GO term in
    # BioMart) being tested
    # (+) robust for GO terms with several genes
    GO_scores <- merge(
        x=aggregate(
            Rank~go_id,
            data=GO_gene_score_all,
            FUN=FUN.GO
            ),
        y=GO_scores,
        by="go_id",
        all.y=TRUE
        )
    colnames(GO_scores)[2] <- "ave_rank"
    # Any GO term without associated gene (NA value) is givne worst rank + 1
    GO_scores[is.na(GO_scores$ave_rank), "ave_rank"] <- max(
        GO_scores$ave_rank, na.rm=TRUE
        ) + 1
    # Notes of other summary metrics tested:
    ## Sum: (-) biased toward general GO terms annotated for many thousands
    ## of genes (e.g. "protein binding")
    ## Max.F.values: (+) insensitive to number of genes annotated for GO term
    #                (-) many GO terms sharing the same gene are tied
    #                (-) biased by outliers
    # Most top ranked GO terms according to the average F value contain a
    # single gene, but this bias can easily be attenuated by filtering for GO
    # terms with a minimal number of genes
    if (rank.by == "rank"){
        # Rank the genes by increasing rank
        genes_score <- genes_score[order(genes_score$Rank),]
        # Rank the GO terms by decreasing average rank
        GO_scores <- GO_scores[order(GO_scores$ave_rank),]
    }
    else if (rank.by == "score"){
        # Rank the genes by increasing rank
        genes_score <- genes_score[order(genes_score$Score, decreasing=TRUE),]
        # Same for the GO terms (by average)
        GO_scores <- GO_scores[order(GO_scores$ave_score, decreasing=TRUE),]
    }
    # Return the results of the analysis
    if (method %in% c("randomForest", "rf")){
        return(
            list(
                GO=GO_scores,
                mapping=GO_genes,
                genes=genes_score,
                factor=f,
                method=method,
                subset=subset,
                rank.by=rank.by,
                FUN.GO=FUN.GO,
                ntree=ntree,
                mtry=mtry
                )
            )
    }
    else if (method %in% c("anova", "a")) {
        return(
            list(
                GO=GO_scores,
                mapping=GO_genes,
                genes=genes_score,
                factor=f,
                method=method,
                subset=subset,
                rank.by=rank.by,
                FUN.GO=FUN.GO
                )
            )
    }
}


mart_from_ensembl <- function(sample_gene){
    # If the gene id starts by "ENS" (most cases, except 3 handled separately
    # below)
    if (length(grep(pattern="^ENS", x=sample_gene))){
        # Extract the full prefix
        prefix <- str_extract(sample_gene, "ENS[[:upper:]]+")
        # If the ENS* prefix is in the table 
        if (prefix %in% prefix2dataset$prefix){
            # load the corresponding biomart dataset
            cat("Looks like Ensembl gene identifier.", fill=TRUE)
            cat(
                "Loading detected dataset",
                prefix2dataset[prefix2dataset$prefix == prefix,]$dataset,
                "...", fill=TRUE
                )
            return(useMart(
                biomart="ensembl",
                dataset=prefix2dataset[
                    prefix2dataset$prefix == prefix,]$dataset
            ))
        }
        # Otherwise return FALSE
        else{
            cat(
                "Did not recognise a valid Ensembl gene identifier.",
                fill=TRUE
                )
            return(FALSE)
        }
    }
    # If the gene id starts with "Y"
    else if (length(grep(pattern="^Y", x=sample_gene))) {
        # load the corresponding biomart dataset
        cat("Looks like Ensembl gene identifier.", fill=TRUE)
        cat("Loading detected dataset scerevisiae_gene_ensembl ...",
            fill=TRUE)
        return(
            useMart(
                biomart="ensembl",
                dataset="scerevisiae_gene_ensembl"
                )
            )
    }
    # If the gene id does not match any known Ensembl gene id prefix, return
    # an error and stop
    else{
        cat("Did not recognise an Ensembl gene identifier.", fill=TRUE)
        return(FALSE)
    }
}

microarray_from_probeset <- function(sample_gene, microarray2dataset.clean){
    matches <- c()
    # For each pattern manually curated as unique to a microarray
    for (pattern in microarray2dataset.clean$pattern[
        microarray2dataset.clean$unique]){
        # if the pattern matches the sample gene
        if (length(which(grep(pattern=pattern, x=sample_gene) == 1))){
            # add the pattern to a vector 
            matches <- c(matches, pattern)
        }
    }
    # if the vector length is at least 1 (patterns were designed to not
    # overlap, so the length cannot be more than 1 here)
    if (length(matches)){
        # return (dataset, microarray) to the main function to build mapping
        # tables
        return(
            microarray2dataset.clean[
                microarray2dataset.clean$pattern == matches[1],
                c("dataset","microarray")
                ]
            )
    }
    # If the sample gene was not recognised in the unique ones,
    # check whether it may be an ambiguous probeset identifier
    # For each pattern manually curated as found in multiple microarrays
    for (pattern in unique(
        microarray2dataset.clean$pattern[!microarray2dataset.clean$unique])
        ){
        # if the pattern matches the sample gene
        if (length(which(grep(pattern=pattern, x=sample_gene) == 1))){
            # add the pattern to a vecto riap
            matches <- c(matches, pattern)
        }
    }
    # if vector contains at least 1 pattern
    if (length(matches)){
        # print sample, first pattern, and list of possible microarray
        cat(
            sample_gene, "matches pattern", matches[1],
            "found in multiple microarrays:", fill=TRUE
            )
        print(
            microarray2dataset.clean[
                microarray2dataset.clean$pattern == matches[1], 
                c("dataset","microarray")
                ],
            row.names=FALSE)
        return(FALSE)
    }
    # if no known microarray pattern matches, return FALSE
    else{
        cat("Did not recognise microarray data.", fill=TRUE)
        return(FALSE)
    }
}
