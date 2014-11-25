# Code borrowed from the web to create a lattice of ggplot2 plots.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # source:
    # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(
            seq(1, cols * ceiling(numPlots/cols)),
            byrow=TRUE,
            ncol=cols, nrow=ceiling(numPlots/cols)
            )
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(
            viewport(
                layout=grid.layout(
                    nrow(layout),
                    ncol(layout)
                    )
                )
            )
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this
            # subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind=TRUE))
            print(
                plots[[i]], vp=viewport(
                    layout.pos.row=matchidx$row,
                    layout.pos.col=matchidx$col
                    )
                )
        }
    }
    return(numPlots)
}

# Splits a string of characters into multiple substrings, each less than 
# a given number of characters. New line characters cannot be inserted within
# words. Words are defined as surrounded by space characters only.
string_Lsplit <- function (string, line.length){
    # Get the (ordered) list of words
    words <- strsplit(x=string, split=" ", )[[1]]
    # Rebuild the original string, while inserting a newline everytime
    # the limit is reached
    # Start with empty title
    newString <- words[1]
    # Count of characters since latest newline
    nc <- nchar(words[1])
    for (word in words[2:length(words)]){
        if (nc + nchar(word) > line.length){
            newString <- paste(newString, word, sep="\n")
            nc <- nchar(word)
        }
        else{
            newString <- paste(newString, word, sep=" ")
            nc <- nc + nchar(word) + 1 # for space character !
        }
    }
    return(newString)
}

# Funtion to build the prefix2dataset table
prefix2dataset.build <- function(){
    # This function requires the libraries biomaRt and RCurl to be preloaded
    # load the RCurl library (used in a loop later below)
    #library("RCurl", quietly=TRUE)
    curlHandle <- getCurlHandle()
    # load the biomaRt package
    #library(biomaRt, quietly=TRUE)
    # Connect to the Ensembl mart
    mart <- useMart(biomart="ensembl")
    # Save the list of datasets available 
    mart.datasets <- listDatasets(mart=mart)
    # Extract the "dataset" column which is the value to access the mart
    mart.dataset <- as.character(mart.datasets$dataset)
    # Extract the species name from the "description" column
    mart.species <- as.character(
        sapply(
            sapply(
                X=mart.datasets$description,
                FUN=strsplit, " genes"
                ),
            "[[", 1
            )
        )
    # For each dataset, fetch a random ensembl_gene_id as an example
    mart.sample <- as.character(
        sapply(
            X=mart.dataset,
            FUN=sampleEnsemblGeneId,
            curl=curlHandle
            )
    )
    # For each sample ensembl_gene_id, identify the prefix defined as the
    # letters starting the ensembl_gene_id
    mart.prefix <- as.character(
        sapply(
            sapply(
                X=mart.sample,
                FUN=strsplit, "[0-9]+"
                ),
            "[[", 1
            )
    )
    # Collate the data in the table
    p2d.table <- data.frame(
        row.names=NULL,
        dataset=mart.dataset,
        species=mart.species,
        prefix=mart.prefix,
        sample=mart.sample,
        stringsAsFactors=FALSE)
    # Sort species names alphabetically for ease of human reading
    p2d.table <- p2d.table[order(p2d.table$species),]
    return(p2d.table)
}

# Function called in prefix2dataset.build sapply statement to fetch
# a sample Ensemblgene identifier given a 
sampleEnsemblGeneId <- function(dataset, curl=getCurlHandle()){
    # User message
    cat("Fetching data for dataset:", dataset, fill=TRUE)
    # connect to the specific mart
    mart.loop <- useMart(biomart="ensembl", dataset=dataset)
    # query the first (automatically non-empty) ensembl_gene_id
    ensembl_gene_id <- getBM(
        attributes="ensembl_gene_id",
        mart=mart.loop,
        curl=curl
        )[1,"ensembl_gene_id"]
    # return the above ensembl_gene_id
    return(ensembl_gene_id)
}

# Funtion to build the microarray2dataset table
microarray2dataset.build <- function(){
    # This function requires the libraries biomaRt and RCurl to be preloaded
    # load the RCurl library (used in a loop later below)
    #library("RCurl", quietly=TRUE)
    curlHandle <- getCurlHandle()
    # load the biomaRt package
    #library(biomaRt, quietly=TRUE)
    # Connect to the Ensembl mart
    mart <- useMart(biomart="ensembl")
    # Save the list of datasets available 
    mart.datasets <- listDatasets(mart=mart)
    # Extract the "dataset" column which is the value to access the mart
    mart.dataset <- as.character(mart.datasets$dataset)
    # Extract the species name from the "description" column
    mart.species <- sapply(
        sapply(
            X=mart.datasets$description,
            FUN=strsplit, " genes"
        ),
        "[[", 1)
    getBM.results <- data.frame(
        dataset=NA, microarray=NA, sample=NA,
        stringsAsFactors=FALSE
        )
    # Count how many species processed
    species_index = 0
    # For each dataset (= species)
    for (dataset.loop in mart.dataset){
        species_index = species_index + 1
        # User message
        cat("Fetching data from dataset: ", dataset.loop,
            " (", species_index, ")",
            sep="", fill=TRUE
            )
        mart.loop <- useMart(biomart="ensembl", dataset=dataset.loop)
        # Query all column header for this dataset
        attributes.loop <- listAttributes(mart=mart.loop, page="feature_page")
        # list all microarray column headers for this dataset
        microarray.headers <- attributes.loop$name[
            grep(
                pattern="probe",
                x=attributes.loop$description,
                ignore.case="TRUE"
                )
            ]
        # For each microarray dataset
        for (microarray.header in microarray.headers){
            # User message
            cat("Fetching data for microarray:", microarray.header, fill=TRUE)
            # Query the first (automatically non-empty) ensembl_gene_id
            probe.set = getBM(
                attributes=microarray.header,
                mart=mart.loop,
                curl=curlHandle
                )[1,microarray.header]
            getBM.results <- rbind(
                getBM.results,
                c(
                    dataset=dataset.loop,
                    microarray=microarray.header,
                    sample=probe.set
                    )
                )
        }
    }
    # Remove the initial blank row
    getBM.results = getBM.results[!is.na(getBM.results$dataset),]
    # Merge the species name with the information collected in the loop above
    m2d.table <- merge(
        y=data.frame(
            dataset=mart.dataset,
            species=mart.species,
            stringsAsFactors=FALSE
            ),
        by="dataset",
        all=TRUE)
    
    # Insert columns for pattern and uniqueness with empty data
    m2d.table <- cbind(m2d.table, pattern=NA, unique=NA)
    
    # for each pattern
    for (pattern in patterns){
        # list the indices of microarray(s) that match
        match.indices  <- which(
            sapply(
                X=m2d.table$sample,
                FUN=grep,
                pattern=pattern
                ) == 1
            )
        # if one unique microarray matches
        if (nrow(m2d.table[match.indices,]) == 1){
            ## if the microarray already has a pattern
            if (!is.na(m2d.table[match.indices,"pattern"])){
                ### stop(microarray matches multiple patterns)
                stop(
                    "microarray already matched a pattern:",
                    m2d.table[match.indices,][i,]
                    )
            }
            else{
                ### set the uniqueness to TRUE
                m2d.table[match.indices,"unique"] <- TRUE
                ### set the pattern
                m2d.table[match.indices,"pattern"] <- pattern
            }
            ### set the uniqueness to TRUE
            ### set the pattern
            # if one unique microarray matches
        } else if (nrow(m2d.table[match.indices,]) > 1) {
            ## for each matching microarray
            for (i in 1:nrow(m2d.table[match.indices,])){
                ### if the microarray already has a pattern
                if (!is.na(m2d.table[match.indices,][i, "pattern"])){
                    #### stop(microarray matches multiple patterns)
                    stop(
                        "microarray already matched a pattern:",
                        m2d.table[match.indices,][i,]
                        )
                }
                else{
                    #### set the uniqueness to TRUE
                    m2d.table[match.indices,][i, "unique"] <- FALSE
                    #### set the pattern
                    m2d.table[match.indices,][i, "pattern"] <- pattern
                }
            }
            
            # if the pattern did not match any microarray, it is useless
        } else {
            warning(
                "Pattern ", pattern,
                " did not match any microarray. Edit or remove.")
        }
    }
    # Organise by species
    m2d.table <- m2d.table[order(m2d.table$species),]
}

# Patterns of microarray probe(set)s used to build the microarray2dataset
# table
patterns <- c(
    "^aa[[:digit:]]+_[a-z]_at$",
    "^A_[[:digit:]]{2}_P[[:digit:]]+$",
    "^AB[[:digit:]]+_at$",
    "^Bt\\..*_at$",
    "^Cf.*_at$",
    "^Dr\\.[[:digit:]]+.*_at$",
    "^GE[[:digit:]]+$",
    "^Gga\\..*_at$",
    "^Hs2\\.[[:digit:]]+\\..*_at$",
    "^Hs\\.[[:digit:]]+\\..*_at$",
    "^ILMN_[[:digit:]]+$",
    "^LOL[[:digit:]]+$",
    "^MKG\\..*_at$",
    "^M[[:digit:]]+_[a-z]_at$",
    "^M[[:digit:]]+_at$",
    "^Mmu\\.[[:digit:]]+.*_at$",
    "^Msa\\.[[:digit:]]+.*_at$",
    "^NM_[[:digit:]]+.*$",
    "^OaE_[[:digit:]]+$",
    "^OTTDART[[:digit:]]+_.*$",
    "^PH_hs_[[:digit:]]+$",
    "^PH_mM_[[:digit:]]+$",
    "^PH_rn_[[:digit:]]+$",
    "^rc_AA[[:digit:]]+_at$",
    "^rc_AI[[:digit:]]+_at$",
    "^S[[:digit:]]+_[a-z]_at$",    
    "^Ssc\\.[[:digit:]]+.*_at$",
    "^Str\\.[[:digit:]]+.*_at$",
    "^TC[[:digit:]]+\\.hg$",
    "^[[:digit:]]+_at$",
    "^[[:digit:]]+_[a-z]_at$",
    "^[[:digit:]]+$"
    )