subEset <- function(eSet, subset=list()){
    # subset should be a list of factor names with the associated values
    # to retain for that factor
    if (!class(subset) == "list"){
        stop("\"subset=\" argument should be a list.")
    }
    if (length(subset) > 0){
        for (f_filter in names(subset)){
            # Check that the name of the list item is a column name in
            # the phenodata
            if (!f_filter %in% colnames(pData(eSet))){
                stop(f_filter, " is not a valid column in pData(eSet).")
            }
            if (length(subset[[f_filter]]) == 0){
                stop("Fo value provided for filter ", f_filter)
            }
            for (v_filter in subset[[f_filter]]){
                # Check that the value is an existing value 
                # in the phenodata
                if(!v_filter %in% pData(eSet)[,f_filter]){
                    stop(v_filter,
                         " is not a valid value of eSet$", f_filter)
                }
            }
            # at this stage, column and values exist
            # Subset the eSet to the
            eSet <- eSet[,pData(eSet)[,f_filter] %in% subset[[f_filter]]]
            # Update the factor levels
            if ("factor" %in% class(pData(eSet)[,f_filter])){
                pData(eSet)[,f_filter] = factor(pData(eSet)[,f_filter])
            }
        }
    }
    else{
        message("Empty list of filters given. Returning the original eSet")
    }
    return(eSet)
}
