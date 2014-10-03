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
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), byrow=TRUE,
                         ncol=cols, nrow=ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow(layout),
                                                 ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this
            # subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind=TRUE))
            print(plots[[i]], vp=viewport(layout.pos.row=matchidx$row,
                                          layout.pos.col=matchidx$col))
        }
    }
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
