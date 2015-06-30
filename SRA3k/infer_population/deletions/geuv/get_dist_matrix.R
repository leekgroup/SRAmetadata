indelsCount <- function(junctionFiles, 
                         sampleNames=names(junctionFiles), 
                         strandSpecific = FALSE, illuminaStranded=FALSE,
                         minCount = 1, maxCores=NULL) {
    
    require(GenomicRanges, quietly=TRUE)
    require(parallel, quietly=TRUE)
    
    if(is.null(maxCores)) {
        maxCores=1
    }
    
    names(junctionFiles) = sampleNames
    cat("Reading in data.\n")
    if(all(is.character(junctionFiles))) {
        theData = mclapply(junctionFiles, function(x) {
            y = read.delim(x, skip = 1, header=FALSE, 
                           col.names = c("chr", "start","end", "sequence", "count"), 
                           colClasses = c("character", "integer", "integer", "character","integer"))
            y = y[y$count >= minCount,] # filter based on min number
            gr = GRanges(y$chr, IRanges(y$start, y$end), 
                         sequence = y$sequence, count = y$count)
            return(gr)
        }, mc.cores=maxCores)
    } else {
        theData = junctionFiles
        stopifnot(all(sapply(theData, class)=="GRanges"))
    }
    cat("Creating master table of junctions.\n")
    
    ## turn into GRangesList
    ### THIS STEP IS SLOW...
    grList = GRangesList(theData)
    # each dataset should be checked
    if(illuminaStranded & strandSpecific) {
        grList = GRangesList(mclapply(grList, function(x) {
            strand(x) = ifelse(strand(x)=="+", "-","+")
            return(x)
        },mc.cores=maxCores))
    }
    
    ## get into GRanges object of unique junctions
    fullGR = unlist(grList)
    if(!strandSpecific) strand(fullGR) = "*"
    
    fullGR = fullGR[!duplicated(fullGR)] # or unique(fullGR)
    fullGR = sort(fullGR)
    fullGR$count = NULL
    
    cat(paste("There are", length(fullGR), "total junctions.\n"))
    
    cat("Populating count matrix.\n")
    
    jNames = paste0(as.character(seqnames(fullGR)),":",
                    start(fullGR),"-",end(fullGR),"(",as.character(strand(fullGR)),")")
    
    ## match GRanges
    options(warn=-1)
    mList = mclapply(grList, match, fullGR, 
                     ignore.strand = !strandSpecific, mc.cores=maxCores)
    options(warn=0)
    
    countList = mList # initiate 
    M = length(jNames)
    
    ## fill in matrix
    for(i in seq(along=grList)) {
        if(i %% 25 == 0) cat(".")
        cc = rep(0,M)
        cc[mList[[i]]] = theData[[i]]$count
        countList[[i]] = Rle(cc)
    }
    countDF = DataFrame(countList, row.names=jNames)
    
    names(fullGR) = jNames
    ## return matrix and GRanges object
    out = list(countDF = countDF, anno = fullGR)
    return(out)
}

#####################@@@MAIN@@@###########################

library('magrittr')
library('stringr')
library('Cairo')

# Set paths to files

path  <- "/dcl01/leek/data/geuvadis/introns_and_indels"
files <- list.files(path,
            pattern="^deletions\\..*\\.bed$", full.names=TRUE)

f         <- str_split(files, "/")
file_list <- lapply(f, `[`, length(f[[1]])) %>% unlist()


# Get population and ids from file names
file_list <- str_split(file_list, "_")
id_1      <- file_list %>% lapply(`[`, 1) %>% unlist() %>% 
    str_split("\\.") %>% lapply(`[`, 2) %>% unlist()
center    <- file_list %>% lapply(`[`, 4) %>% unlist()
id_2      <- file_list %>% lapply(`[`, 5) %>% unlist() %>% 
    str_split("\\.") %>% lapply(`[`, 1) %>% unlist()
pop       <- file_list %>% lapply(`[`, 3) %>% unlist()
ids       <- paste(pop, id_1, center, id_2, sep = "_")


# Get count matrix
count_m <- indelsCount(junctionFiles = files,
               sampleNames = ids,
               maxCores = 10)

save(count_m, file = "deletions_countmatrix_geuvadis.rda")



