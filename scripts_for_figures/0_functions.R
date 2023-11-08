# ----- Loading Packages -----
pkgs <- list("ggplot2", "ggpubr", "ggforce", "immunedeconv", "singscore", "dplyr", 
             "limma", "edgeR", "DESeq2", "pheatmap", "grid")
lapply(pkgs, library, character.only = TRUE)

# ----- Helper Functions -----
# convert the cibersort result from cell type by samples to samples by cell types
#
# file_path  path of the cibersort result
# @return a data frame (22 immune cell types by samples)
ReadCiber <- function(file_path) {
    ciber <- read.table(file_path, head = TRUE, sep = "\t", check.names = FALSE)
    samples <- names(ciber)[-1]
    cell_type <- ciber[, 1]
    ciber <- as.data.frame(t(ciber[, -1]))
    colnames(ciber) <- cell_type
    ciber <- cbind(samples, ciber)
    ciber
}

# calculate immune types in tumor microenvironment using cibersort method
#
# expr  mRNA expression
# method  "cibersort" or "cibersort_abs"
# @return a data frame (22 immune cell types by samples)
CalCiber <- function(expr, method) {
    set_cibersort_binary("~/Documents/Research/WangLab/Projects/IGR/script/CIBERSORT.R")
    set_cibersort_mat("~/Documents/Research/WangLab/Projects/IGR/script/LM22.txt")
    
    if (method != "cibersort" && method != "cibersort_abs") {
        stop("Wrong deconvolute method, should be cibersort or cibersort_abs")
    }
    
    cibersort <- deconvolute(expr, method)
    cibersort <- as.data.frame(cibersort)
    write.table(cibersort, "~/tmp-cibersort.txt", 
                sep = "\t", row.names = FALSE, quote = FALSE)
    ciber <- ReadCiber("~/tmp-cibersort.txt")
    file.remove("~/tmp-cibersort.txt")
    names(ciber) <- gsub(" ", "_", names(ciber))
    ciber
}

# fast read big mRNA expression files
#
# file_path  path of the mRNAseq file
# @return standard expr data frame (genes by samples)
ReadExpr <- function(file_path) {
    expr <- suppressWarnings(data.table::fread(file_path))
    expr <- as.data.frame(expr)
    rownames(expr) <- expr[, 1]
    expr[, 1] <- NULL
    expr
}

# calculate immune signatures based on Dr. Xiaosong Wang's method (singscore)
#
# expr  mRNA expression
# gmt  gmt file that has immune inflamed genesets
# @return a data frame that has different immune inflamed scores
CalImmuSig <- function(expr, 
                       gmt = "~/Documents/Research/WangLab/Projects/IGR/data_for_manuscript/Immune_inflamed_genesets.gmt") {
    read_gmt<-function(fileName,min=5){
        con=file(fileName,"r")
        i<-0
        tmp.posList<-list()
        while (length(line<-readLines(con,n=1))>0){
            #if (grepl("^#", line)){
            #print(paste(line,"skipped",sep=" "))
            #next
            #}
            tmp.line <- unlist(strsplit(line,"\t"))
            tmp.name <- tmp.line[1]
            if (length(tmp.line)>=min+2){
                i=i+1
                tmp.posList[[i]] <- tmp.line[3:length(tmp.line)]
                names(tmp.posList)[i]<-tmp.name  
                if(i %% 2000==0){
                    print(paste("Loaded ",i," features",sep=""))
                }
            }
        }
        close(con)
        return(tmp.posList)
    }
    
    genesets <- read_gmt(gmt, min = 0)
    train <- unique(unlist(genesets))
    ranked <- rankGenes(expr)
    scoredf <- simpleScore(rankData = ranked, upSet = train, 
                           centerScore = T, knownDirection = T)
    scoredf$SampleID <- rownames(scoredf)
    scoredf <- scoredf[, 1:2]
    for (i in 1:length(genesets)) {
        if (length(genesets[[i]] > 0)) {
            tmp.score <- simpleScore(rankData = ranked, 
                                     upSet = genesets[[i]], 
                                     centerScore = T, knownDirection = T)
            tmp.score$SampleID <- rownames(tmp.score)
            colnames(tmp.score) <- paste(names(genesets)[i], 
                                         colnames(tmp.score), sep="_")
            scoredf <- cbind(scoredf,tmp.score[, 1:2])
        }
    }
    totalscore <- data.frame(SampleName=rownames(scoredf),scoredf)
    totalscore
}

# calculate cell cycle expression signatures from mRNA expression
#
# expr  mRNA expression
# @return a data frame that has cell cycle signature values
CalCellCycleSig <- function(expr) {
    gene_list <- c("CENPE", "CCNA2", "CCNB2", "MCM6", "CCNF", 
                   "BUB1", "CDC20", "CDC6", "CDK1", "PLK1")
    sig_mat <- expr[na.omit(match(gene_list, rownames(expr))), ]
    val <- colMeans(sig_mat)
    result <- data.frame(icgc_specimen_id = names(val), 
                         cell_cycle_sig   = as.numeric(val))
    result
}

# convert fpkm to tpm
#
# fpkm  FPKM mRNA expression
# @return TPM mRNA expression
fpkm2tpm <- function(fpkm) {
    apply(fpkm, 2, function (x) x / sum(x)) * 1e6
}

# obtained from https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
    
    # repel.degree = number within [0, 1], which controls how much 
    #                space to allocate for repelling labels.
    ## repel.degree = 0: spread out labels over existing range of kept labels
    ## repel.degree = 1: spread out labels over the full y-axis
    
    heatmap <- pheatmap$gtable
    
    new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
    
    # keep only labels in kept.labels, replace the rest with ""
    new.label$label <- ifelse(new.label$label %in% kept.labels, 
                              new.label$label, "")
    
    # calculate evenly spaced out y-axis positions
    repelled.y <- function(d, d.select, k = repel.degree){
        # d = vector of distances for labels
        # d.select = vector of T/F for which labels are significant
        
        # recursive function to get current label positions
        # (note the unit is "npc" for all components of each distance)
        strip.npc <- function(dd){
            if(!"unit.arithmetic" %in% class(dd)) {
                return(as.numeric(dd))
            }
            
            d1 <- strip.npc(dd$arg1)
            d2 <- strip.npc(dd$arg2)
            fn <- dd$fname
            return(lazyeval::lazy_eval(paste(d1, fn, d2)))
        }
        
        full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
        selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
        
        return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                        to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                        length.out = sum(d.select)), 
                    "npc"))
    }
    new.y.positions <- repelled.y(new.label$y,
                                  d.select = new.label$label != "")
    new.flag <- segmentsGrob(x0 = new.label$x,
                             x1 = new.label$x + unit(0.15, "npc"),
                             y0 = new.label$y[new.label$label != ""],
                             y1 = new.y.positions)
    
    # shift position for selected labels
    new.label$x <- new.label$x + unit(0.2, "npc")
    new.label$y[new.label$label != ""] <- new.y.positions
    
    # add flag to heatmap
    heatmap <- gtable::gtable_add_grob(x = heatmap,
                                       grobs = new.flag,
                                       t = 4, 
                                       l = 4
    )
    
    # replace label positions in heatmap
    heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
    
    # plot result
    grid.newpage()
    grid.draw(heatmap)
    
    # return a copy of the heatmap invisibly
    invisible(heatmap)
}