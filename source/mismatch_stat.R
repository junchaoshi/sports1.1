#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
  stop("No enough input parameters!\n",
       call.=FALSE)
}

###calculate p value
cal.p <- function(ref, mut, err){
  pbinom(round(ref), sum(round(ref), round(mut)), 1 - err, lower.tail = TRUE)
}

###judge significance
judge.sig <- function(x){
  if(x < 0.05){
    "TRUE"
  } else {
    "FALSE"
  }
}

seq.err <- as.numeric(args[2])

mismatch <- read.table(args[1], 
                       header = FALSE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE)
colnames(mismatch) <- c("ID", 
                        "seqs", 
                        "mut.pos", 
                        "ref.base", 
                        "mut.base", 
                        "ref.reads", 
                        "mut.reads", 
                        "adj.mut.reads")

H1.p <- apply(mismatch[, c('ref.reads', 'mut.reads')], 1, function(x) cal.p(x['ref.reads'], x['mut.reads'], seq.err))
H1.padj <- p.adjust(H1.p, method = "fdr") 

H2.p <- apply(mismatch[, c('ref.reads', 'adj.mut.reads')], 1, function(x) cal.p(x['ref.reads'], x['adj.mut.reads'], seq.err))
H2.padj <- p.adjust(H2.p, method = "fdr")

sig <- apply(data.frame(H2.padj), 1, function(x) judge.sig(x))

mismatch <- data.frame(mismatch, H1.p, H1.padj, H2.p, H2.padj, sig)
mismatch <- mismatch[order(mismatch[,1], mismatch[,3]),]

write.table(mismatch, args[1], 
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)