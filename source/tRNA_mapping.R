#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
	stop("No enough input parameters!\n",
	call.=FALSE)
}

library("ggplot2")
library("data.table")
library("grid")
library(stringr)

bar.chart <- function(input, name){
  p <- ggplot(input, aes(x = length, y = RPM)) + 
    geom_area(fill = "#E76BF3") + 
    theme(axis.line = element_line() , 
          panel.grid.minor = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.background = element_blank() ,
          legend.title = element_blank() , 
          axis.text = element_text(size = 15, color = "black") ,
          axis.title = element_text(size = 15) ,
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Length", y = "RPM", title = name)
  return(p)
}

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

tRNA.mapping <- read.table(args[1], sep = "\t")
tRNA.number <- nrow(tRNA.mapping)
pdf(args[2], width = 8, height = 3)
for (i in 1:tRNA.number) {
  rpm <- strsplit(as.character(tRNA.mapping[i,2]), ",")
  title.name <- as.character(tRNA.mapping[i,1])
  len <- length(rpm[[1]])
  length.data <- data.frame(length = c(1:len), RPM = as.numeric(rpm[[1]]))
  pic <- bar.chart(length.data, title.name)
  print(pic)
}
invisible(dev.off())
