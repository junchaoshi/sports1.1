#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
	stop("No enough input parameters!\n",
	call.=FALSE)
}

library("ggplot2")
library("data.table")
library("grid")
library(stringr)

ref.distr <- function(file.address, file.name, rRNA.length){ 
  length.dis <- read.table(paste(file.address
				 , file.name
				 , "_result/"
				 , file.name
				 , "_length_distribution.txt", sep = ""), 
                           header = TRUE, 
                           sep = "\t", 
                           col.names = c("name", "length", "reads")
  )
  total.reads  <- sum(length.dis[which(length.dis$name == "Clean_Reads"), 3])
  
  rRNA.match <- read.table(paste(file.address
                                 , file.name
				 , "_processed/"
				 , file.name
				 , "_output_rRNA_match_genome", sep = ""), 
			   header = FALSE, 
			   sep = "\t", 
                           quote = "",
			   col.names = c("name"
					 , "reads"
					 , "strand"
					 , "reference"
					 , "start.site.1"
					 , "seq"
					 , "quality"
					 , "map.time.1"
					 , "mismatch"
					 )
			  )
									  
  rRNA.unmatch <- read.table(paste(file.address
                                 , file.name
				 , "_processed/"
				 , file.name
				 , "_output_rRNA_unmatch_genome", sep = ""), 
			   header = FALSE, 
			   sep = "\t", 
                           quote = "",
			   col.names = c("name"
					 , "reads"
					 , "strand"
					 , "reference"
					 , "start.site.1"
					 , "seq"
					 , "quality"
					 , "map.time.1"
					 , "mismatch"
					 )
			  )
  
  rRNA.length <- unlist(strsplit(rRNA.length, ","))
  graph <- list()
  for(i in 1:length(rRNA.length)){
    temp.length <- unlist(strsplit(rRNA.length[i], "="))
    len.match <- match.stat(rRNA.match, paste(" ", temp.length[1], sep = ""), temp.length[2], total.reads)
    len.unmatch <- match.stat(rRNA.unmatch, paste(" ", temp.length[1], sep = ""), temp.length[2], total.reads)
    len <- data.frame(length = len.match$length, RPM = (len.match$RPM + len.unmatch$RPM))
    RPM <- sum(len$RPM)
    if(RPM > 0){
  	  if (!grepl("RNY", temp.length[1])){
  		  graph[[temp.length[1]]] <- graph.distr(len, output.color = "#E76BF3", paste("Mapping to ", temp.length[1], " rRNA", sep = ""))
  	  }
    }
  }
  
  if(length(graph) > 0){
    pdf(paste(file.address, file.name, "_result/", file.name, "_rRNA_mapping.pdf", sep=""), width = 8, height = length(graph)*2)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(length(graph)+1 ,1, heights = unit(c(1, rep(4, length(graph))), "null"))))
    pdf.title <- paste("rRNAs Mapping Result of ", file.name, sep = "")
    vplayout <- function(x,y)
      viewport(layout.pos.row = x, layout.pos.col = y)
    grid.text(pdf.title, vp = vplayout(1, 1), gp = gpar(fontsize = 20))
    for (i in 1:length(graph)){
      print(graph[[i]], vp = vplayout(i+1, 1))
    }
    invisible(dev.off())
  }
}

match.stat <- function(input, RNA.name, RNA.length, total.reads){
  rRNA.sub <- input[grep(RNA.name, input$reference, ignore.case = TRUE), 2:6]
  len <- data.frame(name = RNA.name, length = 0:RNA.length, RPM = 0)
  
  if(nrow(rRNA.sub) != 0){ 
    for(i in 1:length(rRNA.sub$reads)){
      start <- rRNA.sub$start.site.1[i]+2
      end <- rRNA.sub$start.site.1[i]+nchar(as.character(rRNA.sub$seq[i]))+1
      len$RPM[start:end] <- len$RPM[start:end] + rRNA.sub$reads[i]/total.reads*10^6
    }
  }
  return(len)
}

graph.distr <- function(input.file, output.color = "black", output.title){
  graph.out <- ggplot(input.file, aes(x = length, y = RPM), fill = name) + 
	            geom_area(fill = output.color) + 
	            theme(axis.line = element_line()
                    , panel.grid.minor = element_blank() 
			        , panel.grid.major = element_blank()
                    , panel.background = element_blank()
                    , legend.title = element_blank()
			        , axis.text = element_text(size = 15, color = "black")
                    , axis.title = element_text(size = 15)
			          ) +
			    labs(x = "", y = "RPM", title = output.title)
  return(graph.out)			
}

length.combine <- function(distr, len){
  if(is.data.frame(distr) && nrow(distr) == 0){
    distr <- data.frame(length = len, reads = 0)
  }
  tmp <- data.table(distr)
  tmp <- as.data.frame(tmp[, lapply(.SD, sum), by=list(length)])
  tmp <- tmp[order(as.numeric(tmp$length)), ]
  if(length(len[-c(match(distr[ ,1], len))]) != 0){
    tmp <- rbind(tmp, data.frame(length = len[-c(match(distr[ ,1], len))], reads = 0))
  }
  distr <- tmp[order(as.numeric(tmp$length)), ]
}

ref.distr(file.address = args[1], file.name = args[2], rRNA.length = args[3])
