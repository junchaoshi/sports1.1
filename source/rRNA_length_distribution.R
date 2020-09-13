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

rRNA.stack <- function(file.address, file.name, rRNA.length){
 
  length.dis <- read.table(paste(file.address
				 , file.name
				 , "_result/"
				 , file.name
				 , "_length_distribution.txt", sep = ""), 
                           header = TRUE, 
                           sep = "\t", 
                           col.names = c("name", "length", "reads")
  )
  length.dis <- length.dis[!grepl("antisense", length.dis$name, ignore.case = TRUE), ]
  len.min    <- min(length.dis$length)
  len.max    <- max(length.dis$length)
  len        <- len.min : len.max
  sum.clean  <- sum(length.dis[which(length.dis$name == "Clean_Reads"), 3])
  length.dis <- data.frame(length.dis[ ,1:2], reads = length.dis$reads * 1000000 / sum.clean)
  
  dis.clean  <- length.dis[which(length.dis$name == "Clean_Reads"), 2:3]
  stack.all  <- c()
						 
  dis.total.rRNA	<- length.dis[grep("-rRNA|_rRNA", length.dis$name, ignore.case = TRUE), ]
  dis.total.rRNA	<- dis.total.rRNA[!grepl("antisense|S-rRNA|other-rRNA", dis.total.rRNA$name, ignore.case = TRUE), 2:3]
						 
  dis.other.rRNA <- data.frame(name = "other.rRNA", length.combine(dis.total.rRNA, len))
  
  rRNA.length <- unlist(strsplit(rRNA.length, ","))
  k <- 0
  for(i in 1:length(rRNA.length)){
    temp.length <- unlist(strsplit(rRNA.length[i], "="))
    if (!grepl("RNY|other", temp.length[1])){
      dis.rRNA <- length.dis[grep(paste("-", temp.length[1], sep = ""), length.dis$name), 2:3]
      dis.rRNA <- data.frame(name = paste(temp.length[1], ".rRNA", sep = ""), length.combine(dis.rRNA, len))
      sum.rRNA <- sum(dis.rRNA$reads)
      if(sum.rRNA > 0){
        k <- k+1
        stack.all <- rbind(stack.all, dis.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.rRNA[3])
      }
    }
  }
  
  stack.all <- rbind(stack.all, dis.other.rRNA)
    
  stack.sep <- ggplot(stack.all, aes(x = length, y = reads, fill = name)) + 
                    geom_bar(stat = "identity") + 
                    facet_grid(. ~ name) +
					labs(x = "length", y = "RPM", title = "") + 
					theme(axis.line = element_line()
					, panel.grid.minor = element_blank() 
					, panel.grid.major = element_blank()
                    , panel.background = element_rect(fill = "transparent", colour = NA)
					, plot.background = element_rect(fill = "transparent", colour = NA)
					, panel.border = element_rect(fill = "transparent")					
					)
					
  
  stack.pic <- ggplot(stack.all, aes(x = length, y = reads, fill = name)) + 
                      geom_bar(stat = "identity") +
				     theme(axis.line = element_line()
				   , panel.grid.minor = element_blank() 
				   , panel.grid.major = element_blank()
				   , panel.background = element_rect(fill = "transparent", colour = NA)
				   , plot.background = element_rect(fill = "transparent", colour = NA)
                   , legend.title = element_blank()
				   , legend.background = element_rect(fill = "transparent", colour = NA)
				   , legend.box.background = element_rect(fill = "transparent", colour = NA)
                   , legend.position = "top"
				   , axis.text = element_text(size = 15, color = "black")
				   , axis.title = element_text(size = 15)
			   	     ) +
	             labs(x = "length", y = "RPM", title = "") 

  pdf(paste(file.address, file.name, "_result/", file.name, "_rsRNA_distribution.pdf", sep=""), width = 8+0.5 * k, height = 12)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 1, heights = unit(c(1, 4, 4, 4), "null"))))
  pdf.title <- paste("rRNAs Length Distribution of ", file.name, sep = "")
  vplayout <- function(x,y)
	viewport(layout.pos.row = x, layout.pos.col = y)
  grid.text(pdf.title, vp = vplayout(1, 1), gp = gpar(fontsize = 20))
  print(stack.pic, vp = vplayout(2:3, 1))
  print(stack.sep, vp = vplayout(4, 1))
  invisible(dev.off())
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

rRNA.stack(file.address = args[1], file.name = args[2], rRNA.length = args[3])
