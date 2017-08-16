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
  len.min    <- min(length.dis$length)
  len.max    <- max(length.dis$length)
  len        <- len.min : len.max
  sum.clean  <- sum(length.dis[which(length.dis$name == "Clean_Reads"), 3])
  length.dis <- data.frame(length.dis[ ,1:2], reads = length.dis$reads * 1000000 / sum.clean)
  
  dis.clean  <- length.dis[which(length.dis$name == "Clean_Reads"), 2:3]
  stack.all  <- c()

  dis.total.rRNA <- rbind(length.dis[grep("rRNAdb-rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                          length.dis[grep("ensembl-Mt_rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                          length.dis[grep("ensembl-rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                          length.dis[grep("Rfam-rRNA", length.dis$name, ignore.case = TRUE), 2:3]
                         )

  dis.other.rRNA <- data.frame(name = "other.rRNA", length.combine(dis.total.rRNA, len))
  
  rRNA.length <- unlist(strsplit(rRNA.length, ","))
  k <- 0
  for(i in 1:length(rRNA.length)){
    temp.length <- unlist(strsplit(rRNA.length[i], "="))
    if(temp.length[1] == "2S"){
      k <- k+1
      dis.2S.rRNA <- length.dis[grep("-2S", length.dis$name), 2:3]
      dis.2S.rRNA <- data.frame(name = "2S.rRNA", length.combine(dis.2S.rRNA, len))
      sum.2S.rRNA <- sum(dis.2S.rRNA$reads)
      if(sum.2S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.2S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.2S.rRNA[3])
      }
    }else if(temp.length[1] == "4.5S"){
      k <- k+1
      dis.4.5S.rRNA <- length.dis[grep("-4.5S", length.dis$name), 2:3]
      dis.4.5S.rRNA <- data.frame(name = "4.5S.rRNA", length.combine(dis.4.5S.rRNA, len))
      sum.4.5S.rRNA <- sum(dis.4.5S.rRNA$reads)
      if(sum.4.5S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.4.5S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.4.5S.rRNA[3])
      }
    }else if(temp.length[1] == "5S"){
      k <- k+1
      dis.5S.rRNA <- length.dis[grep("-5S", length.dis$name), 2:3]
      dis.5S.rRNA <- data.frame(name = "5S.rRNA", length.combine(dis.5S.rRNA, len))
      sum.5S.rRNA <- sum(dis.5S.rRNA$reads)
      if(sum.5S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.5S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.5S.rRNA[3])
      }
    }else if(temp.length[1] == "5.3S"){
      k <- k+1
      dis.5.3S.rRNA <- length.dis[grep("-5.3S", length.dis$name), 2:3]
      dis.5.3S.rRNA <- data.frame(name = "5.3S.rRNA", length.combine(dis.5.3S.rRNA, len))
      sum.5.3S.rRNA <- sum(dis.5.3S.rRNA$reads)
      if(sum.5.3S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.5.3S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.5.3S.rRNA[3])
      }
    }else if(temp.length[1] == "5.8S"){
      k <- k+1
      dis.5.8S.rRNA  <- length.dis[grep("-5.8S", length.dis$name), 2:3]
      dis.5.8S.rRNA  <- data.frame(name = "5.8S.rRNA", length.combine(dis.5.8S.rRNA, len))
      sum.5.8S.rRNA <- sum(dis.5.8S.rRNA$reads)
      if(sum.5.8S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.5.8S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.5.8S.rRNA[3])
      }
    }else if(temp.length[1] == "12S"){
      k <- k+1
      dis.12S.rRNA  <- length.dis[grep("-12S", length.dis$name), 2:3]
      dis.12S.rRNA  <- data.frame(name = "12S.rRNA", length.combine(dis.12S.rRNA, len))
      sum.12S.rRNA <- sum(dis.12S.rRNA$reads)
      if(sum.12S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.12S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.12S.rRNA[3])
      }
    }else if(temp.length[1] == "16S"){
      k <- k+1
      dis.16S.rRNA  <- length.dis[grep("-16S", length.dis$name), 2:3]
      dis.16S.rRNA  <- data.frame(name = "16S.rRNA", length.combine(dis.16S.rRNA, len))
      sum.16S.rRNA <- sum(dis.16S.rRNA$reads)
      if(sum.16S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.16S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.16S.rRNA[3])
      }
    }else if(temp.length[1] == "17S"){
      k <- k+1
      dis.17S.rRNA  <- length.dis[grep("-17S", length.dis$name), 2:3]
      dis.17S.rRNA  <- data.frame(name = "17S.rRNA", length.combine(dis.17S.rRNA, len))
      sum.17S.rRNA <- sum(dis.17S.rRNA$reads)
      if(sum.17S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.17S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.17S.rRNA[3])
      }
    }else if(temp.length[1] == "18S"){
      k <- k+1
      dis.18S.rRNA  <- length.dis[grep("-18S", length.dis$name), 2:3]
      dis.18S.rRNA  <- data.frame(name = "18S.rRNA", length.combine(dis.18S.rRNA, len))
      sum.18S.rRNA <- sum(dis.18S.rRNA$reads)
      if(sum.18S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.18S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.18S.rRNA[3])
      }
    }else if(temp.length[1] == "25S"){
      k <- k+1
      dis.25S.rRNA  <- length.dis[grep("-25S", length.dis$name), 2:3]
      dis.25S.rRNA  <- data.frame(name = "25S.rRNA", length.combine(dis.25S.rRNA, len))
      sum.25S.rRNA <- sum(dis.25S.rRNA$reads)
      if(sum.25S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.25S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.25S.rRNA[3])
      }
    }else if(temp.length[1] == "26S"){
      k <- k+1
      dis.26S.rRNA  <- length.dis[grep("-26S", length.dis$name), 2:3]
      dis.26S.rRNA  <- data.frame(name = "26S.rRNA", length.combine(dis.26S.rRNA, len))
      sum.26S.rRNA <- sum(dis.26S.rRNA$reads)
      if(sum.26S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.26S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.26S.rRNA[3])
      }
    }else if(temp.length[1] == "28S"){
      k <- k+1
      dis.28S.rRNA  <- length.dis[grep("-28S", length.dis$name), 2:3]
      dis.28S.rRNA  <- data.frame(name = "28S.rRNA", length.combine(dis.28S.rRNA, len))
      sum.28S.rRNA <- sum(dis.28S.rRNA$reads)
      if(sum.28S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.28S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.28S.rRNA[3])
      }
    }else if(temp.length[1] == "45S"){
      k <- k+1
      dis.45S.rRNA  <- length.dis[grep("-45S", length.dis$name), 2:3]
      dis.45S.rRNA  <- data.frame(name = "45S.rRNA", length.combine(dis.45S.rRNA, len))
      sum.45S.rRNA <- sum(dis.45S.rRNA$reads)
      if(sum.45S.rRNA > 0){
        stack.all <- rbind(stack.all, dis.45S.rRNA)
        dis.other.rRNA <- data.frame(dis.other.rRNA[1:2], dis.other.rRNA[3] - dis.45S.rRNA[3])
      }
    }
  }
  stack.all <- rbind(stack.all, dis.other.rRNA)
    
  stack.sep <- ggplot(stack.all, aes(x = length, y = reads)) + 
                    geom_bar(stat = "identity") + 
                    facet_grid(. ~ name) +
					labs(x = "length", y = "RPM", title = "") + 
					theme(axis.line = element_line()
					, panel.grid.minor = element_blank() 
					, panel.grid.major = element_blank()
                    , panel.background = element_blank()
					, panel.border = element_rect(fill = "transparent")					
					)
					
  
  stack.pic <- ggplot(stack.all, aes(x = length, y = reads, fill = name)) + 
                                     geom_bar(stat = "identity") +
				     theme(axis.line = element_line()
				   , panel.grid.minor = element_blank() 
				   , panel.grid.major = element_blank()
                                   , panel.background = element_blank()
                                   , legend.title = element_blank()
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
  dev.off()
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


