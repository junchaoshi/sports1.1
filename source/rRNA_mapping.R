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
  RPM.2S <- 0
  RPM.4.5S <- 0
  RPM.5S <- 0
  RPM.5.3S <- 0
  RPM.5.8S <- 0
  RPM.12S <- 0
  RPM.16S <- 0
  RPM.17S <- 0
  RPM.18S <- 0
  RPM.25S <- 0
  RPM.26S <- 0
  RPM.28S <- 0
  RPM.45S <- 0
  rRNA.length <- unlist(strsplit(rRNA.length, ","))
  j <- 0
  for(i in 1:length(rRNA.length)){
    temp.length <- unlist(strsplit(rRNA.length[i], "="))
    if(temp.length[1] == "2S"){
      len.2S.match <- match.stat(rRNA.match, "2S", temp.length[2], total.reads)
      len.2S.unmatch <- match.stat(rRNA.unmatch, "2S", temp.length[2], total.reads)
      len.2S <- data.frame(length = len.2S.match$length, RPM = (len.2S.match$RPM + len.2S.unmatch$RPM))
      RPM.2S <- sum(len.2S$RPM)
      if(RPM.2S > 0){
        j <- j+1
        graph.2S <- graph.distr(len.2S, output.color = "#E76BF3", "Mapping to 2S rRNA")
      }
    }else if(temp.length[1] == "4.5S"){
      len.4.5S.match <- match.stat(rRNA.match, "4.5S", temp.length[2], total.reads)
      len.4.5S.unmatch <- match.stat(rRNA.unmatch, "4.5S", temp.length[2], total.reads)
      len.4.5S <- data.frame(length = len.4.5S.match$length, RPM = (len.4.5S.match$RPM + len.4.5S.unmatch$RPM))
      RPM.4.5S <- sum(len.4.5S$RPM)
      if(RPM.4.5S > 0){
        j <- j+1
        graph.4.5S <- graph.distr(len.4.5S, output.color = "#E76BF3", "Mapping to 4.5S rRNA")
      }
    }else if(temp.length[1] == "5S"){
      len.5S.match <- match.stat(rRNA.match, " 5S", temp.length[2], total.reads)
      len.5S.unmatch <- match.stat(rRNA.unmatch, " 5S", temp.length[2], total.reads)
      len.5S <- data.frame(length = len.5S.match$length, RPM = (len.5S.match$RPM + len.5S.unmatch$RPM))
      RPM.5S <- sum(len.5S$RPM)
      if(RPM.5S > 0){
        j <- j+1
        graph.5S <- graph.distr(len.5S, output.color = "#E76BF3", "Mapping to 5S rRNA")
      }
    }else if(temp.length[1] == "5.3S"){
      len.5.3S.match <- match.stat(rRNA.match, " 5.3S", temp.length[2], total.reads)
      len.5.3S.unmatch <- match.stat(rRNA.unmatch, " 5.3S", temp.length[2], total.reads)
      len.5.3S <- data.frame(length = len.5.3S.match$length, RPM = (len.5.3S.match$RPM + len.5.3S.unmatch$RPM))
      RPM.5.3S <- sum(len.5.3S$RPM)
      if(RPM.5.3S > 0){
        j <- j+1
        graph.5.3S <- graph.distr(len.5.3S, output.color = "#E76BF3", "Mapping to 5.3S rRNA")
      }
    }else if(temp.length[1] == "5.8S"){
      len.5.8S.match <- match.stat(rRNA.match, "5.8S", temp.length[2], total.reads)
      len.5.8S.unmatch <- match.stat(rRNA.unmatch, "5.8S", temp.length[2], total.reads)
      len.5.8S <- data.frame(length = len.5.8S.match$length, RPM = (len.5.8S.match$RPM + len.5.8S.unmatch$RPM))
      RPM.5.8S <- sum(len.5.8S$RPM)
      if(RPM.5.8S > 0){
        j <- j+1
        graph.5.8S <- graph.distr(len.5.8S, output.color = "#E76BF3", "Mapping to 5.8S rRNA")
      }
    }else if(temp.length[1] == "12S"){
      len.12S.match <- match.stat(rRNA.match, "12S", temp.length[2], total.reads)
      len.12S.unmatch <- match.stat(rRNA.unmatch, "12S", temp.length[2], total.reads)
      len.12S <- data.frame(length = len.12S.match$length, RPM = (len.12S.match$RPM + len.12S.unmatch$RPM))
      RPM.12S <- sum(len.12S$RPM)
      if(RPM.12S > 0){
        j <- j+1
        graph.12S <- graph.distr(len.12S, output.color = "#E76BF3", "Mapping to 12S rRNA")
      }
    }else if(temp.length[1] == "16S"){
      len.16S.match <- match.stat(rRNA.match, "16S", temp.length[2], total.reads)
      len.16S.unmatch <- match.stat(rRNA.unmatch, "16S", temp.length[2], total.reads)
      len.16S <- data.frame(length = len.16S.match$length, RPM = (len.16S.match$RPM + len.16S.unmatch$RPM))
      RPM.16S <- sum(len.16S$RPM)
      if(RPM.16S > 0){
        j <- j+1
        graph.16S <- graph.distr(len.16S, output.color = "#E76BF3", "Mapping to 16S rRNA")
      }
    }else if(temp.length[1] == "17S"){
      len.17S.match <- match.stat(rRNA.match, "17S", temp.length[2], total.reads)
      len.17S.unmatch <- match.stat(rRNA.unmatch, "17S", temp.length[2], total.reads)
      len.17S <- data.frame(length = len.17S.match$length, RPM = (len.17S.match$RPM + len.17S.unmatch$RPM))
      RPM.17S <- sum(len.17S$RPM)
      if(RPM.17S > 0){
        j <- j+1
        graph.17S <- graph.distr(len.17S, output.color = "#E76BF3", "Mapping to 17S rRNA")
      }
    }else if(temp.length[1] == "18S"){
      len.18S.match <- match.stat(rRNA.match, "18S", temp.length[2], total.reads)
      len.18S.unmatch <- match.stat(rRNA.unmatch, "18S", temp.length[2], total.reads)
      len.18S <- data.frame(length = len.18S.match$length, RPM = (len.18S.match$RPM + len.18S.unmatch$RPM))
      RPM.18S <- sum(len.18S$RPM)
      if(RPM.18S > 0){
        j <- j+1
        graph.18S <- graph.distr(len.18S, output.color = "#E76BF3", "Mapping to 18S rRNA")
      }
    }else if(temp.length[1] == "25S"){
      len.25S.match <- match.stat(rRNA.match, "25S", temp.length[2], total.reads)
      len.25S.unmatch <- match.stat(rRNA.unmatch, "25S", temp.length[2], total.reads)
      len.25S <- data.frame(length = len.25S.match$length, RPM = (len.25S.match$RPM + len.25S.unmatch$RPM))
      RPM.25S <- sum(len.25S$RPM)
      if(RPM.25S > 0){
        j <- j+1
        graph.25S <- graph.distr(len.25S, output.color = "#E76BF3", "Mapping to 25S rRNA")
      }
    }else if(temp.length[1] == "26S"){
      len.26S.match <- match.stat(rRNA.match, "26S", temp.length[2], total.reads)
      len.26S.unmatch <- match.stat(rRNA.unmatch, "26S", temp.length[2], total.reads)
      len.26S <- data.frame(length = len.26S.match$length, RPM = (len.26S.match$RPM + len.26S.unmatch$RPM))
      RPM.26S <- sum(len.26S$RPM)
      if(RPM.26S > 0){
        j <- j+1
        graph.26S <- graph.distr(len.26S, output.color = "#E76BF3", "Mapping to 26S rRNA")
      }
    }else if(temp.length[1] == "28S"){
      len.28S.match <- match.stat(rRNA.match, "28S", temp.length[2], total.reads)
      len.28S.unmatch <- match.stat(rRNA.unmatch, "28S", temp.length[2], total.reads)
      len.28S <- data.frame(length = len.28S.match$length, RPM = (len.28S.match$RPM + len.28S.unmatch$RPM))
      RPM.28S <- sum(len.28S$RPM)
      if(RPM.28S > 0){
        j <- j+1
        graph.28S <- graph.distr(len.28S, output.color = "#E76BF3", "Mapping to 28S rRNA")
      }
    }else if(temp.length[1] == "45S"){
      len.45S.match <- match.stat(rRNA.match, "45S", temp.length[2], total.reads)
      len.45S.unmatch <- match.stat(rRNA.unmatch, "45S", temp.length[2], total.reads)
      len.45S <- data.frame(length = len.45S.match$length, RPM = (len.45S.match$RPM + len.45S.unmatch$RPM))
      RPM.45S <- sum(len.45S$RPM)
      if(RPM.45S > 0){
        j <- j+1
        graph.45S <- graph.distr(len.45S, output.color = "#E76BF3", "Mapping to 45S rRNA")
      }
    }
  }

  if(j > 0){
    pdf(paste(file.address, file.name, "_result/", file.name, "_rRNA_mapping.pdf", sep=""), width = 8, height = j*2)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(j+1 ,1, heights = unit(c(1, rep(4, j)), "null"))))
    pdf.title <- paste("rRNAs Mapping Result of ", file.name, sep = "")
    vplayout <- function(x,y)
    	viewport(layout.pos.row = x, layout.pos.col = y)
  grid.text(pdf.title, vp = vplayout(1, 1), gp = gpar(fontsize = 20))
    k <- 2
    if(RPM.2S > 0){
      print(graph.2S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.4.5S > 0){
      print(graph.4.5S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.5S > 0){
      print(graph.5S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.5.3S > 0){
      print(graph.5.3S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.5.8S > 0){
      print(graph.5.8S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.12S > 0){
      print(graph.12S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.16S > 0){
      print(graph.16S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.17S > 0){
      print(graph.17S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.18S > 0){
      print(graph.18S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.25S > 0){
      print(graph.25S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.26S > 0){
      print(graph.26S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.28S > 0){
      print(graph.28S, vp = vplayout(k, 1))
      k <- k+1
    }
    if(RPM.45S > 0){
      print(graph.45S, vp = vplayout(k, 1))
      k <- k+1
    }
  }
  invisible(dev.off())
}

match.stat <- function(input, RNA.name, RNA.length, total.reads){
  rRNA.sub <- input[grep(RNA.name, input$reference, ignore.case = TRUE), 2:6]
  len <- data.frame(name = RNA.name, length = 0:RNA.length, RPM = 0)
  
  if(nrow(rRNA.sub) != 0){ 
    for(i in 1:length(rRNA.sub$reads)){
      for(j in 1:nchar(as.character(rRNA.sub$seq[i]))){
        len$RPM[rRNA.sub$start.site.1[i]+j+1] <- len$RPM[rRNA.sub$start.site.1[i]+j+1] + rRNA.sub$reads[i]/total.reads*1000000
      }
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