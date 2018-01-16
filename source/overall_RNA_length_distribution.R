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

length.stack <- function(file.address, file.name){
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
  dis.miRNA  <- length.dis[grep("mir", length.dis$name, ignore.case = TRUE), 2:3]
  dis.tRNA   <- rbind(length.dis[grep("tRNA_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                      length.dis[grep("tRNA_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3]
                      )

  dis.tRNA.5.end    <- length.dis[grep("tRNA_5_end", length.dis$name, ignore.case = TRUE), 2:3]
  dis.tRNA.3.end    <- length.dis[grep("tRNA_3_end", length.dis$name, ignore.case = TRUE), 2:3]
  dis.tRNA.CCA.end  <- length.dis[grep("tRNA_CCA_end", length.dis$name, ignore.case = TRUE), 2:3]

  dis.rRNA   <- rbind(length.dis[grep("rRNAdb-rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                      length.dis[grep("ensembl-Mt_rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                      length.dis[grep("ensembl-rRNA", length.dis$name, ignore.case = TRUE), 2:3],
                      length.dis[grep("Rfam-rRNA", length.dis$name, ignore.case = TRUE), 2:3]
                     )
  dis.piRNA  <- length.dis[grep("piRNA", length.dis$name, ignore.case = TRUE), 2:3]
  dis.unanno.match <- length.dis[grep("Unannotated_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3]
  dis.unanno.unmatch <- length.dis[grep("Unannotated_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3]
  
  dis.clean  <- length.combine(dis.clean, len)
  dis.miRNA  <- data.frame(name = "miRNA", length.combine(dis.miRNA, len))
  dis.tRNA   <- data.frame(name = "tsRNA", length.combine(dis.tRNA, len))

  dis.tRNA.5.end   <- data.frame(name = "tsRNA-5'end", length.combine(dis.tRNA.5.end, len))
  dis.tRNA.3.end   <- data.frame(name = "tsRNA-3'end", length.combine(dis.tRNA.3.end, len))
  dis.tRNA.CCA.end <- data.frame(name = "tsRNA-CCA end", length.combine(dis.tRNA.CCA.end, len))
  dis.tRNA.other   <- data.frame(name = "tsRNA-other", length = len, reads = dis.tRNA$reads - dis.tRNA.5.end$reads - dis.tRNA.3.end$reads - dis.tRNA.CCA.end$reads)


  dis.rRNA   <- data.frame(name = "rsRNA", length.combine(dis.rRNA, len))
  dis.piRNA  <- data.frame(name = "piRNA", length.combine(dis.piRNA, len))
  dis.unanno.match <- data.frame(name = "unanno MG", length.combine(dis.unanno.match, len))
  dis.unanno.unmatch <- data.frame(name = "unanno UMG", length.combine(dis.unanno.unmatch, len))
  dis.other  <- data.frame(name = "other", length = len, 
                           reads = dis.clean$reads - dis.miRNA$reads - dis.tRNA$reads - dis.rRNA$reads - dis.piRNA$reads - dis.unanno.match$reads - dis.unanno.unmatch$reads)

  dis.rRNA.match <- rbind(length.dis[grep("rRNAdb-rRNA_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("ensembl-Mt_rRNA_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("ensembl-rRNA_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("Rfam-rRNA_Match_Genome", length.dis$name, ignore.case = TRUE), 2:3]
                           )
  dis.rRNA.match <- data.frame(name = "rsRNA.match", length.combine(dis.rRNA.match, len))
 
  dis.rRNA.unmatch <- rbind(length.dis[grep("rRNAdb-rRNA_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("ensembl-Mt_rRNA_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("ensembl-rRNA_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3],
                            length.dis[grep("Rfam-rRNA_Unmatch_Genome", length.dis$name, ignore.case = TRUE), 2:3]
                           )
  dis.rRNA.unmatch <- data.frame(name = "rsRNA.unmatch", length.combine(dis.rRNA.unmatch, len))
  
  sum.rRNA.match   <- data.frame(name = "rRNA MG", RPM = sum(dis.rRNA.match[ ,3]))
  sum.rRNA.unmatch <- data.frame(name = "rRNA UMG", RPM = sum(dis.rRNA.unmatch[ ,3]))
  sum.rRNA <- rbind(sum.rRNA.match, sum.rRNA.unmatch)


  myLabel.rRNA     <- paste(round(sum.rRNA$RPM / sum(sum.rRNA$RPM) * 100, 1), "%", sep = "")
  sum.rRNA.match   <- data.frame(name = paste("rRNA MG (", myLabel.rRNA[1], ")" , sep = ""), RPM = sum(dis.rRNA.match[ ,3]))
  sum.rRNA.unmatch <- data.frame(name = paste("rRNA UMG (", myLabel.rRNA[2], ")" , sep = ""), RPM = sum(dis.rRNA.unmatch[ ,3]))
  sum.rRNA         <- rbind(sum.rRNA.match, sum.rRNA.unmatch)



  sum.tRNA.5.end   <- data.frame(name = "tsRNA 5' end", RPM = sum(dis.tRNA.5.end[ ,3]))
  sum.tRNA.3.end   <- data.frame(name = "tsRNA 3' end", RPM = sum(dis.tRNA.3.end[ ,3]))
  sum.tRNA.CCA.end <- data.frame(name = "tsRNA CCA end", RPM = sum(dis.tRNA.CCA.end[ ,3]))
  sum.tRNA.other   <- data.frame(name = "tsRNA other", RPM = sum(dis.tRNA.other[ ,3]))
  sum.tRNA         <- rbind(sum.tRNA.5.end, sum.tRNA.3.end, sum.tRNA.CCA.end, sum.tRNA.other)


  myLabel.tRNA     <- paste(round(sum.tRNA$RPM / sum(sum.tRNA$RPM) * 100, 1), "%", sep = "")
  sum.tRNA.5.end   <- data.frame(name = paste("tRNA-5' (", myLabel.tRNA[1], ")" , sep = ""), RPM = sum(dis.tRNA.5.end[ ,3]))
  sum.tRNA.3.end   <- data.frame(name = paste("tRNA-3' (", myLabel.tRNA[2], ")" , sep = ""), RPM = sum(dis.tRNA.3.end[ ,3]))
  sum.tRNA.CCA.end <- data.frame(name = paste("tRNA-CCA (", myLabel.tRNA[3], ")" , sep = ""), RPM = sum(dis.tRNA.CCA.end[ ,3]))
  sum.tRNA.other   <- data.frame(name = paste("tRNA-other (", myLabel.tRNA[4], ")" , sep = ""), RPM = sum(dis.tRNA.other[ ,3]))
  sum.tRNA         <- rbind(sum.tRNA.5.end, sum.tRNA.3.end, sum.tRNA.CCA.end, sum.tRNA.other)


  sum.rsRNA <- sum(dis.rRNA$reads)
  sum.piRNA <- sum(dis.piRNA$reads)
  sum.miRNA <- sum(dis.miRNA$reads)
  sum.tsRNA <- sum(dis.tRNA$reads)
  sum.other <- sum(dis.other$reads)
  sum.unanno.match <- sum(dis.unanno.match$reads)
  sum.unanno.unmatch <- sum(dis.unanno.unmatch$reads)

  stack.all  <- c()
  if(sum.miRNA > 0){
	stack.all <- rbind(stack.all, dis.miRNA)
  }
  if(sum.piRNA > 0){
	stack.all <- rbind(stack.all, dis.piRNA) 
  }
  if(sum.tsRNA > 0){
	stack.all <- rbind(stack.all, dis.tRNA)
  }  
  if(sum.rsRNA > 0){
	stack.all <- rbind(stack.all, dis.rRNA)
  }
  if(sum.other > 0){
	stack.all <- rbind(stack.all, dis.other)
  }
  if(sum.unanno.match > 0){
	stack.all <- rbind(stack.all, dis.unanno.match)
  }
  if(sum.unanno.unmatch > 0){
	stack.all <- rbind(stack.all, dis.unanno.unmatch)
  }



  stack.sep <- ggplot(stack.all, aes(x = length, y = reads)) + 
                    geom_bar(stat = "identity") + 
                    facet_grid(. ~ name) +
		    labs(x = "length", y = "RPM", title = "") + 
		    theme(axis.line = element_line()
			, panel.grid.minor = element_blank() 
			, panel.grid.major = element_blank()
			, panel.background = element_blank()
			, panel.border = element_rect(fill = "transparent")		
			, text = element_text(size = 15, color = "black")
			, axis.text.x = element_text(size = 11, color = "black")
			, axis.text.y = element_text(size = 15, color = "black")					
			)

				
  
  stack.pic <- ggplot(stack.all, aes(x = length, y = reads, fill = name)) + 
			geom_bar(stat = "identity") +
			theme(axis.line = element_line()
			, panel.grid.minor = element_blank() 
			, panel.grid.major = element_blank()
			, panel.background = element_rect(fill = "transparent")
			, plot.background = element_rect(fill = "transparent")
			, legend.title = element_blank(), legend.position = "top"
			, legend.direction =  "horizontal"
			, text = element_text(size = 15, color = "black")
			, axis.text = element_text(size = 15, color = "black")
			) +
			labs(x = "length", y = "RPM", title = "") +
			guides(fill=guide_legend(ncol=3, byrow=TRUE)) + 
			scale_fill_brewer(palette = "Paired")



						 
  myLabel.rRNA  <- paste(round(sum.rRNA$RPM / sum(sum.rRNA$RPM) * 100, 1), "%", sep = "")
  rRNA.pie <- ggplot(sum.rRNA, aes(x = "", y = RPM, fill = name)) + 
                     geom_bar(stat = "identity", width = 1) + 
                     coord_polar(theta = "y") + 
			        labs(x = "", y = "", title = "") + 
                                theme(panel.background = element_rect(fill = "transparent")
				, plot.background = element_rect(fill = "transparent")
			    	, axis.ticks = element_blank() 
				, panel.grid = element_blank() 
				, panel.border = element_blank()
                                , axis.text.x = element_blank()
				, text = element_text(size = 15, color = "black")
                    		, legend.title = element_blank()
				, legend.position = "bottom"
				, legend.direction =  "vertical"
				, plot.margin = margin(-80, 10, -30, 10)
				, legend.margin = margin(-50, 0, -30, 0)
				) 

					
  myLabel.tRNA  <- paste(round(sum.tRNA$RPM / sum(sum.tRNA$RPM) * 100, 1), "%", sep = "")
  tRNA.pie <- ggplot(sum.tRNA, aes(x = "", y = RPM, fill = name)) + 
                     geom_bar(stat = "identity", width = 1) + 
                     coord_polar(theta = "y") + 
			        labs(x = "", y = "", title = "") + 
                                theme(panel.background = element_rect(fill = "transparent")
				, plot.background = element_rect(fill = "transparent")
			    	, axis.ticks = element_blank() 
				, panel.grid = element_blank() 
				, panel.border = element_blank()
                                , axis.text.x = element_blank()
				, text = element_text(size = 15, color = "black")
                    		, legend.title = element_blank()
				, legend.position = "bottom"
				, legend.direction =  "vertical"
				, plot.margin = margin(-80, 10, 0, 10)
				, legend.margin = margin(-50, 0, 0, 0)
				) 
  
  
  pdf(paste(file.address, file.name, "_result/", file.name, "_sncRNA_distribution.pdf", sep=""), width = 9, height = 12)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 3, heights = unit(c(1, 4, 4, 4), "null"))))
  pdf.title <- paste("Small RNAs Length Distribution of ", file.name, sep = "")
  vplayout <- function(x,y)
	viewport(layout.pos.row = x, layout.pos.col = y)
  grid.text(pdf.title, vp = vplayout(1, 1:3), gp = gpar(fontsize = 20))
  if((sum.rsRNA) > 0 && (sum.rsRNA)){
	print(stack.pic, vp = vplayout(2:3, 1:2))	
	print(rRNA.pie, vp = vplayout(2, 3))
	print(tRNA.pie, vp = vplayout(3, 3))
	print(stack.sep, vp = vplayout(4, 1:3))
  }else if ((sum.rsRNA) > 0){
	print(stack.pic, vp = vplayout(2:3, 1:3))
	subvp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.8)
	print(rRNA.pie, vp = subvp)
	print(stack.sep, vp = vplayout(4, 1:3))
  }else if ((sum.tsRNA) > 0){
	subvp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.8)
	print(tRNA.pie, vp = subvp)
	print(stack.sep, vp = vplayout(4, 1:3))
  }else{
	print(stack.pic, vp = vplayout(2:3, 1:3))
	print(stack.sep, vp = vplayout(4, 1:3))
  }
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

length.stack(file.address = args[1], file.name = args[2])
