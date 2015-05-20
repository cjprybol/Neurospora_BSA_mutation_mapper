library("ggplot2")
library("gridExtra")

args <- commandArgs(trailingOnly=TRUE)

data = read.table(args[1],sep='\t',colClasses=c(rep("numeric",5)),header=TRUE, stringsAsFactors=FALSE)
save_path = args[2]
file_base = args[3]

# the supercontig/chromsome #, position on that contig, # of reads matching reference snp, # matching alternate snp, # matching neither (mismatch)
colnames(data) <- c("CONTIG", "POS", "REF", "ALT", "MIS")

data$TOTAL = data$REF+data$ALT
data$RATIO = ( data$REF / data$TOTAL ) * 100
data$KB = data$POS/1000
data = na.omit(data)

chrom1 <- data[ which(data$CONTIG==1),]
chrom2 <- data[ which(data$CONTIG==2),]
chrom3 <- data[ which(data$CONTIG==3),]
chrom4 <- data[ which(data$CONTIG==4),]
chrom5 <- data[ which(data$CONTIG==5),]
chrom6 <- data[ which(data$CONTIG==6),]
chrom7 <- data[ which(data$CONTIG==7),]

alpha_max = max(data$REF + data$ALT)

p1 <- ggplot(chrom1, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 1") +
		stat_smooth()
		
		
p2 <- ggplot(chrom2, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 2") +
		stat_smooth()


p3 <- ggplot(chrom3, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 3") +
		stat_smooth()

		
p4 <- ggplot(chrom4, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 4") +
		geom_vline(xintercept = 3481, linetype = 2, colour="blue") +
		stat_smooth()
		
p5 <- ggplot(chrom5, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 5") +
		stat_smooth()
		
p6 <- ggplot(chrom6, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 6") +
		stat_smooth()
		
p7 <- ggplot(chrom7, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
		labs(x = "Position (kb)", title="Chromosome 7") +
		stat_smooth()

ggsave(filename=paste(save_path, file_base, "_contig_1.png", sep=""), plot=p1, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_2.png", sep=""), plot=p2, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_3.png", sep=""), plot=p3, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_4.png", sep=""), plot=p4, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_5.png", sep=""), plot=p5, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_6.png", sep=""), plot=p6, width=16, height=9, units="in")
ggsave(filename=paste(save_path, file_base, "_contig_7.png", sep=""), plot=p7, width=16, height=9, units="in")


point_size <- 0.7

p1 <- ggplot(chrom1, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 1") +
		stat_smooth()
		
p1_legend <- ggplot(chrom1, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point() + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.box = "horizontal") +
		labs(x = "Position (kb)", title="Chromosome 1") +
		stat_smooth()
		
p2 <- ggplot(chrom2, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 2") +
		stat_smooth()


p3 <- ggplot(chrom3, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 3") +
		stat_smooth()

		
p4 <- ggplot(chrom4, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 4") +
		geom_vline(xintercept = 3481, linetype = 2, colour="blue") +
		stat_smooth()
		
p5 <- ggplot(chrom5, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 5") +
		stat_smooth()
		
p6 <- ggplot(chrom6, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 6") +
		stat_smooth()
		
p7 <- ggplot(chrom7, aes(x = KB, y = RATIO, colour = RATIO, alpha = TOTAL)) + 
		geom_point(size = point_size) + 
		ylim(0, 100) +
		scale_colour_gradient2(name = "% Oak Ridge", low = "#0000FF", mid = "#404040", high ="#FF0000", midpoint = median(data$RATIO), space = "rgb", limits = c(0,100)) +
		scale_alpha_continuous(name = "Read Count ~\nTransparency", limits = c(0,alpha_max), guide = guide_legend(reverse=TRUE)) +
		theme(axis.ticks = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
		labs(x = "Position (kb)", title="Chromosome 7") +
		stat_smooth()


g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

legend <- g_legend(p1_legend)
multiplot <- arrangeGrob(p1,p2,p3,p4,p5,p6,p7,legend, ncol = 2)
ggsave(filename=paste(save_path, file_base, "_multiplot.png", sep=""), plot=multiplot, width=8.5, height=11, units="in")
