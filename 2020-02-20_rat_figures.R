# this is a file to plot all the figures for the paper
# it should be zipped with all the data necessary to generate the the figures
# and any additional scripts required to process the data
# make sure we aren't relying on data or variables stored in memory
#rm(list=ls())

### this is your working directory, which should contain a "fig_data" folder
setwd("~/Documents/Manuscripts/Current/rat_diets/figures")

### change these to booleans to suppress data loading or plotting for specific figures
preprocess <- 0 # preprocess / load data for plotting or don't

### if you want to use raw data (unfiltered) set this to zero
filter.fp <- 0

fig1 <- 0 # maps of locations
fig2 <- 0 # hist of read lengths amd barcode splits
fig3 <- 0 # biplot of alignment length and e-value and blast hits as qual changes
fig4 <- 0 # hist of readlength:alignlength hist for rats, bacteria, diet items
fig5 <- 0 # boxplot of diet item diversity, must have tables produced for this figure to be produced.
fig6 <- 0 # heatmap of diet item frequency, must have tables produced for this figure to be produced.
fig7 <- 0 # scatter plot of CAP analysis for filtered families
fig8 <- 0 # plot of genera with fraction of sequences that are genomic vs. rrna/mtdna/chloro/microsatellite
figS1 <- 0 # biplots of read quals / lengths per barcode
figS2 <- 1 # read quality vs identity and alignment length
print.tables <- 0 # tables of taxa numbers (order, family, and genus)
print.data <- 0 # table of all blast hits (datafile)
new.figure <- 1

### if you want to change the parameters, e.g. cutoffs
### to say this is a genus-level classification
### these are here:
genus.pid <- 82.5
genus.ratio <- 0.55
family.pid <- 77.5
family.ratio <- 0.1
family.length <- 150

###################################
### should be able to ignore below here
###################################
library(plyr)
library(RColorBrewer)

### preprocessing of data
if(preprocess) {
	
	data <- read.csv(file="./fig_data/MegaBestHit_17_6_18.csv", header=T)
	jandata <- subset(data, SequenceRun=="January")
	janbcdata <- subset(jandata, Barcode!="J Unbarc")
	mardata <- subset(data, SequenceRun=="March")
	marbcdata <- subset(mardata, Barcode!="M Unbarc")
	cat("jan hits: ", dim(janbcdata)[1])
	cat("\nmar hits: ", dim(marbcdata)[1])
	allbcdata <- rbind(janbcdata, marbcdata)
	hq <- subset(allbcdata, (evalue<1e-20 & alignment.length>=100))
	
	jan.q <- read.table(file="./fig_data/Jan_qual_scores.txt", sep="\t")
	mar.q <- read.table(file="./fig_data/Mar_qual_scores.txt", sep="\t")
	all.q <- read.table(file="./fig_data/all_reads_lengthqual.txt", sep="\t")
		
	hq$Ratio <- hq$alignment.length/hq$ReadLength
	hq$LogReadLength <- log10(hq$ReadLength)
	hq$Logalignment.length <- log10(hq$alignment.length)
	
	### need to parse kingdom and assign taxonid only, not names
	k <- read.table(file="./fig_data/MEGAN_kingdom.txt")
	colnames(k) <- c("query.def2", "Megan.Kingdom")
	hq <- join(hq, k, by="query.def2", match="first")
	# kingdom
	p <- read.table(file="./fig_data/MEGAN_phylum.txt")
	colnames(p) <- c("query.def2", "Megan.Phylum")
	hq <- join(hq, p, by="query.def2", match="first")
	# class
	c <- read.table(file="./fig_data/MEGAN_class.txt")
	colnames(c) <- c("query.def2", "Megan.Class")
	hq <- join(hq, c, by="query.def2", match="first")
	# order
	o <- read.table(file="./fig_data/MEGAN_order.txt")
	colnames(o) <- c("query.def2", "Megan.Order")
	hq <- join(hq, o, by="query.def2", match="first")
	#family
	f <- read.table(file="./fig_data/MEGAN_family.txt")
	colnames(f) <- c("query.def2", "Megan.Family")
	hq <- join(hq, f, by="query.def2", match="first")
	#genus
	g <- read.table(file="./fig_data/MEGAN_genus.txt")
	colnames(g) <- c("query.def2", "Megan.Genus")
	hq <- join(hq, g, by="query.def2", match="first")
	
	# exclude rodents, humans, and bacteria
	# Euks are 2759
	hq.diet <- subset(hq, Megan.Kingdom==2759 )
	# Bacterial kingdom is 2
	hq.micro <- subset(hq, Megan.Kingdom==2)
	cat("\nfract. bacterial", round(dim(hq.micro)[1]/(dim(hq.micro)[1]+dim(hq.diet)[1] ),3), "\n" )
	hq.lactob <- subset(hq, Megan.Genus==1578)
	cat("fract. Lactobacillus", round(dim(hq.lactob)[1]/dim(hq.micro)[1],3), "\n" )
	hq.lactoc <- subset(hq, Megan.Genus==1357)
	cat("fract. Lactococcus", round(dim(hq.lactoc)[1]/dim(hq.micro)[1],3), "\n" )
	hq.pseudom <- subset(hq, Megan.Genus==286)
	cat("fract. Pseudomonas", round(dim(hq.pseudom)[1]/dim(hq.micro)[1],3), "\n" )

	# Rodentia Order is 9989, Mammalia is 40674
	# orders and classes that are NA are not excluded
	hq.diet <- subset(hq.diet, !(Megan.Order==9989) | is.na(Megan.Order))
	hq.diet <- subset(hq.diet, !(Megan.Class==40674) | is.na(Megan.Class))
	# for cases in which the Classes or orders are not defined
	hq.diet <- subset(hq.diet, !(Genus %in% c("Rattus", "rat", "Mus", "Mouse", "Peromyscus", "Myotis", "Castor")) )
	# fix some cases in which Families or Orders are not identified
	# do we need the equals sign here?
	if(1){
		r <- which(hq.diet$Genus=="Syphacia")
		hq.diet$Megan.Family[r] <- 51026
		hq.diet$Megan.Order[r] <- 70426
		hq.diet$Megan.Class[r] <- 119089
		r <- which(hq.diet$Genus=="Hymenolepis")
		hq.diet$Megan.Family[r] <- 6214
		hq.diet$Megan.Order[r] <- 6201
		hq.diet$Megan.Class[r] <- 6199
		hq.diet$Megan.Phylum[r] <- 6157
		r <- which(hq.diet$Megan.Order==70426)
		hq.diet$Megan.Phylum[r] <- 6231
	}

	# Rattus genus is 10114
	hq.rat <- subset(hq, Megan.Genus==10114)
	cat("Rattus reads:", dim(hq.rat)[1], " Median pid:", median(hq.rat$X..identity), " Median ratio:", round(median(hq.rat$Ratio),3))
	# Mus genus is 10088
	hq.mouse <- subset(hq, Megan.Genus==10088)
	cat("\nMus reads:", dim(hq.mouse)[1], " Median pid:", median(hq.mouse$X..identity), " Median ratio:", round(median(hq.mouse$Ratio),3))
	cat("\nincluded mouse: ", length(which(hq.mouse$Ratio>genus.ratio & hq.mouse$X..identity>genus.pid)))
	cat("\nincluded rat: ", length(which(hq.rat$Ratio>genus.ratio & hq.rat$X..identity>genus.pid)), "\n")
	hq.diet$Megan.Genus <- as.factor(hq.diet$Megan.Genus)
	hq.diet$Megan.Family <- as.factor(hq.diet$Megan.Family)
	hq.diet$Megan.Order <- as.factor(hq.diet$Megan.Order)
	hq.diet$Megan.Class <- as.factor(hq.diet$Megan.Class)
	hq.diet$Megan.Phylum <- as.factor(hq.diet$Megan.Phylum)
	hq.diet$Megan.Kingdom <- as.factor(hq.diet$Megan.Kingdom)
	
	if(filter.fp) {
		# get rid of FP genera and families
		incorrect.genus <- which(hq.diet$X..identity<genus.pid | hq.diet$Ratio<genus.ratio)
		hq.diet$Megan.Genus[incorrect.genus] <- NA
		
		incorrect.family <- which(hq.diet$X..identity<family.pid | hq.diet$Ratio<family.ratio | hq.diet$alignment.length<family.length)
		hq.diet$Megan.Family[incorrect.family] <- NA
	}
	tg <- round(1 - sum(is.na(hq.diet$Megan.Genus))/dim(hq.diet)[1],3)
	tf <- round(1 - sum(is.na(hq.diet$Megan.Family))/dim(hq.diet)[1],3)
	to <- round(1 - sum(is.na(hq.diet$Megan.Order))/dim(hq.diet)[1],3)
	tc <- round(1 - sum(is.na(hq.diet$Megan.Class))/dim(hq.diet)[1],3)
	tp <- round(1 - sum(is.na(hq.diet$Megan.Phylum))/dim(hq.diet)[1],3)
	cat("Genus, Family, Order, Class, Phylum, Kingdom (in order):", tg, tf, to, tc, tp)
	# fastq files are here, load full paths and short names
	jan.qual <- dir("./fig_data/2018_06_14_albacore/2017_01_24_albacore_fastq", full.names=T, pattern="qual.txt$")
	jan.qual.names <- dir("./fig_data/2018_06_14_albacore/2017_01_24_albacore_fastq", full.names=F, pattern="qual.txt$")
	mar.qual <- dir("./fig_data/2018_06_14_albacore/2017_03_17_albacore_fastq", full.names=T, pattern="qual.txt$")
	mar.qual.names <- dir("./fig_data/2018_06_14_albacore/2017_03_17_albacore_fastq", full.names=F, pattern="qual.txt$")
	qual.files <- c(jan.qual, mar.qual)
	jan.qual.names <- gsub("_qual.txt","",jan.qual.names)
	mar.qual.names <- gsub("_qual.txt","",mar.qual.names)
	qual.names <- c(jan.qual.names, mar.qual.names)
	
	taxa.table <- read.table(file="./fig_data/abridged_taxonomy_norank.tab", sep="\t")
	colnames(taxa.table) <- c("taxid", "taxlevel", "name")
	all.rat.bc <- c("JNB01","JNB02","JNB03","JNB04","JNB05","JNB06","JNB07","JNB08","JNB09","JNB10","JNB11","JNB12", "MNB01","MNB02","MNB03","MNB04","MNB05","MNB06","MNB07","MNB08","MNB09","MNB10","MNB11","MNB12")
		
	all.rat.locs <- c("OS2","OS13","OS7","OS5","PW7","PW6","PW3","PW5","AP11","AP9","AP5","AP8","AP10","OS10","PW1","AP12","OS11","PW2","AP2","OS12","PW4","AP3","OS14","PW8")
	all.rat.locs.names <- c("OB2","OB13","OB7","OB5","LB7","LB6","LB3","LB5","WP11","WP9","WP5","WP8","WP10","OB10","LB1","WP12","OB11","LB2","WP2","OB12","LB4","WP3","OB14","LB8")
	all.rat.loc.lookup <- matrix(c(all.rat.locs,all.rat.locs.names),ncol=2,byrow=F)
	colnames(all.rat.loc.lookup) <- c("Sample","paper.names")

}

###############################################
###############################################
### Figure 1 is a map of the locations
###############################################
###############################################
if(fig1) {

library(maps)
library(mapdata)
library(maptools)
library(ggplot2)
library(ggmap)

### check here for styling adjustments
### https://developers.google.com/maps/documentation/static-maps/styling#features

### location of samples
rat.locs <- read.table(file="./fig_data/rat_locations.txt", header=T)

### Get map of whole area
aklmap <- get_googlemap(center = 'riverhead, new zealand', maptype="roadmap", zoom=11, style="feature:administrative|element:labels|visibility:off&style=feature:road|element:labels|visibility:off&style=feature:road|element:geometry|visibility:off")

### set zoom levels for maps and rectangle drawing
LB_zoom <- c(174.702, 174.749, -36.685, -36.663)
WP_zoom <- c(174.496, 174.520, -36.895, -36.87)

### render and save overview map
ggmap(aklmap) + geom_point(data = rat.locs, mapping = aes(x = lon, y = lat, colour = factor(place)), size=1.5, show.legend=F) + labs(x="Longitude", y="Latitude") + xlim(174.4, 174.755) + ylim(-36.905,-36.62) + scale_color_manual(values=c("red3","blue3","green4")) + annotate("text", x=c(174.475, 174.704, 174.73), y=c(-36.862, -36.649, -36.71), label=c("Waitakere Ranges\nRegional Park", "Okura Bush\nScenic Reserve", "Long Bay\nRegional Park\nWetlands"), colour="black", size=4.5) + annotate("text", x=174.74, y=-36.843, label="Auckland", colour="green4", size=5) + annotate("rect", xmin=LB_zoom[1], xmax=LB_zoom[2], ymin=LB_zoom[3], ymax=LB_zoom[4], colour="black", fill=NA, linetype=3) + annotate("rect", xmin=WP_zoom[1], xmax=WP_zoom[2], ymin=WP_zoom[3], ymax=WP_zoom[4], colour="black", fill=NA, linetype=3)
ggsave("./Fig1.pdf", useDingbats=FALSE)

### Long Bay area
### get data
zoom1map <- get_googlemap(center = 'long bay, new zealand', maptype="terrain", zoom=13, style="feature:administrative|element:labels|visibility:off&style=feature:road|element:labels|visibility:off&style=feature:poi|element:labels|visibility:off")

### render and save LB zoomed map
ggmap(zoom1map) + geom_point(data = rat.locs, mapping = aes(x = lon, y = lat, colour = factor(place)), size=4.5, show.legend=F) + annotate("text", x=rat.locs[6:21,4]-0.0021, y=rat.locs[6:21,3], label=rat.locs[6:21,1], colour="black", size=6.5) + xlim(LB_zoom[1], LB_zoom[2]) + ylim(LB_zoom[3], LB_zoom[4]) + scale_color_manual(values=c("red3","blue3","green4")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, -1, -1), 'lines')) + xlab('') + ylab('')
ggsave("./Fig1_inset1_small.pdf", useDingbats=FALSE)

### Waitakere Area
### get data
zoom2map <- get_googlemap(center = 'Waitakere Reservoir, new zealand', maptype="terrain", zoom=13, style="feature:administrative|element:labels|visibility:off&style=feature:road|element:labels|visibility:off")

### render and save WP zoomed map
ggmap(zoom2map) + 
geom_point(data = rat.locs, mapping = aes(x = lon, y = lat, colour = factor(place)), size=1, show.legend=F) + 
annotate("text", x=rat.locs[1:5,4]-0.0021, y=rat.locs[1:5,3], label=rat.locs[1:5,1], colour="black", size=6.5) + 
xlim(WP_zoom[1], WP_zoom[2]) + 
ylim(WP_zoom[3], WP_zoom[4]) + 
scale_color_manual(values=c("red3","blue3","green4")) + 
theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, -1, -1), 'lines')) + 
xlab('') + 
ylab('')

ggsave("./Fig1_inset2_small.pdf", useDingbats=FALSE)
cat("Figure 1\n")
}
###############################################
###############################################
### End Fig1 plotting
###############################################
###############################################


###############################################
###############################################
### Figure 2 is a histogram of the read lengths
### and BC numbers for ALL reads from Jan and March
### in constrast to S1, which is lengths and quals 
### for each rat
###############################################
###############################################
if(fig2) {
	legend.cex <- 1.6
	pdf(file="./Fig2_reverse.pdf", height=3.5, width=10.5)
	par(mfrow=c(1,3))
	par(las=1)
	# subset and hist of Jan lengths
	janlen <- read.table(file="./fig_data/2018_06_14_albacore/Jan_qual_scores.txt", sep="\t")
	
	## new
	par(mar=c(5,5,4,4))
	bcdata <- read.table(file="./fig_data/jan_read_numbers.txt", header=T)
	ylim <- c(0.5,16)
	
	barplot(sort(bcdata$total_reads), names.arg=bcdata$rat[order(bcdata$total_reads)], ylim=ylim, xlim=c(0,2e4),xaxt="n", xlab="Number of reads (thousands)", horiz=T, cex.lab=1.2, cex.axis=1.2)
	axis(1, at=c(0,5e3,10e3,15e3,20e3), labels=c("0","5", "10", "15", "20"), cex.axis=1.2)
	mtext("a", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)
	cat("\njan reads:", sum(bcdata$total_reads), "\njan no bc:", bcdata[13,3]/sum(bcdata$total_reads))
	
	bcdata <- read.table(file="./fig_data/mar_read_numbers.txt", header=T)
	ylim <- c(0.5,16)
	barplot(sort(bcdata$total_reads), names.arg=bcdata$rat[order(bcdata$total_reads)], ylim=ylim, xlim=c(0,3e4), xaxt="n", xlab="Number of reads (thousands)", horiz=T, cex.lab=1.2, cex.axis=1.2)
	axis(1, at=c(0,10e3,20e3,30e3), labels=c("0", "10", "20", "30"), cex.axis=1.2)
	mtext("b", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)
	
	 cat("\nmar reads:", sum(bcdata$total_reads), "\nmar no bc:", bcdata[13,3]/sum(bcdata$total_reads))
	 
	 ## new
	 par(mar=c(5,6,4,3))
	hl <- hist(log10(janlen[,2]), breaks=seq(0,6,by=0.05), plot=F)
	ylim <- c(0,0.14)
	plot(-1,-1, xlim=c(2.1,4.3), ylim=ylim, xaxt="n", bty="n", xlab="Read length (base pairs)", ylab="", cex.lab=1.2, cex.axis=1.2)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(0,0.5,1,0.2), border="dark blue")
	axis(1,at=log10(c(200,2000,20000)),labels=c("200", "2K", "20K"), cex=1.2)
	mtext("c", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	mtext("Fraction of reads", side=2, cex=0.9, adj=0.6, padj=-4.5, las=0)
	text(3.3, 0.035, labels="January", col="dark blue", cex=1.1)
	cat("\njan median length: ", median(janlen[,2]))
	# subset and hist of Mar lengths
	marlen <- read.table(file="./fig_data/2018_06_14_albacore/Mar_qual_scores.txt", sep="\t")
	
	hl <- hist(log10(marlen[,2]), breaks=seq(0,7,by=0.05), plot=F)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(1,0.5,0,0.2), border="dark orange")
	text(2.3, 0.055, labels="March", col="dark orange", cex=1.1)
	cat("\nmar median length: ", median(marlen[,2]))
	
	dev.off()
	cat("\nFigure 2\n")
}
###############################################
###############################################
### End Fig2 plotting
###############################################
###############################################

###############################################
###############################################
### Figure 3 is a biplot of alignment length 
### and e-value
###############################################
###############################################
if(fig3) {
	pdf(file="./Fig3_upper.pdf", height=4, width=8)
	par(mfrow=c(1,2))
	layout(matrix(c(1,0,4,4,2,3,4,4),nrow=2, byrow=T), rep(c(1,0.3),4), rep(c(0.3,1),4) )
	legend.cex <- 1.5
	par(las=1)
	
	par(mar=c(0,7,2,0))
	hl <- hist(log10(allbcdata$alignment.length), breaks=seq(0,5,by=0.025), plot=F)
	plot(-1,-1, xlim=c(1.5,4), ylim=c(0,0.12), xaxt="n", bty="n", yaxt="n", xlab="", ylab="")
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col="light blue")
	
	par(mar=c(5,7,0,0))
	
	plot(allbcdata$alignment.length, (allbcdata$evalue+1e-190), log="xy", cex=0.3, col=rgb(0,0,0,0.05), pch=19, xlim=c(10^1.5,1e4), ylim=c(1e-190,1), xaxt="n", yaxt="n", xlab="Alignment length (base pairs)", ylab="")
	segments(100,1e-200,100,1e-20, lty=2, col="red", lwd=0.8)
	segments(100,1e-20,1e5,1e-20, lty=2, col="red", lwd=0.8)
	
	axis(1,at=c(50,500,5000))
	axis(2,at=c(1e-0, 1e-40, 1e-80, 1e-120, 1e-160, 1e-200))
	par(las=3)
	mtext("Blast E-value", side=2, line=4.5, cex=0.7)
	
	par(las=1)
	par(mar=c(5,0,0,2))
	hl <- hist(log10(allbcdata$evalue+1e-190), breaks=seq(-191,1,by=2), plot=F)
	plot(-1,-1, xlim=c(0,0.2), ylim=c(-190,0), xaxt="n", bty="n", yaxt="n", xlab="", ylab="")
	polygon( c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), c(sort(rep(hl$breaks,2)), hl$breaks[1]), col="light blue")
	
	mtext("A", side=2, cex=legend.cex, at=50, adj=18)
### End Fig3A plotting

### Figure 3B is plot of blast hits as qual changes

	par(mar=c(5,7,7,3))
	bin.size <- 0.02
	jan.bin.hits <- vector()
	jan.bin.total <- vector()
	mar.bin.hits <- vector()
	mar.bin.total <- vector()

	# get long reads
	jan.blast.hits <- janbcdata[c('evalue','quals_mean')]
	mar.blast.hits <- marbcdata[c('evalue','quals_mean')]
	jan.blast.hits[,1] <- jan.blast.hits[,1] + 1e-190
	mar.blast.hits[,1] <- mar.blast.hits[,1] + 1e-190
	jan.blast.hits[,3] <- jan.blast.hits$evalue<1e-20
	mar.blast.hits[,3] <- mar.blast.hits$evalue<1e-20
	qual.bins <- seq(0.75,0.95,by=bin.size)
	#i <- 1
	for (i in 1:(length(qual.bins)-1)) {
		jan.blast.sub <- subset(jan.blast.hits, quals_mean >= qual.bins[i] & quals_mean < qual.bins[i+1])
		jan.bin.total[i] <- length(which(jan.q[,3] >= qual.bins[i] & jan.q[,3] < qual.bins[i+1]))
		jan.bin.hits[i] <- sum(jan.blast.sub[,3])
		
		mar.blast.sub <- subset(mar.blast.hits, quals_mean >= qual.bins[i] & quals_mean < qual.bins[i+1])
		mar.bin.total[i] <- length(which(mar.q[,3] >= qual.bins[i] & mar.q[,3] < qual.bins[i+1]))
		mar.bin.hits[i] <- sum(mar.blast.sub[,3])

	}
	plot(qual.bins[1:(length(qual.bins)-1)]-bin.size/2, jan.bin.hits/jan.bin.total, pch=21, cex=1.2, bg="light blue",xlab="Mean accuracy of read (binned)", ylab="Fraction of reads in bin with\na high quality Blast hit", ylim=c(-0.02,0.4))
	points(qual.bins[1:(length(qual.bins)-1)]-bin.size/2, mar.bin.hits/mar.bin.total, pch=21, cex=1.2, bg="orange")
	text(qual.bins[1:(length(qual.bins)-1)]-bin.size/2, jan.bin.hits/jan.bin.total-0.02,labels=round(jan.bin.total/1000,1), cex=0.9)
	text(qual.bins[1:(length(qual.bins)-1)]-bin.size/2, mar.bin.hits/mar.bin.total+0.02,labels=round(mar.bin.total/1000,1), cex=0.9)

	mtext("B", side=2, cex=legend.cex, at=0.5, adj=2)
	dev.off()
	
	if(0) {
	### to get numbers
	lowqual <- 0.75
	total <- all.q[which(all.q[,3]>0),3]
	blast.hits <- data[which(data$quals_mean>0),]
	blast.hits <- blast.hits[which(blast.hits$quals_mean<lowqual),]
	fractbad <- sum(blast.hits$evalue<1e-20)/sum(total<=lowqual)
	
	hiqual <- 0.92
	total <- all.q[which(all.q[,3]>hiqual),3]
	blast.hits <- data[which(data$quals_mean>hiqual),]
	blast.hits <- blast.hits[which(blast.hits $quals_mean<1),]
	fractgood <- sum(blast.hits$evalue<1e-20)/sum(total<=1)	
	cat("fract low qual seqs (", lowqual, ") with hits: ", fractbad, "\n", "fract high qual (", hiqual, ") seqs with hits:", fractgood, "\n")
	}
	cat("Figure 3\n")
	
}
###############################################
###############################################
### End Fig3 plotting
###############################################
###############################################

###############################################
###############################################
### Figure 4 is a histogram of the read 
### lengths/alignment ratio and percent 
### identity for bacteria, rats, and diet
###############################################
###############################################
if(fig4) {
	legend.cex <- 2
	pdf(file="./Fig4.pdf", height=5.5, width=11.5)
	par(mfrow=c(1,2))
	par(mar=c(5,7,4,1.5))
	par(las=1)
		
	### identities
	hl <- hist(hq.diet$X..identity, breaks=seq(60,100,by=0.5), plot=F)
	ylim <- c(0,0.065)
	plot(-1,-1, xlim=c(70,100), ylim=c(0,0.065), bty="n", xlab="Percent identity", ylab="", cex.lab=1.2, cex.axis=1.2)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(0,0.5,1,0.2), border="dark blue")
	text(96.5, 0.015, labels="Diet items", col="dark blue", cex=1.1)
	mtext("Fraction of reads", side=2, cex=1.2, adj=0.5, padj=-5.2, las=0)

	hl <- hist(hq.rat$X..identity, breaks=seq(60,100,by=0.5), plot=F)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(0,0,0,0.2), border="black")
	text(92, 0.044, labels="Rattus", col="black", cex=1.1, font=3)
	
	hl <- hist(hq.mouse$X..identity, breaks=seq(60,100,by=0.5), plot=F)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(1,0,0,0.2), border="dark orange")
	text(74.5, 0.03, labels="Mus", col="dark orange", cex=1.1, font=3)
	abline(v=c(genus.pid, family.pid), lty=2)
	mtext("a", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)

	### ratios
	hl <- hist(hq.diet$Ratio, breaks=seq(0,1.2,by=0.02), plot=F)
	ylim <- c(0,0.055)
	plot(-1,-1, xlim=c(0,1.2), ylim=ylim, bty="n", xlab="Ratio of alignment length to read length", ylab="", cex.lab=1.2, cex.axis=1.2)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(0,0.5,1,0.2), border="dark blue")
	text(0.935, 0.039, labels="Diet items", col="dark blue", cex=1.1)
	mtext("Fraction of reads", side=2, cex=1.2, adj=0.5, padj=-5.2, las=0)

	hl <- hist(hq.rat$Ratio, breaks=seq(0,1.2,by=0.02), plot=F)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(0,0,0,0.2), border="black")
	text(0.5, 0.048, labels="Rattus", col="black", cex=1.1, font=3)
	
	hl <- hist(hq.mouse$Ratio, breaks=seq(0,1.2,by=0.02), plot=F)
	polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col=rgb(1,0,0,0.2), border="dark orange")
	text(0.19, 0.02, labels="Mus", col="dark orange", cex=1.1, font=3)
	abline(v=genus.ratio, lty=2)
	mtext("b", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)

	dev.off()
	cat("Figure 4\n")
}	
###############################################
###############################################
### End Fig4 plotting
###############################################
###############################################

###############################################
###############################################
### New figure of family-level quantifications
###############################################
###############################################

if(new.figure) {
	
	pdf(file="./Fignew.pdf", height=5, width=5.5)
	par(las=1)
	
	hq.filt <- hq[which(hq$X..identity>77.5),]
	hq.filt <- hq[which(hq.filt$Ratio>0.1),]
	hq.rodent <- subset(hq, Megan.Order==9989)
	hq.murid <- subset(hq.rodent, Megan.Family==10066)
	hq.cricet <- subset(hq.rodent, Megan.Family==337677)
	hq.spalac <- subset(hq.rodent, Megan.Family==337664)
	
	if(0) {
		hq.mus <- subset(hq.rodent, Megan.Genus==10088)
		hq.rat <- subset(hq.rodent, Megan.Genus==10114)
		plot(hq.rat$X..identity, hq.rat$Ratio, pch=19, cex=0.5, col=rgb(0,0,0,0.3), xlab="Percent identity", ylab="Ratio of read length to alignment length")
		points(hq.mus$X..identity, hq.mus$Ratio, pch=19, cex=0.5, col=rgb(1,0.4,0,0.2))

	}
	hq.filt.rodent <- subset(hq.filt, Megan.Order==9989)
	hq.filt.murid <- subset(hq.filt.rodent, Megan.Family==10066)
	hq.filt.cricet <- subset(hq.filt.rodent, Megan.Family==337677)
	hq.filt.spalac <- subset(hq.filt.rodent, Megan.Family==337664)

	
	plot(hq.murid$X..identity, hq.murid$Ratio, pch=19, cex=0.5, col=rgb(0,0,0,0.3), xlab="Percent identity", ylab="Ratio of read length to alignment length")
	points(hq.cricet$X..identity, hq.cricet$Ratio, pch=19, cex=0.8, col=rgb(1,0.4,0,0.8))
	points(hq.spalac $X..identity, hq.spalac $Ratio, pch=19, cex=0.8, col=rgb(0,0.9,0.5,0.8))
	tp <- (dim(hq.cricet)[1]+dim(hq.spalac)[1])/(dim(hq.cricet)[1]+dim(hq.spalac)[1]+dim(hq.murid)[1])
	tp.filt <- (dim(hq.filt.cricet)[1]+dim(hq.filt.spalac)[1])/(dim(hq.filt.cricet)[1]+dim(hq.filt.spalac)[1]+dim(hq.filt.murid)[1])
	dev.off()
	cat("New figure\n")
}
###############################################
###############################################
### END NEW FIG
###############################################
###############################################

###############################################
###############################################
### These tables contain info on which rats
### ate which taxa, can be used for PCA
### For ease, this section is bundled with Fig 5 and 6
###############################################
###############################################

if(print.tables) {
	all.rat.family.diets <- matrix()
	all.rat.genus.diets <- matrix()
	all.rat.family.ngen <- matrix()
	all.rat.order.nfam <- matrix()
	all.rat.order.ngen <- matrix()
	total.genus <- vector()
	total.family <- vector()
	total.order <- vector()
	total.rat <- vector()
	
	tax.hier <- unique(hq.diet[c("Megan.Family","Megan.Order","Megan.Phylum")])
	tax.hier[,4] <- taxa.table[match(tax.hier$Megan.Family, taxa.table$taxid),3]
	tax.hier[,5] <- taxa.table[match(tax.hier$Megan.Order, taxa.table$taxid),3]
	tax.hier[,6] <- taxa.table[match(tax.hier$Megan.Phylum, taxa.table$taxid),3]
	colnames(tax.hier)[4:6] <- c("Family.Name","Order.Name","Phylum.Name")
	
	for(i in 1:24) {
		hq.diet.rat <- subset(hq.diet, Rat==all.rat.locs[i])
		
		# get table of numbers in each family (or other taxon level)
		family.taxa <- as.data.frame(table(hq.diet.rat$Megan.Family))
		genus.taxa <- as.data.frame(table(hq.diet.rat$Megan.Genus))
		order.taxa <- as.data.frame(table(hq.diet.rat$Megan.Order))

		colnames(family.taxa) <- c("taxid",all.rat.locs[i])
		colnames(genus.taxa) <- c("taxid",all.rat.locs[i])
		colnames(order.taxa) <- c("taxid",all.rat.locs[i])
		
		total.rat[i] <- dim(hq.diet.rat)[1]
				
		if(i==1) {
			all.rat.family.diets <- merge(family.taxa, taxa.table, by="taxid")
			all.rat.genus.diets <- merge(genus.taxa, taxa.table, by="taxid")
			all.rat.order.diets <- merge(order.taxa, taxa.table, by="taxid")
		}
		else {
			all.rat.family.diets <- merge(family.taxa, all.rat.family.diets, by="taxid")
			all.rat.genus.diets <- merge(genus.taxa, all.rat.genus.diets, by="taxid")
			all.rat.order.diets <- merge(order.taxa, all.rat.order.diets, by="taxid")
		}
	}
	
	# we only print out genera, families, orders, that are present in any one taxa
	total <- apply(all.rat.family.diets[,2:25], 1, sum)
	all.rat.family.diets <- all.rat.family.diets[which(total>0),]
	total.rat.print <- as.data.frame(t(c("NA",total.rat,"NA","NA")))
	colnames(total.rat.print) <- c("taxid", all.rat.locs, "taxlevel", "name")
	
	hmf <- all.rat.family.diets[,2:25]
	row.names(hmf) <- all.rat.family.diets$name
	total <- apply(all.rat.genus.diets[,2:25], 1, sum)
	all.rat.genus.diets <- all.rat.genus.diets[which(total>0),]
	
	total <- apply(all.rat.order.diets[,2:25], 1, sum)
	all.rat.order.diets <- all.rat.order.diets[which(total>0),]
		
	if(filter.fp) {
		write.table(rbind(all.rat.family.diets, total.rat.print), file="./family_numbers_filt.txt", sep="\t", quote=F, row.names=F)
		write.table(rbind(all.rat.genus.diets, total.rat.print), file="./genus_numbers_filt.txt", sep="\t", quote=F, row.names=F)
		write.table(rbind(all.rat.order.diets, total.rat.print), file="./order_numbers_filt.txt", sep="\t", quote=F, row.names=F)
		write.table(tax.hier, file="./tax_hierarchy.txt", sep="\t", quote=F, row.names=F)
		cat("\nTables, filtered data\n")
	}
	if(filter.fp==0) {
		write.table(all.rat.family.diets, file="./family_numbers_unfilt.txt", sep="\t", quote=F)
		write.table(all.rat.genus.diets, file="./genus_numbers_unfilt.txt", sep="\t", quote=F)
		write.table(all.rat.order.diets, file="./order_numbers_unfilt.txt", sep="\t", quote=F)
		cat("Tables, UNfiltered data\n")
	}
	
	
	###############################################
	###############################################
	### Figure 5 is boxplot of the number of families, etc.
	###############################################
	###############################################
	if(fig5) {
		
		pdf(file="./Fig5.pdf", height=4, width=4.5)
		par(mfrow=c(1,1))
		par(las=1)
		par(mar=c(5,4,2,2))
		par(xpd=T)
		order.per.rat <- apply(all.rat.order.diets[2:25]>0, 2, sum)
		family.per.rat <- apply(all.rat.family.diets[2:25]>0, 2, sum)
		taxa.matrix <- matrix(c(family.per.rat,order.per.rat), ncol=6, byrow=F)
		colnames(taxa.matrix) <- c("LB.fam", "OB.fam", "WP.fam", "LB.ord", "OB.ord", "WP.ord")
		cat("\nmean families:", round(mean(taxa.matrix[,1:3]),1), "min:", min(taxa.matrix[,1:3]), "max:", max(taxa.matrix[,1:3]))
		cat("\nmean orders:", round(mean(taxa.matrix[,4:6]),1), "min:", min(taxa.matrix[,4:6]), "max:", max(taxa.matrix[,4:6]))
		boxplot(taxa.matrix, names=c("LB", "OB", "WP", "LB", "OB", "WP"), col="light blue", at=c(1:3,5:7), ylab="Number of taxa per rat")
		points(sort(rep(c(1:3,5:7),8))+runif(48,-0.1,0.1), taxa.matrix, cex=1, pch=21, bg="dark grey")
		segments(0.9,-4.6,3.1,-4.6, col=1, lwd=1.5)
		segments(4.9,-4.6,7.1,-4.6, col=1, lwd=1.5)
		mtext("Family", 1, line=3, at=2)
		mtext("Order", 1, line=3, at=6)
		
		dev.off()
		
		cat("\nfigure 5\n")
	}
	
	###############################################
	###############################################
	### End fig 5
	###############################################
	###############################################	
	
	
	###############################################
	###############################################
	### Figure 6 is plot of the frequency of genera,
	### families, and orders in rat diets
	###############################################
	###############################################
	if(fig6) {
		library(gplots)
				
		hmf <- all.rat.family.diets[,2:25]
		hmo <- all.rat.order.diets[,2:25]
		row.names(hmf) <- all.rat.family.diets$name
		row.names(hmo) <- all.rat.order.diets$name
		hmf <- data.matrix(hmf)
		hmo <- data.matrix(hmo)
		
		# get presence / absence and order by that and then mean
		hmf.scale <- scale(hmf, scale=rev(total.rat), center=F)
		hmf.means <- apply(hmf.scale,1,mean)
		o <- order(apply(hmf.scale>0,1,sum), hmf.means, decreasing=T)
		hmf.scale <- hmf.scale[o,]
		
		hmo.scale <- scale(hmo, scale=rev(total.rat), center=F)
		hmo.means <- apply(hmo.scale,1,mean)
		o <- order(apply(hmo>0,1,sum), hmo.means, decreasing=T)
		hmo.scale <- hmo.scale[o,]
		
		# should we manually filter out some obvious FP taxa? Proabably not.
		# families: Salmonidae Poeciliidae Octopodidae Buthidae
		# orders: Salmoniformes Scorpiones Octopoda
		
		# append unassigned families by checking how many reads
		# there are compared to what is expected
		hmf.abs <- t(data.frame(1-apply(hmf.scale,2,sum, na.rm=T)))
		hmf.scale <- rbind(hmf.scale,hmf.abs)
		rownames(hmf.scale)[dim(hmf.scale)[1]] <- "Not assigned"
		
		phyla <- c("Arthropoda", "Mollusca", "Nematoda", "Platyhelminthes", "Chordata", "Ascomycota", "Basidiomycota", "Mucoromycota", "Streptophyta")
		color.spec <- c(463,616,8,114,132,503,499,569,139)
		
		phylum.colors <- matrix(c(phyla,color.spec),byrow=F,ncol=2)
		legend.y <- c(rep(-8.5,5),rep(-10,4))
		legend.x <- c(seq(10,45,by=8),seq(10,35,by=8))
		
		hmf.fams <- rownames(hmf.scale)
		hmf.phy <- tax.hier$Phylum.Name[match(hmf.fams, tax.hier$Family.Name)]
		hmf.col <- as.numeric(phylum.colors[match(hmf.phy,phylum.colors[,1]),2])
		
		# append unassigned orders
		hmo.abs <- t(data.frame(1-apply(hmo.scale,2,sum, na.rm=T)))
		hmo.scale <- rbind(hmo.scale,hmo.abs)
		rownames(hmo.scale)[dim(hmo.scale)[1]] <- "Not assigned"
		
		hmo.ords <- rownames(hmo.scale)
		hmo.phy <- tax.hier$Phylum.Name[match(hmo.ords, tax.hier$Order.Name)]
		hmo.col <- as.numeric(phylum.colors[match(hmo.phy,phylum.colors[,1]),2])
		
		# order cols by name
		oc <- order(colnames(hmf.scale))
		hmf.scale <- t(hmf.scale[,oc])
		hmf.scale[hmf.scale==0] <- NA
		
		oc <- order(colnames(hmo.scale))
		hmo.scale <- t(hmo.scale[,oc])
		hmo.scale[hmo.scale==0] <- NA
		
		# a pain to reorder the names, do so by hand here
		ordered.names <- sort(c("OB2","OB13","OB7","OB5","LB7","LB6","LB3","LB5","WP11","WP9","WP5","WP8","WP10","OB10","LB1","WP12","OB11","LB2","WP2","OB12","LB4","WP3","OB14","LB8"))
		ordered.names <- ordered.names[c(17:24,9:16,1:8)]
		
		num.fam <- dim(hmf.scale)[2]
		num.ord <- dim(hmo.scale)[2]
		num.fam <- 56
		num.ord <- 56
		
		heat.breaks <- c(0,0.0001,0.0002,0.01,0.02,0.05,0.1,0.2,0.4,0.8,1)
		legend.breaks <- c("< 0.01","0.01-0.02","0.02-0.05","0.05-0.1","0.1-0.2","0.2-0.4","0.4-0.8","> 0.8")

		pdf(file="./Fig6a.pdf", height=6, width=13)
		par(xpd=T)
		#draw heatmap
		heatmap.2(hmf.scale[,c(1:(num.fam-1),dim(hmf.scale)[2])], Rowv=F, Colv=F, dendrogram="none", labRow=ordered.names, key=F, trace="none", colsep=1:num.fam, rowsep=1:24, sepcolor="white", sepwidth=c(0.03,0.03), col=colorRampPalette(brewer.pal(8,"YlGnBu")), na.color="grey95", margins=c(12.,12), breaks=heat.breaks, offsetRow=-0.2, offsetCol=0.8, srtCol=70, lhei=c(0.02,0.98), lwid=c(0.02,0.98), cexCol=0.9, add.expr=list(
		segments(x0=c(0,0), y0=c(8.46,16.46), x1=rep(num.fam,2), lwd=2), # dividing lines
		points(1:(num.fam-1), rep(0.02,(num.fam-1)), pch=21, bg=colors()[hmf.col], cex=2.1), # phylum points
		points(legend.x,legend.y, pch=21, bg=colors()[color.spec], cex=2.1), # phylum legend
		text(legend.x, legend.y, labels=phyla, pos=4, offset=0.6, cex=0.9), # phylum legend text
		rect(rep(num.fam+4.5,8),9:16,rep(num.fam+5.5,8),10:17, col=colorRampPalette(brewer.pal(8,"YlGnBu"))(9)[2:9] , border="white"), # legend
		text(rep(num.fam+5.3,8), seq(9.5,16.5,by=1), labels=legend.breaks, pos=4) ) ) # legend text
		
		dev.off()
		
		pdf(file="./Fig6b.pdf", height=6., width=13)
		par(xpd=T)
		#draw heatmap
		heatmap.2(hmo.scale[,c(1:(num.ord-1),dim(hmo.scale)[2])], Rowv=F, Colv=F, dendrogram="none", labRow=ordered.names, key=F, trace="none", colsep=1:num.ord, rowsep=1:24, sepcolor="white", sepwidth=c(0.03,0.03), col=colorRampPalette(brewer.pal(8,"YlGnBu")), na.color="grey95", margins=c(12,12), breaks=heat.breaks, offsetRow=-0.2, offsetCol=.8, srtCol=70, lhei=c(0.02,0.98), lwid=c(0.02,0.98), cexCol=0.9, add.expr=list(
		segments(x0=c(0,0), y0=c(8.46,16.46), x1=rep(num.ord,2), lwd=2), # dividing lines
		points(1:(num.ord-1), rep(0.02,(num.ord-1)), pch=21, bg=colors()[hmo.col], cex=2.1), # phylum points
		points(legend.x,legend.y, pch=21, bg=colors()[color.spec], cex=2.1), # phylum legend
		text(legend.x, legend.y, labels=phyla, pos=4, offset=0.6, cex=0.9), # phylum legend text
		rect(rep(num.ord+4.5,8),9:16,rep(num.ord+5.5,8),10:17, col=colorRampPalette(brewer.pal(8,"YlGnBu"))(9)[2:9], border="white"),  # legend
		text(rep(num.ord+5.3,8), seq(9.5,16.5,by=1), labels=legend.breaks, pos=4) ) ) # legend text
		
		dev.off()

		cat("figure 6 \n")

	}
	###############################################
	###############################################
	### End fig 6
	###############################################
	###############################################	
	

}
###############################################
###############################################
### End of table making
###############################################
###############################################

###############################################
###############################################
### Figure 7 is a plot of the CAP analysis
### 
###############################################
###############################################

if(fig7) {
	library(shape)
	legend.cex <- 1.4
	pdf(file="./Fig7.pdf", height=3.5, width=10)
	par(mfrow=c(1,3))
	par(las=1)
	par(mar=c(5,5,4,2))
	
	mds.fam <- read.table(file="./fig_data/filter_family_nMDS.txt", header=T)
	mds.fam <- merge(mds.fam, all.rat.loc.lookup, by="Sample")
	
	cap.fam <- read.table(file="./fig_data/filter_family_CAP.txt", header=T)
	cap.fam <- merge(cap.fam, all.rat.loc.lookup, by="Sample")
	
	cap.pearson <- read.table(file="./fig_data/CAP_pearson.txt", header=T)
	
	### order: LB, OB, WP
	o <- order(mds.fam$paper.names)
	mds.fam.col <- c(rep("light blue",8),rep("orange",8),rep("purple",8))

	ylim <- c(-1.6,1.6)
	plot(mds.fam$MDS1, mds.fam$MDS2, bg=mds.fam.col, cex=1.1, pch=23, xlim=c(-1.4,1.8), ylim=ylim, xlab="nMDS1", ylab="nMDS2")
	### Takes some manual adjustment :( :(
	text(mds.fam$MDS1, mds.fam$MDS2, pos=4, offset=0.3, labels=mds.fam$paper.names, cex=0.9)
	mtext("a", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)
	
	o <- order(cap.fam$paper.names)
	cap.fam.col <- c(rep("light blue",8),rep("orange",8),rep("purple",8))
	ylim <- c(-0.5,0.3)
	plot(cap.fam$CAP1, cap.fam$CAP2, bg=cap.fam.col, cex=1.1, pch=23, xlim=c(-0.4,0.5), ylim=ylim, xlab="CAP1", ylab="CAP2")
	
	### Takes some manual adjustment :( :(
	text(cap.fam$CAP1, cap.fam$CAP2, pos=4, offset=0.3, labels=cap.fam$paper.names, cex=0.9)
	mtext("b", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)
	
	ylim <- c(-0.5,0.3)
	arrow.scale <- 2.3
	plot(-1,-1, xlim=c(-0.4,0.5), ylim=ylim, xlab="CAP1", ylab="CAP2")
	num.taxa <- 8
	o <- order(apply(abs(cap.pearson[2:3]),1,sum), decreasing=T)
	cap.pearson <- cap.pearson[o,]
	Arrows(rep(0,num.taxa),rep(0,num.taxa), cap.pearson$CAP1[1:num.taxa]/arrow.scale, cap.pearson$CAP2[1:num.taxa]/arrow.scale, col="dark blue", lwd=1., arr.type="triangle", arr.length=0.12, arr.width=0.12)
	
	### Takes some manual adjustment :( :(
	text(cap.pearson$CAP1[1:num.taxa]/arrow.scale*1.15, cap.pearson$CAP2[1:num.taxa]/arrow.scale*1.15, offset=0.3, labels=cap.pearson$Taxon[1:num.taxa], cex=0.9)
	mtext("c", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=4)
	
	dev.off()
	cat("Figure 7\n")
}

###############################################
###############################################
### End Fig 7 plotting
###############################################
###############################################

###############################################
###############################################
### Figure 8 is a plot of genera with fraction 
### of sequences that are genomic vs. rrna/mtdna/chloro/microsatellite 
### rrna/mtdna/chloro/microsatellite
###############################################
###############################################
if(fig8) {
	legend.cex <- 2
	pdf(file="./Fig8.pdf", height=9.5, width=5.5)
	par(mfrow=c(1,1))
	par(mar=c(5,4.5,5,1.5))
	par(las=1)
	
	genus.num <- as.data.frame(table(hq.diet$Genus))
	genus.freq <- genus.num[which(genus.num$Freq > 20),]
	genus.freq <- rbind(genus.freq, c("Rattus",dim(subset(hq, Genus=="Rattus"))[1]))
	non.gen <- vector()
	genomic <- vector()
	for (i in 1:dim(genus.freq)[1]) {
		g <- subset(hq, Genus==genus.freq[i,1])
		non.genom.subset <- g[grep("rrna|mtDNA|microsat|chloro|rdna|mitoch|mrna|ribosom|subunit|cytochrom|complete cds|exon|transcript|ncrna", g$subject.description,ignore.case=T),]
		non.gen[i] <- 1 - dim(non.genom.subset)[1]/dim(g)[1]
		
		genom.subset <- g[grep("rrna|mtDNA|microsat|chloro|rdna|mitoch|mrna|ribosom|subunit|cytochrom|complete cds|exon|transcript|ncrna", g$subject.description,ignore.case=T, invert=T),]
		genomic[i] <- length(grep("genom|BAC.+complete|chromosom|clone|fosmid", genom.subset$subject.description, perl=T))/dim(g)[1]
	}
	
	genus.freq <- cbind(genus.freq,non.gen)
	genus.freq <- cbind(genus.freq,genomic)
	colnames(genus.freq) <- c("genus","total","non.genomic", "genomic")
	o <- order(genus.freq$non.genomic, decreasing=T)	
	genus.freq <- genus.freq[o,]
	plot(genus.freq$non.genomic,1:dim(genus.freq)[1], xlab="", pch=21, bg="light blue", xlim=c(0,1.4), xaxt="n", yaxt="n", cex=1.3, ylab="", bty="n")
	text(genus.freq$non.genomic,1:dim(genus.freq)[1], labels=paste(genus.freq$genus,genus.freq$total,sep=" "), cex=0.9, pos=4, offset=0.45, font=3)
	axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
	axis(3, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
	mtext("Fraction of reads matching\ngenomic sequence",side=1,line=3.5, at=0.5)
	abline(v=c(0,0.2,0.4,0.6,0.8,1), lty=2, col="grey", lwd=1)

	dev.off()
	cat("Figure 8\n")
}	
###############################################
###############################################
### End Fig8 plotting
###############################################
###############################################


###############################################
###############################################
### Figure S1 are biplots of qual/length for each barcode
###  this is now abbreviated ot mouse and rat only
###############################################
###############################################
if(figS1) {
	library(MASS)
	library(RColorBrewer)

	# open PDF
	pdf(file="./FigS1.pdf", height=4, width=4)
	par(mfrow=c(4,13))
	par(mar=c(4,5,2,2))
	par(las=1)
	run.dates <- c(rep("jan",13),rep("mar",13))
	#layout(matrix(1:52, nrow=13, ncol=4, byrow=F))
	for (i in 1:length(qual.files)) {
		cat(qual.files[i], "\n", run.dates[i], "\n\n")
		stats <- read.table(qual.files[i], sep="\t")
		tot.bp <- round(sum(stats[,2])/1e6,2)
		tot.reads <- length(stats[,2])
		k10 <- length(which(stats[,2]>1e4))
		k50 <- length(which(stats[,2]>5e4))
		hl <- hist(log10(stats[,2]), breaks=seq(0,7,by=0.02), plot=F)
		stats[,3] <- 10*(-log10((1-stats[,3])))
		hist.qual <- hist(stats[,3], breaks=seq(0,60,by=0.08), plot=F)
		
		
		nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(1,0.3), c(0.3,1), TRUE)
		k <- kde2d(log10(stats[,2]), stats[,3], n=100, lims=c(1,4,2,15))
		#crp.fun <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
		crp.fun <- colorRampPalette((brewer.pal(9, 'BuGn')))
		crp <- crp.fun(100)
		par(mar=c(5,5,0,0))
		image(k, col=crp, xlim=c(2.,4), ylim=c(3,14), xaxt="n", xlab="Length", ylab="Average quality score")
		
		# add some data about the fastq read numbers, etc.
		#text(3.5,y.text,labels=c(qual.names[i],paste("Mbp", tot.bp, sep=" "), paste("Reads",tot.reads,sep=" "), paste(">10K",k20,sep=" "), paste(">50K",k50,sep=" ")), pos=1, cex=1)
		
		# put on a reasonable axis
		axis(1, at=2:4, labels=c("100", "1K", "10K"))
		
		### plot length
		par(mar=c(0,5,2,0))
		plot(-1,-1, xlim=c(2,4), ylim=c(0,0.075), xaxt="n", bty="n", yaxt="n", xlab="", ylab="", main=paste(run.dates[i], qual.names[i], " Mbp",tot.bp, "  Reads",tot.reads, sep=" "), cex.main=0.8, font.main=1)
		
		polygon(c(sort(rep(hl$breaks,2)), hl$breaks[1]), c(0,hl$counts[sort(rep(1:length(hl$counts),2))]/sum(hl$counts),0,0), col="light blue")
		
		### plot quality
		par(mar=c(5,0,0,2))
		plot(-1,-10, xlim=c(0,0.07), ylim=c(3,14), xaxt="n", bty="n", yaxt="n", xlab="", ylab="")
		polygon(c(0, hist.qual$counts[sort(rep(1:length(hist.qual$counts),2))]/sum(hist.qual$counts),0,0), c(sort(rep(hist.qual$breaks,2)), hist.qual$breaks[1]), col="grey")

	}
	dev.off()
	cat("Figure S1\n")

}
###############################################
###############################################
### END FigS1 plotting
###############################################
###############################################

###############################################
###############################################
### FigS2 plotting
### Read quality vs identity and alignment length and ratio
### this is no abridged
###############################################
###############################################
if(figS2) {
	pdf(file="./FigS2_3.pdf", height=4, width=8)
	par(mfrow=c(1,2))
	par(las=1)
	legend.cex <- 1.4
	window <- 801
	# what?
	od <- order(hq.diet$quals_mean)
	or <- order(hq.rat$quals_mean)
	om <- order(hq.mouse$quals_mean)
	qd <- hq.diet$quals_mean[od]
	qr <- hq.rat$quals_mean[or]
	qm <- hq.mouse$quals_mean[om]
	idd <- hq$X..identity[od]
	ld <- hq.diet$alignment.length[od]
	rd <- hq.diet$Ratio[od]
	idr <- hq.rat$X..identity[or]
	lr <- hq.rat$alignment.length[or]
	rr <- hq.rat$Ratio[or]
	idm <- hq.mouse$X..identity[om]
	lm <- hq.mouse$alignment.length[om]
	rm <- hq.mouse$Ratio[om]
	
	iddrm <- runmed(idd,window)
	ldrm <- runmed(ld,window)
	rdrm <- runmed(rd,window)
	
	idrrm <- runmed(idr,window)
	lrrm <- runmed(lr,window)
	rrrm <- runmed(rr,window)
	
	idmrm <- runmed(idm,window)
	lmrm <- runmed(lm,window)
	rmrm <- runmed(rm,window)
	
	
	ylim <- c(70,98)
	xlim=c(0.8,0.95)
	if(0) {
	plot(qd, idd, xlab="Mean read accuracy", ylab="Percent identity of top BLAST hit", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qd, iddrm, col="orange", lwd=2)
	mtext("a", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	abline(0,100,lty=2,col="red", lwd=2)
	}
	
	plot(qr, idr, xlab="Mean read accuracy", ylab="Percent identity of top BLAST hit", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qr, idrrm, col="orange", lwd=2)
	mtext("a", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	abline(0,100,lty=2,col="red", lwd=2)
	
	plot(qm, idm, xlab="Mean read accuracy", ylab="Percent identity of top BLAST hit", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qm, idmrm, col="orange", lwd=2)
	mtext("b", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	abline(0,100,lty=2,col="red", lwd=2)
	
	if(0) {
	ylim <- c(0.1,2)
	plot(qd, ld/1e3, xlab="Mean read accuracy", ylab="Alignment length (Kbp)", pch=19, cex=0.5, col=rgb(0,0,0,0.4), log="y", xlim=xlim, ylim=ylim)
	lines(qd, ldrm/1e3, col="orange", lwd=2)
	mtext("d", side=2, cex=legend.cex, at=10^(log10(ylim[1]) + (log10(ylim[2])-log10(ylim[1]))*1.12), adj=3)
	
	plot(qr, lr/1e3, xlab="Mean read accuracy", ylab="Alignment length (Kbp)", pch=19, cex=0.5, col=rgb(0,0,0,0.4), log="y", xlim=xlim, ylim=ylim)
	lines(qr, lrrm/1e3, col="orange", lwd=2)
	mtext("e", side=2, cex=legend.cex, at=10^(log10(ylim[1]) + (log10(ylim[2])-log10(ylim[1]))*1.12), adj=3)

	plot(qm, lm/1e3, xlab="Mean read accuracy", ylab="Alignment length (Kbp)", pch=19, cex=0.5, col=rgb(0,0,0,0.4), log="y", xlim=xlim, ylim=ylim)
	lines(qm, lmrm/1e3, col="orange", lwd=2)
	mtext("f", side=2, cex=legend.cex, at=10^(log10(ylim[1]) + (log10(ylim[2])-log10(ylim[1]))*1.12), adj=3)

	ylim <- c(0,1.1)
	plot(qd, rd, xlab="Mean read accuracy", ylab="Ratio of read length to alignment length", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qd, rdrm, col="orange", lwd=2)
	mtext("g", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	
	ylim <- c(0,1.1)
	plot(qr, rr, xlab="Mean read accuracy", ylab="Ratio of read length to alignment length", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qr, rrrm, col="orange", lwd=2)
	mtext("h", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	
	ylim <- c(0,1.1)
	plot(qm, rm, xlab="Mean read accuracy", ylab="Ratio of read length to alignment length", pch=19, cex=0.5, col=rgb(0,0,0,0.4), xlim=xlim, ylim=ylim)
	lines(qm, rmrm, col="orange", lwd=2)
	mtext("i", side=2, cex=legend.cex, at=(ylim[1] + (ylim[2]-ylim[1])*1.12), adj=3)
	}
	cat("Fig S2\n")

dev.off()
}
###############################################
###############################################
### END FigS2 plotting
###############################################
###############################################


###############################################
###############################################
### Tables S1 and S2 are read numbers and bp for
### each run
###############################################
###############################################

###############################################
###############################################
### End of S1 and S2 construction
### 
###############################################
###############################################


###############################################
###############################################
### Datafile is a table of all blast hits and taxa
###############################################
###############################################
if(print.data) {
	

		hq.edit <- hq[,c(4:25,28:33)]
		names <- c("seq.name","subject.id","percent.ident","align.length","mismatches","gap.opens","query.start","query.end","subject.start","subject.end","evalue","bit.score","subject.descr","rat.id","barcode","species","blast.genus","read.length","mean.read.accur","median.read.accur","sequence.run","readl.alignl.ratio","megan.kingdom","megan.phylum","megan.class","megan.order","megan.family","megan.genus")
		colnames(hq.edit) <- names
		hq.edit <- hq.edit[,c(1,14,15,18,19,20,21,2,3,4,5,6,7,8,9,10,12,13,16,17,22:28)]
		write.table(hq.edit, file="./datafile_S1.tab", sep="\t",row.names=F)
	
	if(filter.fp) {
		hq.edit <- hq.diet[,c(4:25,28:33)]
		names <- c("seq.name","subject.id","percent.ident","align.length","mismatches","gap.opens","query.start","query.end","subject.start","subject.end","evalue","bit.score","subject.descr","rat.id","barcode","species","blast.genus","read.length","mean.read.accur","median.read.accur","sequence.run","readl.alignl.ratio","megan.kingdom","megan.phylum","megan.class","megan.order","megan.family","megan.genus")
		colnames(hq.edit) <- names
		hq.edit <- hq.edit[,c(1,14,15,18,19,20,21,2,3,4,5,6,7,8,9,10,12,13,16,17,22:28)]
		write.table(hq.edit, file="./datafile_S2.tab", sep="\t",row.names=F)
	}
	
}

###############################################
###############################################
### End of S1 and S2 construction
### 
###############################################
###############################################
