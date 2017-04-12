#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts, "i=s", "o=s", "l=i", "p=f", "pie=s", "bar=s", "mb=i", "ml=i", "h!");

my $usage = <<"USAGE";
	Program: $0
	Usage: perl $0 [options]
		-i	GO.list.level*.txt/level2.go.txt
		-o	output R script prefix, defaults: rcmd
		-l	GO level, 2/3/4, defaults: 2
		-p	if num/all <= p push to others, defaults: 0.01
		-bar	[T/F] plot bar or not, defaults: T
		-pie	[T/F] plot pie or not, defaults: T
		-mb	bottom space for bar plot, defaults: 12
		-ml	left space for bar plot, defaults: 6
USAGE

die $usage if ( !( $opts{i} ) || $opts{h} );

#define defaults
$opts{o}=$opts{o}?$opts{o}:"rcmd";
$opts{bar}=$opts{bar}?$opts{bar}:"T";
$opts{pie}=$opts{pie}?$opts{pie}:"T";
$opts{l}=$opts{l}?$opts{l}:2;
$opts{mb}=$opts{mb}?$opts{mb}:12;
$opts{ml}=$opts{ml}?$opts{ml}:6;
$opts{p}=$opts{p}?$opts{p}:0.01;

### plot pie/bar
open RCMD, ">$opts{o}.r" || die $!;
print RCMD "

library(\"RColorBrewer\")
ColorList<-brewer.pal(11,\"Paired\")
gostat_level2<-read.delim(\"$opts{i}\",header=T,check.names=F,skip=1,stringsAsFactors=F)
p<-$opts{p}

##BP level2
BP_2<-gostat_level2[which(gostat_level2[[1]]==\"biological_process\"),]
BP_2_sort<-BP_2[order(BP_2[[3]],decreasing=T),][,2:3]
percent<-matrix(nrow=nrow(BP_2_sort), ncol=1)
for (i in 1:nrow(BP_2_sort)){
	percent[i,1]<-BP_2_sort[i,2]/sum(BP_2_sort[,2])
}
BP_2_sort<-cbind(BP_2_sort, percent)

if (length(BP_2_sort[BP_2_sort<p])>1) {
	BP_2_sort_over_p<-BP_2_sort[which(BP_2_sort\$percent>p),]
	BP_2_sort_below_p<-BP_2_sort[which(BP_2_sort\$percent<=p),]
	others<-c(\"others\",sum(BP_2_sort_below_p[,2]),sum(BP_2_sort_below_p[,3]))
	BP_2_sort<-rbind(BP_2_sort_over_p,others)
} else {
	BP_2_sort<-BP_2_sort
}

BP_2_splice<-as.numeric(as.character(BP_2_sort[[2]]))
#BP_2_label<-paste(100*as.numeric(format(as.numeric(BP_2_sort[[3]]),digits=2)),\"%\",sep=\"\")
BP_2_label<-as.character(BP_2_sort[[2]])
BP_2_number<-as.character(BP_2_sort[[2]])
BP_2_legend<-paste(as.character(BP_2_sort[[1]]),\" (\",BP_2_number,\")\",sep=\"\")

#plot BP level2 pie
if (T == T) {
	pdf0 <- paste(\"Biological_process.level\", 2, sep=\"\")
	pdf <- paste(pdf0, \".pie.pdf\", sep=\"\")
	pdf(pdf, w=12, h=8)
	par(mfrow=c(1,2))
	pie(BP_2_splice,labels=BP_2_label,col=ColorList,clockwise=T,cex=1,radius=0.85,border=\"white\")
	legend(\"bottom\",\"Biological Process\",bty=\"n\",cex=1.5,text.font=2)
	
	#plot BP level2 legend
	plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
	legend(\"left\",legend=BP_2_legend,fill=ColorList,cex=1,bty=\"n\")
	dev.off()
}

#plot BP level2 bar
if ($opts{bar} == T) {
	BP_2_sort<-BP_2[order(BP_2[[3]],decreasing=T),][,2:3]
	
	pdf0 <- paste(\"Biological_process.level\", $opts{l}, sep=\"\")
	pdf <- paste(pdf0, \".bar.pdf\", sep=\"\")
	pdf(pdf, w=10, h=7)
	par(mar=c($opts{mb},$opts{ml},6,3))
	bp <- barplot(BP_2_sort[,2], col=\"#4F81BD\",border=FALSE, axisnames=F,las=0,width=1,space=0.5, plot=T,beside=TRUE,main=\"Biological process\",cex.main=1.2,ylab=\"Number of proteins\",cex.axis=0.8,las=0)
	text(x=bp,y=-0.015*max(BP_2_sort[,2]),srt=45,adj=1,labels=BP_2_sort[,1],xpd=T,cex=0.8,cex.axis=0.8)
	jump<-max(BP_2_sort[,2])/10
	text(bp,(BP_2_sort[,2]+jump/3),labels=BP_2_sort[,2],xpd=T,cex=0.8)

	dev.off()
}


##CC level2
CC_2<-gostat_level2[which(gostat_level2[[1]]==\"cellular_component\"),]
CC_2_sort<-CC_2[order(CC_2[[3]],decreasing=T),][,2:3]
percent<-matrix(nrow=nrow(CC_2_sort), ncol=1)
for (i in 1:nrow(CC_2_sort)){
	percent[i,1]<-CC_2_sort[i,2]/sum(CC_2_sort[,2])
}
CC_2_sort<-cbind(CC_2_sort, percent)

if (length(CC_2_sort[CC_2_sort<p])>1) {
	CC_2_sort_over_p<-CC_2_sort[which(CC_2_sort\$percent>p),]
	CC_2_sort_below_p<-CC_2_sort[which(CC_2_sort\$percent<=p),]
	others<-c(\"others\",sum(CC_2_sort_below_p[,2]),sum(CC_2_sort_below_p[,3]))
	CC_2_sort<-rbind(CC_2_sort_over_p,others)
} else {
	CC_2_sort<-CC_2_sort
}

CC_2_splice<-as.numeric(as.character(CC_2_sort[[2]]))
#CC_2_label<-paste(100*as.numeric(format(as.numeric(CC_2_sort[[3]]),digits=2)),\"%\",sep=\"\")
CC_2_label<-as.character(CC_2_sort[[2]])
CC_2_number<-as.character(CC_2_sort[[2]])
CC_2_legend<-paste(as.character(CC_2_sort[[1]]),\" (\",CC_2_number,\")\",sep=\"\")

#plot CC level2 pie
if (T == T) {
	pdf0 <- paste(\"Cellular_component.level\", 2, sep=\"\")
	pdf <- paste(pdf0, \".pie.pdf\", sep=\"\")
	pdf(pdf, w=12, h=8)
	par(mfrow=c(1,2))
	pie(CC_2_splice,labels=CC_2_label,col=ColorList,clockwise=T,cex=1,radius=0.85,border=\"white\")
	legend(\"bottom\",\"Cellular Component\",bty=\"n\",cex=1.5,text.font=2)
	
	#plot CC level2 legend
	plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
	legend(\"left\",legend=CC_2_legend,fill=ColorList,cex=1,bty=\"n\")
	dev.off()
}

#plot CC level2 bar
if ($opts{bar} == T) {
	CC_2_sort<-CC_2[order(CC_2[[3]],decreasing=T),][,2:3]
	
	pdf0 <- paste(\"Cellular_component.level\", $opts{l}, sep=\"\")
	pdf <- paste(pdf0, \".bar.pdf\", sep=\"\")
	pdf(pdf, w=10, h=7)
	par(mar=c($opts{mb},$opts{ml},6,3))
	bp <- barplot(CC_2_sort[,2], col=\"#4F81BD\",border=FALSE, axisnames=F,las=0,width=1,space=0.5, plot=T,beside=TRUE,main=\"Biological process\",cex.main=1.2,ylab=\"Number of proteins\",cex.axis=0.8,las=0)
	text(x=bp,y=-0.015*max(CC_2_sort[,2]),srt=45,adj=1,labels=CC_2_sort[,1],xpd=T,cex=0.8,cex.axis=0.8)
	jump<-max(CC_2_sort[,2])/10
	text(bp,(CC_2_sort[,2]+jump/3),labels=CC_2_sort[,2],xpd=T,cex=0.8)

	dev.off()
}



##MF level2
MF_2<-gostat_level2[which(gostat_level2[[1]]==\"molecular_function\"),]
MF_2_sort<-MF_2[order(MF_2[[3]],decreasing=T),][,2:3]
percent<-matrix(nrow=nrow(MF_2_sort), ncol=1)
for (i in 1:nrow(MF_2_sort)){
	percent[i,1]<-MF_2_sort[i,2]/sum(MF_2_sort[,2])
}
MF_2_sort<-cbind(MF_2_sort, percent)

if (length(MF_2_sort[MF_2_sort<p])>1) {
	MF_2_sort_over_p<-MF_2_sort[which(MF_2_sort\$percent>p),]
	MF_2_sort_below_p<-MF_2_sort[which(MF_2_sort\$percent<=p),]
	others<-c(\"others\",sum(MF_2_sort_below_p[,2]),sum(MF_2_sort_below_p[,3]))
	MF_2_sort<-rbind(MF_2_sort_over_p,others)
} else {
	MF_2_sort<-MF_2_sort
}

MF_2_splice<-as.numeric(as.character(MF_2_sort[[2]]))
#MF_2_label<-paste(100*as.numeric(format(as.numeric(MF_2_sort[[3]]),digits=2)),\"%\",sep=\"\")
MF_2_label<-as.character(MF_2_sort[[2]])
MF_2_number<-as.character(MF_2_sort[[2]])
MF_2_legend<-paste(as.character(MF_2_sort[[1]]),\" (\",MF_2_number,\")\",sep=\"\")

#plot CC level2 pie
if (T == T) {
	pdf0 <- paste(\"Molecular_function.level\", 2, sep=\"\")
	pdf <- paste(pdf0, \".pie.pdf\", sep=\"\")
	pdf(pdf, w=12, h=8)
	par(mfrow=c(1,2))
	pie(MF_2_splice,labels=MF_2_label,col=ColorList,clockwise=T,cex=1,radius=0.85,border=\"white\")
	legend(\"bottom\",\"Molecular Function\",bty=\"n\",cex=1.5,text.font=2)
	
	#plot CC level2 legend
	plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
	legend(\"left\",legend=MF_2_legend,fill=ColorList,cex=1,bty=\"n\")
	dev.off()
}

#plot MF level2 bar
if ($opts{bar} == T) {
	MF_2_sort<-MF_2[order(MF_2[[3]],decreasing=T),][,2:3]
	
	pdf0 <- paste(\"Molecular_function.level\", $opts{l}, sep=\"\")
	pdf <- paste(pdf0, \".bar.pdf\", sep=\"\")
	pdf(pdf, w=10, h=7)
	par(mar=c($opts{mb},$opts{ml},6,3))
	bp <- barplot(MF_2_sort[,2], col=\"#4F81BD\",border=FALSE, axisnames=F,las=0,width=1,space=0.5, plot=T,beside=TRUE,main=\"Biological process\",cex.main=1.2,ylab=\"Number of proteins\",cex.axis=0.8,las=0)
	text(x=bp,y=-0.015*max(MF_2_sort[,2]),srt=45,adj=1,labels=MF_2_sort[,1],xpd=T,cex=0.8,cex.axis=0.8)
	jump<-max(MF_2_sort[,2])/10
	text(bp,(MF_2_sort[,2]+jump/3),labels=MF_2_sort[,2],xpd=T,cex=0.8)

	dev.off()
}

";

`Rscript $opts{o}.r`;
