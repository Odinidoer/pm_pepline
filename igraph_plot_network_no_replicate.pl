#!/usr/bin/perl -w
use warnings;
use strict;

use Getopt::Long;
my %opts;
my $VERSION = "v1.20161219"; 

GetOptions( \%opts,"iLinks=s","iExp=s","n=s","o=s","s=s","h");
my $usage = <<"USAGE";
	Program:	$0   
	Discription:	plot network for different expressed proteins from STRING database (sample with no replicate)
	Version:	$VERSION
	Contact:	ting.kuang\@majorbio.com
	Usage:	perl $0 [options]		
		-iLinks*	input links file, e.g. *_vs_*.network.txt 
			###############################################
			#Protein1	Protein2	Combined_score
			#Q61315	A0A0A0MQ73	900
			#Q02248	Q61315	999
			#Q02248	Q3USK2	872
			###############################################
		-iExp*	input different express file, e.g. *_vs_*.DE.exp.xls (no head)
			##################################################################
			#Q61315	3.90E+09	6.20E+09	1.575655114	0.655951786	up
			#Q02248	3.90E+06	2.30E+06	0.605968858	-0.722684442	down
			#Q3USK2	1.80E+07	3.00E+07	1.640144665	0.71382307	up
			##################################################################
		-o*	output file prefix, eg. *_vs_*
		-n	whether label edges, T/F, default: T
		-s	whether remove the single node, T/F, default: T
		-h	show this information
	Example:	$0 -iLinks A_vs_B.network.txt -iExp A_vs_B.DE.exp.xls -o A_vs_B -n T -s T

USAGE

die $usage if((!($opts{iLinks} && $opts{iExp} && $opts{o})) || $opts{h});

$opts{n}=defined $opts{n}?$opts{n}:"T";
$opts{s}=defined $opts{s}?$opts{s}:"T";

open CMD,">$opts{o}.network.r";

print CMD "
library(\'igraph\')

links <- read.table(\"$opts{iLinks}\", header=T,as.is=T, sep = \"\\t\")
nodes <- read.table(\"$opts{iExp}\", header=F, as.is=T, sep = \"\\t\", check.names=F)
head(links)
head(nodes)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)

#根据FC值的4分位设置点的大小,大小分别为：3,4,5,6,7,9
count <- nodes\$V5
count <- abs(count)
names(count) <- nodes\$V1

gene_sep <- quantile(count, probs = seq(0, 1, by = 0.2))
gene_size <- count
base_size <- 3
gene_size[which(count < gene_sep[2])] <- base_size
gene_size[which(count < gene_sep[3] & count >= gene_sep[2])] <- base_size + 1
gene_size[which(count < gene_sep[4] & count >= gene_sep[3])] <- base_size + 2
gene_size[which(count < gene_sep[5] & count >= gene_sep[4])] <- base_size + 3
gene_size[which(count < gene_sep[6] & count >= gene_sep[5])] <- base_size + 4
gene_size[which(count >= gene_sep[6])] <- base_size + 5
V(net)\$size <- gene_size
V(net)\$size

##另选：根据FC值设置点的大小,直接设置，会出现有的点特别大，有的点特别小
#deg <- V(net)\$V4
#V(net)\$size <- deg*4
#V(net)\$size

##另选：根据连线多少设置点的大小，会出现有的点特别大，有的点特别小
#deg <- degree(net, mode=\"all\")
#V(net)\$size <- deg+4
#V(net)\$size

#根据combined_score设置连线的宽度，0, 0.5, 1, 1.5, 2, 2.5
egde_wide <- links\$Combined_score
names(egde_wide) <- seq(1:length(egde_wide))
egde_wide[which(egde_wide < 500)] <- 0
egde_wide[which(egde_wide < 600 & egde_wide >= 500)] <- 0.5
egde_wide[which(egde_wide < 700 & egde_wide >= 600)] <- 1
egde_wide[which(egde_wide < 800 & egde_wide >= 700)] <- 1.5
egde_wide[which(egde_wide < 900 & egde_wide >= 800)] <- 2
egde_wide[which(egde_wide < 1000 & egde_wide >= 900)] <- 2.5
E(net)\$width <- egde_wide

#根据上下调设置点的颜色
colrs <- c(\"green\", \"red\")
aa <- V(net)\$V6
aa[which(aa==\"up\")] <- 2
aa[which(aa==\"down\")] <- 1
aa
V(net)\$color <- colrs[as.numeric(aa)]
V(net)\$color

#处理孤立的点,选择：去掉或保留
single <- $opts{s}
if (single == T) {
	dg <- degree(net, mode=\"all\")
	net <- induced.subgraph(net, which(dg>0))
} else if (single == F) {
	net <- net
}


#判断是否加上蛋白的名字
lab <- $opts{n}
edge_label <- NA
if(lab == T) {
	edge_label <- E(net)\$label
} else if(lab == F) {
	edge_label <- NA
}

#根据点的多少设置图图片大小
summary(net)
node_num <- length(V(net))
if (node_num < 20) {
  h <- 8
  w <- 8
} else if (node_num >= 20 & node_num < 50) {
  h <- 12
  w <- 12
} else if (node_num >= 50 & node_num < 100) {
  h <- 20
  w <- 20
} else if (node_num >= 100 & node_num < 200) {
  h <- 30
  w <- 30
} else if (node_num >= 200 & node_num < 300) {
  h <- 40
  w <- 40
} else if (node_num >= 300) {
  h <- 60
  w <- 60
}

outname <- paste(\"$opts{o}\", \"network.pdf\", sep = \".\")
pdf(outname, h=h, w=w)
par(mar=c(2,2,4,2))

#绘图
plot(net, layout = layout_with_fr,
     edge.color=\"black\", 
     edge.width=E(net)\$width,
     edge.arrow.mode=0,
     vertex.label=edge_label,
     #edge.arrow.size=.2,
     
     vertex.size=V(net)\$size, 
     vertex.frame.color=\"gray60\", 
     vertex.label.cex=1,
     vertex.label.font=1,
     
     #main=\"$opts{o}\"
	 )

dev.off()
";

close CMD;

`Rscript $opts{o}.network.r`;
