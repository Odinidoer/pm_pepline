#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Thread::Semaphore;
use Cwd 'abs_path';
use threads::shared;
use File::Spec;
use File::Spec::Functions qw(rel2abs);


my %opts;
my $VERSION = "2013-05-30";
GetOptions (\%opts,"i=s","log=i","o=s","k=i","w=f","h=f","cls=f","rls=f","ml=s","mr=s","clt=s");
my $usage = <<"USAGE";
Program : $0
Contact : $VERSION
Author : yan.wang\@majorbio.com
Description: ploting the heatmap and the subcluster line
Usage:perl $0 [options]	     
Options: 
         -i   *         STRING  matrix for heatmap plot
	                        ################################
			        # gene_id sample1 sample2 ...
                                # gene1 12 24
                                # gene2 45 7
                                # ...
         -o   *          STRING output file prefix,default: heatmap
         -log            INT    converted the input data matrix by log function (log2 or log10,default: ignore this step)		 
         -k              INT    number of subclusters(defaul: 10)  
         -w              FLOAT  plot width(defalt will be 12)                                                                      
         -h              FLOAT  plot height(defalt will be:16)
	 -cls            FLOAT  colum text size,default 1
         -rls            FLOAT  row text size,default 1
         -ml             FLOAT  low margin,default 8
         -mr             FLOAT  right margin,default 8
         -clt	         STRING  cluster type:both,row,colum or none (default,both)
USAGE
	 
die $usage if (!$opts{i});

#define defalts
$opts{i}=abs_path($opts{i});
$opts{o}=$opts{o}?rel2abs($opts{o}):"heatmap";
$opts{log}=$opts{log}?$opts{log}:0;
$opts{k}=$opts{k}?$opts{k}:10;
$opts{w}=$opts{w}?$opts{w}:12;
$opts{h}=$opts{h}?$opts{h}:16;
$opts{cls}=$opts{cls}?$opts{cls}:1;
$opts{rls}=$opts{rls}?$opts{rls}:0;
$opts{ml}=$opts{ml}?$opts{ml}:12;
$opts{mr}=$opts{mr}?$opts{mr}:12;
$opts{clt}=$opts{clt}?$opts{clt}:"both";

my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

my $line = `wc -l $opts{i}`;
$line =~ s/\s.*//g;
if($line<20){
	$opts{h} = 8;
}
if($line<10){
	$opts{h} = 7;
}
if($line<5){
	$opts{h} = 6 ;
}




   
open RCMD, ">$opts{o}.r";
print RCMD "
options(warn=-100)
input_matrix<-\"$opts{i}\"
logNorm<-$opts{log}
subclustNum<-$opts{k}
width<-$opts{w}
height<-$opts{h}
clsize<-$opts{cls}
rlsize<-$opts{rls}
mlow<-$opts{ml}
mright<-$opts{mr}
cltype<-\"$opts{clt}\" #### both row column none

library(cluster)
library(gplots)
library(Biobase)

data = read.delim(input_matrix, header=T, check.names=F, sep=\"\t\")
rownames(data) = data[,1] # set rownames to gene identifiers
data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
myheatcol = redgreen(75)[75:1]

if(logNorm!=0){
data = log(data+1,base=logNorm)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric=\"euclidean\") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method=\"spearman\")), method=\"complete\") # cluster conditions
final_data<-centered_data
}
if(logNorm==0){
hc_genes = agnes(data,diss=FALSE, metric=\"euclidean\") # cluster genes
hc_samples = hclust(as.dist(1-cor(data, method=\"spearman\")), method=\"complete\") # cluster conditions
final_data<-data
}
if(cltype==\"both\"){Rowv=as.dendrogram(hc_genes);Colv=as.dendrogram(hc_samples)}
if(cltype==\"row\"){Rowv=as.dendrogram(hc_genes);Colv=NA}
if(cltype==\"column\"){Rowv=NV;Colv=as.dendrogram(hc_samples)}
if(cltype==\"none\"){Rowv=NA;Colv=NA}

#gene_partition_assignments <- cutree(as.hclust(hc_genes), k=subclustNum);
#partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
#gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file=\"all.RData\")

### cexRow cexCol

### heatmap-plot
#heatmap_filename<-\"Heatmap.pdf\"
heatmap_filename<-paste(\"$opts{o}\",\".pdf\",sep=\"\")
pdf(file=heatmap_filename, width=width,height=height, paper=\"special\")
heatmap.2(final_data, dendrogram=cltype,Rowv=Rowv,Colv=Colv,col=myheatcol, scale=\"none\", density.info=\"none\", trace=\"none\",cexCol=1, cexRow=1,lhei=c(1,3), lwid=c(2,4),margins=c(mlow,mright))
dev.off()

## output the odered matrix after clustered
	if(cltype==\"both\"){
		order_mat<-data[hc_genes\$order[nrow(data):1],hc_samples\$order]
		write.table(order_mat,paste(\"$opts{o}\",\".pdf.orderedmat\",sep=\"\"),sep=\"\t\",col.names=T,row.names=T,quote=F)
	}
	if(cltype==\"row\"){
		order_mat<-data[hc_genes\$order[nrow(data):1],]
		write.table(order_mat,paste(\"$opts{o}\",\".pdf.orderedmat\",sep=\"\"),sep=\"\t\",col.names=T,row.names=T,quote=F)
	}
";

#system ("R --restore --no-save < $opts{o}.r");
`R --restore --no-save < $opts{o}.r`;
#`R --restore --quiet --no-save <cmd.r`
system('rm *.r');
if(-e "all.RData"){system("rm all.RData");}



