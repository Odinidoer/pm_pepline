#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"t=s","o=s","col=s","mb=f","ml=f");
my $usage = <<"USAGE";
Usage:		perl /home/wangyan/script/kegg_pathway_top20bars.pl -t pathways/pathway_table.xls -o pathway
Options:
	-t*	STRING	KEGG annotation file(pathway_table.xls)
	-o*	STRING	output file prefix, eg: pathway
	-col	STRING	barplot color(default: #4F81BD)
	-mb	FLOAT	margin for bottom white space(default:15)
	-ml	FLOAT	margin for left white space(default:10)
USAGE

die $usage if ( !defined $opts{t});
#define defalts
open (my $fh, $opts{t}) or die $!;
open (my $fo, ">bar.temp") or die $!;

while(<$fh>){
	chomp;
	my @line = split(/\t/,$_);
	if($line[0] =~ /01\d\d\d/){
		next;
	}else{
		print $fo $_."\n";
	}
}


close $fh;
close $fo;

$opts{o}=defined$opts{o}?$opts{o}:"pathway";
$opts{col}=defined$opts{col}?$opts{col}:"#4F81BD";
$opts{mb}=defined$opts{mb}?$opts{mb}:15;
$opts{ml}=defined$opts{ml}?$opts{ml}:10;




open RCMD, ">$opts{t}.r";
print RCMD "
options(warn=-1)
pathway_annot<-read.delim(\"bar.temp\",sep=\"\t\",header=T,check.names=F)
pathway_annot1<-pathway_annot[,1:3]
pathway_annot_sort_top<-pathway_annot1[order(pathway_annot1[[3]],decreasing=T),][1:20,]
numbers<-as.numeric(pathway_annot_sort_top[[3]])
labels<-as.character(pathway_annot_sort_top[[2]])
jump<-max(numbers)/10

pdfname <-\"$opts{o}.top20.pdf\"
pdf(file=pdfname,width=12,height=8)
par(mar = c($opts{mb},$opts{ml},3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col=\"$opts{col}\",border=F,horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(\" Number of Proteins \",sep=\"\"),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=45,adj=1,xpd=T,cex=1,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=1,font=1)
dev.off()

";

system ("R --restore --no-save < $opts{t}.r");
#`R --restore --no-save < $opts{t}.r`;
`rm -rf $opts{t}.r`;



