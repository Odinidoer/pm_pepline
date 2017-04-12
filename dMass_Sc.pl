#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","o=s","h!");
my $usage = <<"USAGE";
		Program : $0
		Version : 1.0
		Discription: Mass spectrometry data quality evaluation
		Usage :perl $0 [options]
			-i	PetitdeSummary.txt
			-o	outfile prefix, defaults: dMass
		example: $0 -i PetitdeSummary.txt -o dMass
USAGE

$opts{o}= defined $opts{o}?$opts{o}:"dMass";

die $usage if ( !$opts{i} || $opts{h});

open RCMD, ">$opts{o}.r";
print RCMD " options(warn=-100)
library(\"ggplot2\")
dat<-read.table(\"$opts{i}\",sep=\"\t\",header=T)
df<-data.frame(dat\$dMass, dat\$Sc)
p<-ggplot(data=df,aes(x=df[,1],y=df[,2]))
p<-p+geom_jitter(color=\"red\",size=1)+ scale_x_continuous(limits = c(-10,10))+theme_bw() 
p<-p+theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())
p<-p+xlab(\"Mass delta(ppm)\")+ylab(\"Ion Score\")
p<-p+theme(axis.title= element_text(size = 15))
ggsave('$opts{o}.pdf', p)
";

`R --restore --no-save < $opts{o}.r`; 

system ("rm $opts{o}.r");
