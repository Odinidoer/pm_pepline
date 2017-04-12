#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","o=s");
my $usage = <<"USAGE";
		Program : $0
		Version : 1.0
		Discription: Mass spectrometry data quality evaluation
		Usage :perl $0 [options]
			-i	ProteinSummary.txt
			-o	outfile prefix
		example: $0 -i ProteinSummary.txt -o Peptides_distrubution
USAGE

die $usage if ( !$opts{i});
$opts{o}= defined $opts{o}?$opts{o}:"Peptides_distrubution";

system("awk -F \"\t\" '{print \$7 }' $opts{i}\| sed '1d'> $opts{o}.txt ");

open RCMD, ">$opts{o}.r";
print RCMD " options(warn=-100)
df<-read.table(\"$opts{o}.txt\",sep=\"\t\")
df1<-length(which(df>0&df<6))
df2<-length(which(df>=6&df<11))
df3<-length(which(df>=11&df<16))
df4<-length(which(df>=16&df<21))
df5<-length(which(df>=21&df<26))
df6<-length(which(df>=26&df<31))
df7<-length(which(df>=31&df<36))
df8<-length(which(df>=36&df<41))
df9<-length(which(df>=41&df<46))
df10<-length(which(df>=46&df<50))
df11<-length(which(df>=50))
number_protein<-c(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11)
number_class<-as.character(c(\"1-5\",\"6-10\",\"11-15\",\"16-20\",\"21-25\",\"26-30\",\"31-35\",\"36-40\",\"41-45\",\"46-50\",\">50\"))
number_ratio <-round(number_protein/sum(number_protein)*100,2)
number_ratio<- paste(number_ratio,\"%\",sep=\"\")
number_plot<-data.frame(number_class,number_protein,number_ratio)
number_plot<-number_plot[order(number_plot\$number_protein),]

library(\"RColorBrewer\")
pdf(\"$opts{o}.pie.pdf\",width=10, height=8)
pie(number_plot\$number_protein,labels = number_plot\$number_ratio,col=brewer.pal(length(number_plot\$number_protein), \"Set3\"), border=F, radius=0.8,cex=0.3)
pie_legend<-paste(number_plot\$number_class,number_plot\$number_protein,sep=\"\(\")
pie_legend<-paste(pie_legend,\"\",sep=\"\)\")
legend(\"left\",legend=pie_legend,bty=\"n\",fill=brewer.pal(length(number_plot\$number_protein),\"Set3\"),cex=0.8)

dev.off()

#plot bar of the cover of peptides
library(\"ggplot2\")
df_bar<-data.frame(number_class,number_protein)

#df_bar<-df_bar[order(-df_bar['number_protein']),]
p<-ggplot(data=df_bar)

p<-p+geom_bar(aes(x=factor(df_bar\$number_class,level=df_bar\$number_class), y=df_bar\$number_protein), stat='identity',fill='\#436EEE')
p<-p+theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())
p<-p+xlab(\"Peptide number\")+ylab(\"Protein number\")
p<-p+theme(axis.text.x=element_text(angle=45,colour=\"black\",hjust=1))

jump<-max(df_bar\$number_protein)/10
p<-p+geom_text(aes(x=df_bar\$number_class,y=df_bar\$number_protein+jump/3),label=(df_bar\$number_protein))
p<-p+labs(x=\"Peptide number\",y=\"Number of protein in class\",title=\"Distribution of peptide number\")
p<-p+theme(plot.title=element_text(size=15, face=\"bold\"))
ggsave(\"$opts{o}.bar.pdf\",p, w=12, h=8)

";

`R --restore --no-save < $opts{o}.r`; 

system ("rm $opts{o}.r ");
 
