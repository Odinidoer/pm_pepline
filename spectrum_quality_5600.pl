#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","d=s","b=s","c=s");
my $usage = <<"USAGE";
		Program : $0
		Version : 1.0
		Discription: Mass spectrometry data quality evaluation
		Usage :perl $0 [options]
			-i	PetitdeSummary.txt
			-d	out of dmass_score.plot prefix
			-b	base.txt  #base informations of mass spectra
			-c	ProteinSummary.txt
		example: $0 -i PetitdeSummary.txt -d dmass_score -b base.txt -c ProteinSummary.txt
USAGE
die $usage if ( !$opts{i});

#plot dmass_score 
system("awk -F \"\t\" \'{print \$16 \"\t\"\$22}\' $opts{i}> $opts{d}.txt ");

open RCMD, ">$opts{d}.r";
print RCMD " options(warn=-100)
library(\"ggplot2\")
df<-read.table(\"$opts{d}.txt\",sep=\"\t\",header=T)
p<-ggplot(data=df,aes(x=df[,1],y=df[,2]))
p<-p+geom_jitter(color=\"red\",size=1)+ scale_x_continuous(limits = c(-4,4))+theme_bw() 
p<-p+theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())
p<-p+xlab(\"Mass delta(ppm)\")+ylab(\"Ion Score\")
p<-p+theme(axis.title= element_text(size = 15))
ggsave('$opts{d}.pdf', p)
";

#plot base informations 
open RCMD, ">$opts{b}.r";
print RCMD " options(warn=-100)
library(\"ggplot2\")
df<-read.table(\"$opts{b}\",sep=\"\t\",header=T)
df<-df[order(-df[,'count']),]
p<-ggplot(data=df)
p<-p+geom_bar(aes(x=factor(attribute,level=attribute), y=count), stat='identity',fill='\#436EEE')
p<-p+labs(title=\"base information\")+theme(plot.title=element_text(size=15), face=\"bold\")

p<-p+theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())
p<-p+xlab(\"\")+ylab(\"number\")
p<-p+theme(axis.text.x=element_text(angle=45,colour=\"black\",hjust=1))

jump<-max(df\$count)/10
p<-p+geom_text(aes(x=df\$attribute,y=df\$count+jump/3),label=(df\$count))
p<-p+labs(title=\"Base information statistics\")
p<-p+theme(plot.title=element_text(size=15, face=\"bold\"))


ggsave('$opts{b}.pdf', p)
";

#plot pie of the cover of peptides
system("awk -F \"\t\" '{print \$10 }' $opts{c}\| sed '1d'> peptide.txt ");

open RCMD, ">$opts{c}.r";
print RCMD " options(warn=-100)
df<-read.table(\"peptide.txt\",sep=\"\t\")
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
pdf(\"$opts{c}.pie.pdf\",width=8, height=8)
pie(number_plot\$number_protein,labels = number_plot\$number_ratio,col=brewer.pal(length(number_plot\$number_protein), \"Set3\"), border=F, radius=0.8,cex=0.3,main=\"Pie number distribution\")
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
p<-p+labs(x=\"Peptide number\",y=\"Protein number\",title=\"Peptide number distrution\")
p<-p+theme(plot.title=element_text(size=15, face=\"bold\"))
ggsave(\"$opts{c}.bar.pdf\",p)

";

`R --restore --no-save < $opts{d}.r`;

`R --restore --no-save < $opts{b}.r`;

`R --restore --no-save < $opts{c}.r`; 

system ('rm *.r');
system('rm peptide.txt');
system("rm $opts{d}.txt");
 
