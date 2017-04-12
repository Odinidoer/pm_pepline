#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s");
my $usage = <<"USAGE";
	Usage: perl $0 [options]
	Description: This program is used for ploting GO enrichment enrichment analysis results.
	Enrichment file name format: *enrichment.detail.xls
	Contact: ting.kuang\@majorbio.com
	Version: v1.2016.04.26
	Options:
		-i	input file, GO enrichment result table[required]
			(must be obtained from exact_goatools.pl)
USAGE

die $usage if ( !defined $opts{i});
open GO, $opts{i} or die $!;
open GOO, ">$opts{i}.go" or die $!;
while(<GO>){
	chomp;
	my @tmp = split /\t/;
	if(/^id/){
		print GOO "description\tratio\tp_bonferroni\ttype\tpvalue\n";
	}else{
		my $fenzi = (split /\//, $tmp[3])[0];
		my $fenmu = (split /\//, $tmp[4])[0];
		my $ratio = $fenzi / $fenmu;
		next if($tmp[1] =~  /p/);
		print GOO "$tmp[2]\t$ratio\t$tmp[6]\t$tmp[7]\t$tmp[5]\n";
	}
}
close GO;
close GOO;
system("head -81 $opts{i}.go > $opts{i}.go.draw");

open RCMD, ">$opts{i}.r";
print RCMD "
GO_enrichraw <- read.delim(\"$opts{i}.go.draw\",header=T,sep=\"\t\",check.names=F)
GO_enrich <- GO_enrichraw[order(GO_enrichraw[,4],GO_enrichraw[,3],GO_enrichraw[,5],-GO_enrichraw[,2],decreasing=F),]

#define p<0.001 color: 
GO_enrich.p1 <- subset(GO_enrich, (pvalue < 0.001))
GO_enrich.p1 <- cbind(GO_enrich.p1, rep(\"***\", nrow(GO_enrich.p1)))
colnames(GO_enrich.p1)[6] <- \"color\"

##define 0.001<=p<0.01 color: 
GO_enrich.p2 <- subset(GO_enrich, (pvalue >= 0.001 & pvalue < 0.01))
GO_enrich.p2 <- cbind(GO_enrich.p2, rep(\"**\", nrow(GO_enrich.p2)))
colnames(GO_enrich.p2)[6] <- \"color\"

##define 0.01<=p<0.05 color: 
GO_enrich.p3 <- subset(GO_enrich, (pvalue >= 0.01 & pvalue < 0.05))
GO_enrich.p3 <- cbind(GO_enrich.p3, rep(\"*\", nrow(GO_enrich.p3)))
colnames(GO_enrich.p3)[6] <- \"color\"

##define p>=0.05 color: 
GO_enrich.p4 <- subset(GO_enrich, (pvalue >= 0.05))
GO_enrich.p4 <- cbind(GO_enrich.p4, rep(4, nrow(GO_enrich.p4)))
colnames(GO_enrich.p3)[6] <- \"color\"

GO_enrich.t <- rbind(GO_enrich.p1, GO_enrich.p2, GO_enrich.p3, GO_enrich.p4)
GO_enrich.t\$color <- as.factor(GO_enrich.t\$color)

GO_enrich.t <- GO_enrich.t[order(GO_enrich.t[,4],GO_enrich.t[,3],GO_enrich.t[,5],decreasing=F),]

### pathway class:
name2typeI<-GO_enrich.t[,c(1,4)]
shortName<-gsub(\"biological_process\",\"BP\",gsub(\"cellular_component\",\"CC\",gsub(\"molecular_function\",\"MF\",name2typeI[[2]])))

xlabs<-paste(GO_enrich.t[[1]],shortName,sep=\": \")
typeII<-c(\"GO Class:\",\"BP: Biological Process\",\"CC: Cellular Component\",\"MF : Molecular Function\")

pdf(\"$opts{i}.go.pdf\", w=16, h=10)

library(\"ggplot2\")
#windowsFonts(myFont=windowsFont(\"Arial\")) 
p <- ggplot(data=GO_enrich.t) + geom_bar(aes(x=factor(xlabs,level=xlabs), y=ratio, fill = color), stat=\"identity\", position=position_dodge(0.4),width=0.8) + 
  geom_text(aes(x = xlabs, y = ratio, label = color), colour = \"black\", size = 3) +
  ggtitle(\"\") + xlab(\"\") + ylab(\"EnrichmentRatio: Sample_number/Background_number\") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x=element_text(angle=45, colour=\"black\", hjust=1)) + 
  theme(axis.title.y= element_text(hjust=0.5, vjust=0.5)) + scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0.0, 1.0, 0.2)) + 
  scale_fill_manual(values = c(\"#228b22\", \"#3cb371\", \"#00cd66\", \"#00fa9a\"), name=\"Pvalue\", breaks=c(\"***\", \"**\", \"*\", \"\"), labels=c(\"p<0.001\",\"0.001<p<0.01\",\"0.01<p<0.05\",\"p>0.05\"))  

library(\"grid\")
subvp <- viewport(x = 0.58, y = 0.45, width = 0.85, height = 0.85)
print(p, vp = subvp)

dev.off()

";

`R --restore --no-save < $opts{i}.r`;

`rm $opts{i}.r`;

