#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","w=i","h=i");
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: plot bar chart
       Usage :perl $0 [options]
                   -i input cog file
                   -w  pdf width defalt:10
                   -h  pdf heigth defalt:6
                 
       example:
USAGE
die $usage if ( !$opts{i});

#define defalts
$opts{w}=$opts{w}?$opts{w}:10;
$opts{h}=$opts{h}?$opts{h}:6;


open RCMD, ">cmd.r";
print RCMD " def.par <- par(no.readonly = TRUE)

cog <-matrix(c(\"A\" , \"RNA processing and modification\",
               \"B\" , \"Chromatin structure and dynamics\",
               \"C\" , \"Energy production and conversion\",
               \"D\" , \"Cell cycle control, cell division, chromosome partitioning\",
               \"E\" , \"Amino acid transport and metabolism\",
               \"F\" , \"Nucleotide transport and metabolism\",
               \"G\" , \"Carbohydrate transport and metabolism\",
               \"H\" , \"Coenzyme transport and metabolism\",
               \"I\" , \"Lipid transport and metabolism\",
               \"J\" , \"Translation, ribosomal structure and biogenesis\",
               \"K\" , \"Transcription\",
               \"L\" , \"Replication, recombination and repair\",
               \"M\" , \"Cell wall/membrane/envelope biogenesis\",
               \"N\" , \"Cell motility\",
               \"O\" , \"Posttranslational modification, protein turnover, chaperones\",
               \"P\" , \"Inorganic ion transport and metabolism\",
               \"Q\" , \"Secondary metabolites biosynthesis, transport and catabolism\",
               \"R\" , \"General function prediction only\",
               \"S\" , \"Function unknown\",
               \"T\" , \"Signal transduction mechanisms\",
               \"U\" , \"Intracellular trafficking, secretion, and vesicular transport\",
               \"V\" , \"Defense mechanisms\",
               \"W\" , \"Extracellular structures\",
               \"Y\" , \"Nuclear structure\",
               \"Z\" , \"Cytoskeleton\"),nrow=25,ncol=2,byrow=TRUE)
cog

dat <- read.table(file=\"$opts{i}\",sep=\"\t\")
dat
x <- c()
cog[1,2]
for (i in 1:25){
     if(length(which(dat[,1] \%in\% cog[i,2]))){
         x[i] <- dat[ which(dat[,1] \%in\% cog[i,2]),2]
     }else{ x[i] <- 0 }

}
x <-as.matrix(x)
row.names(x) <- cog[,1]
x
pdf(\"COG.class.catalog.pdf\",width=$opts{w},height=$opts{h})

#layout(matrix(c(1,2),1,2),width=c(5:6.2))
par(mar=c(3,3,2,0))
xx <-barplot(x,beside=TRUE,col=rainbow(25),space=0.4,las=1,xlim=c(-4,70),ylim=c(-0.3*max(x),max(x)*1.15),cex.names=0.6,main=\"COG Function Classification\",cex.main=0.8,cex.lab=0.8,font.lab=3,font.names=3,axes=FALSE,adj.lab=0)
text(xx,y=-0.05*max(x),labels=rownames(x),xpd=T,cex=0.8)
text(x=-5,y=0.5*max(x),labels=c(\"Number of Proteins\"), cex=0.9, font=3,srt=90)
text(x=20,y=-0.15*max(x),labels=c(\"Function Class\"), cex=0.9, font=3)
rect(0,0,36,max(x)*1.02)
at <-seq(0,max(x),ceiling(max(x)/5))
axis(side=2,las=1,cex=0.8,pos=c(0,0),at=at)
#plot.new()
#par(mar=c(4,0,2,0))
legend(x=38,y=max(x)*1.1,legend=paste(row.names(x),\": \",cog[,2]),cex=0.6,bty=\"n\",y.intersp=1.5)
 
par(def.par)
dev.off()
";

system ('R --restore --no-save < cmd.r');
system ('rm cmd.r Rplots.pdf');
