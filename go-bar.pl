#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","a=i","w=i","h=i","off=f");
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: plot bar chart
       Usage :perl $0 [options]
                   -i  input go file
                   -a  all proteins number  (you can provide the number in first line of input file)
                   -w  pdf width  defalt:12
                   -h  pdf heigth  defalt:10
                   -off offset defalt:1.005
                 
       example:
USAGE
die $usage if ( !$opts{i});

#define defalts
$opts{w}=$opts{w}?$opts{w}:12;
$opts{h}=$opts{h}?$opts{h}:10;
$opts{off}=$opts{off}?$opts{off}:1.005;
my $skip =0;
if(!$opts{a}){
      open IN,"<$opts{i}";
      my $l1=<IN>;
      if($l1=~/\s+(\d+)\s/){
           $opts{a}=$1;$skip =1;
      }else{
           print "Error:no proteins number provided !\n";
           exit;
      }      
}
print "Total gene number is $opts{a}\n.\n\n";

open RCMD, ">cmd.r";
print RCMD " def.par <- par(no.readonly = TRUE)
mycol <-c(\"#6B8E23\",\"#7EC0EE\",\"#FF69B4\")


dat <-read.table(file=\"$opts{i}\",sep=\"\t\",head=T,skip=$skip)
x <- as.matrix(dat[,c(1,3)])
rownames(x) <-dat[,2]
x <- x[order(x[,1]),]

g <-table(x[,1])
bcol <-c(rep(1,g[1]),rep(2,g[2]),rep(3,g[3]))

pdf(\"$opts{i}.pdf\",width=$opts{w},height=$opts{h})
par(mar=c(5,2,1,2)+0.1)

x0 <- as.numeric(x[,2])*100/$opts{a}
x0 <- sapply(x0,function(x) if(x<0.1)x=0.1 else x=x)

h <-$opts{h}
ylim0=-5*$opts{h}
ylim1=4*$opts{h}

lpx <-log(x0,base=10)+1
lpx <-h*lpx


xa <- barplot(lpx,beside=TRUE,space=1.2,col=mycol[bcol],axes=FALSE,asp=1,ylim=c(ylim0,ylim1),border=NA,plot=FALSE)
xlim1 <-max(xa)*1.02 
xlim0 <--xlim1*0.15
xa <- barplot(lpx,beside=TRUE,space=1.2,col=mycol[bcol],axes=FALSE,asp=1,ylim=c(ylim0,ylim1),xlim=c(xlim0,xlim1*1.1),border=NA)
text(xa,y=-1,srt=70,adj=1,labels=rownames(x),xpd=T,col=mycol[bcol])

axis(side=2, at=c(0,h,2*h,3*h),labels = c(0,1,10,100),las=1,cex=0.8,pos=c(0,0))
axis(side=4, at=c(0,h,2*h,3*h),labels = c(0,as.integer($opts{a}/100),as.integer($opts{a}/10),as.integer($opts{a})),las=1,cex=0.8,pos=c(xlim1,0))

axis(side=2, at=h,labels =\"Percent of Proteins\" ,las=3,cex.axis=0.9,pos=c(0,0),mgp=c(10,4,0),font.axis=3)
axis(side=4, at=h,labels =\"Number of Proteins\",las=3,cex.axis=0.9,pos=c(xlim1,0),mgp=c(10,4,0),font.axis=3)
#text(x=c(-xlim1*0.05,xlim1*1.05),y=c(15,15),labels=c(\"Percent of Proteins\",\"Number of Proteins\"), cex=0.9, font=3,srt=90)
#text(x=c(xlim0*0.8,xlim1*1.1),y=c(15,15),labels=c(\"Percent of Protein\",\"Number of Protein\"), cex=0.9, font=3,srt=90)

labw <-strwidth(rownames(x),units=\"user\")
ylow <- -max(labw) 
#abline(h=c(0,10,20,30,ylow ))
pos <- ylow*tan(pi*21.5/180)
tx0 <- c(rep(0,4),pos)
tx1 <- c(rep(xlim1,4),xlim1+pos)
ty0 <- ty1 <- c(0,h,2*h,3*h,ylow)
segments(tx0,ty0,tx1,ty1)

lg <-c(0,g)
lx0 <-c(0)
for(i in 2:length(g)){
      lg[i] =lg[i-1]+lg[i]
      lx0[i]=0.5*(xa[lg[i]]+xa[lg[i]+1])*$opts{off}
}
lx0 <- c(lx0,xlim1)
ly0 <- rep(0,length(lx0))
lx1 <- lx0+pos
ly1 <- rep(ylow,length(lx0))
segments(lx0,ly0,lx1,ly1)

ttx <-c()
for(i in 1:3){
      ttx[i] =0.5*(lx1[i]+lx1[i+1])
}

text(x=ttx,y=rep(ylow*1.01,3),labels=names(g), cex=0.9,col=mycol[c(1,2,3)],pos=1)


par(def.par)
dev.off()
";

system ('R --restore --no-save < cmd.r');
system ('rm cmd.r');

print "output: $opts{i}.pdf\n";









