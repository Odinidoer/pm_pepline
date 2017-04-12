#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"spe=s","o=s","t=s","tp=s","scal=f","dp=s","g=s","gl=s","lp=s","cex=f","pcex=f","w=f","h=f","ps=s","td=s","aa=s","ag=f");
my $usage = <<"USAGE";
	Program : $0
	Discription: plot pca & ca from exp table or taxlevel file.
	Usage:perl $0 [options]
		-spe*	input spe tables
			########################################
			#     groups   species1 species2 species3...
			#     sampleA      1        3       3
			#     sampleB      3        4       5
			#     sampleC      3        2       1
			#        .
			#########################################
		-o	output file predix
		-t	T/F    default:T
			if input table is exp table or taxlevel table, then T;                                                                                        
		-tp	pca   plot pda from spe table   default:pca                                                                 
			ca    plot ca from spe table 
		-scal    
			1  relationship between species
			2  relationship between samples  default:2
			3  relationship between samples and species
		-dp       display the sites or the species in the plot, 1;2;12
			1 sp      plot species
			2 si      plot sites     default:2
		-g    groups file for samples eg. [  a1 A
                                             b1 B
                                             b2 B ...]
		-gl  [l/t/b]  group legend or point text or both . defalt:t             
		-w       420
		-h       450
		-cex     0.8
		-pcex    0.9
		-td  [T/F] defalt:F plot 3d pca   if T,dp=2
		-ag  3d plot angle defalt:30 du
		-aa  00-FF (Where the color is "#RRGGBBAA" and the AA portion is the trasparency)defalt :FF
	Example:$0 -spe exp.xls (2d pca)
		$0 -spe exp.xls -td T (3d pca)
                                   
USAGE
die $usage if (!($opts{spe} ));
my $pdf;
my $spe;
my $gp =$opts{g}?1:0;
$opts{t}=$opts{t}?$opts{t}:"T";
$opts{scal}=$opts{scal}?$opts{scal}:2;
$opts{tp}=$opts{tp}?$opts{tp}:"pca";
$opts{dp}=$opts{dp}?$opts{dp}:"2";
$opts{pcex}=$opts{pcex}?$opts{pcex}:1.0;
$opts{gl}=$opts{gl}?$opts{gl}:"t";
$opts{lp}=$opts{lp}?$opts{lp}:"topright";
$opts{ag}=$opts{ag}?$opts{ag}:30;
$opts{aa}=$opts{aa}?$opts{aa}:"FF";
$opts{td}=$opts{td}?$opts{td}:"F";
$opts{cex}=$opts{cex}?$opts{cex}:0.8;
$opts{w}=$opts{w}?$opts{w}:420;
$opts{h}=$opts{h}?$opts{h}:450;

if($opts{td}=~/T/){$pdf ="3d.".$pdf;}

open RCMD, ">cmd.r";
print RCMD "
mycol <-c(\"#000000\",\"#BE0027\",\"#3A89CC\",\"#769C30\",\"#D99536\",\"#7B0078\",\"#BFBC3B\",\"#4C8B35\",\"#3C68AE\",\"#C10077\",\"#CAAA76\",\"#2E165B\",\"#458B00\",\"#8B4513\",\"#008B8B\",\"#6E8B3D\",\"#8B7D6B\",\"#7FFF00\",\"#FF1493\",\"#FF69B4\")
mypch <-c(21:25,3,4,7,9,8,10,15:18,0:14)
col <- 1
gl <-\"$opts{gl}\"
gp <-$gp
if(gp==1){
        group <- read.table(\"$opts{g}\")
        co <-1
        gcol <- as.matrix(unique(group[,2]))
        pcol <- sapply(1:nrow(group),function(x) which(gcol[,] \%in% group[x,2])+1)
        ppch <- pcol-1 
        pcol <- mycol[pcol]     
        ppch <- mypch[ppch]               
}else{
        gcol <-\"#0000FF\"   #blue
        pcol <-\"#0000FF\"  #blue
        ppch <-20
}

specol <- \"#008B8B\"
spepch <-25

td <-\"$opts{td}\";

	spe<-read.table(file=\"$opts{spe}\",header=T,check.names=FALSE,sep=\"\t\")
	name0_spe <- spe[,1]
	rownames(spe) <-spe[,1]
	rownames(spe) <-sapply(rownames(spe),function(x) gsub(\"_*{\.+}\",\" \",x,perl = TRUE))
	write.table(rownames(spe),\"rowname.xls\",sep=\"\t\") 
	spe <-spe[,-1]
	spe_head <-colnames(spe)
	#write.table(spe,\"spe.xls\",sep=\"\t\")

	library(vegan)
	method<-\"$opts{tp}\";
	displ<-\"$opts{dp}\";
	tur<-\"$opts{t}\";

############# PCA ###################

if(method ==\"pca\")
	{
	    if(tur ==\"T\")
		{
	    	spe.pca <-rda(t(spe),scaling=\"$opts{scal}\")
	    }
		else 
		{
    		spe.pca <-rda(spe,scaling=\"$opts{scal}\")
  		  }

    	pc.sites<-summary(spe.pca)\$sites
    	pc.spe<-summary(spe.pca)\$species
   		pc.cont<-summary(spe.pca)\$cont\$importance
		
		basename <- \"$opts{o}\"
		sites <- paste(basename,\".PCA.sites.xls\",sep=\"\")
		spe <- paste(basename,\".PCA.spe.xls\",sep=\"\")
		cont <- paste(basename,\".PCA.cont.xls\",sep=\"\")
		pdf <- paste(basename,\".PCA.pdf\",sep=\"\")
    	write.table(pc.sites,sites,sep=\"\t\",row.names=TRUE)
    	write.table(pc.spe,spe,sep=\"\t\",row.names=TRUE)
    	write.table(pc.cont,cont,sep=\"\t\",row.names=TRUE,col.names=NA)
   
    	pc <- pc.cont[2,]*100

    	axis<-read.table(file=sites,header=T,check.names=FALSE,sep=\"\t\") 
	#write.table(axis,\"axis.xls\")      
        pca <-axis           

   ###### display samples 2D-pca ####### 
    
    
     if(displ==2){

        if(td ==\"F\"){  
		#tiff(\"$opts{spe}.pca.si.tiff\",w=$opts{w},h=$opts{h})
        pdf(pdf)		
        
	pca <-pca[,1:2]
	me<-0.1*abs(max(pca)-min(pca))  

	par(mar=c(5,5,5,2))
	plot(pca,pch=ppch,col=pcol,bg=paste(pcol,\"$opts{aa}\",sep=\"\") ,xlim=c(min(pca[,1])-me,max(pca[,1])+me*3),xlab=paste(\"PC1  \",round(pc[1],2),\"%\",sep=\"\"),ylab=paste(\"PC2  \",round(pc[2],2),\"%\",sep=\"\"),main=\"PCA\",cex=$opts{pcex},las=1)
        #abline(h=0)
        #abline(v=0)

        library(\"maptools\")
          if(gl==\"t\"||gl==\"b\"){
                  pointLabel(x=pca[,1],y=pca[,2],labels=rownames(pca),cex=$opts{cex},col=pcol)
          }
          if(gl==\"l\"||gl==\"b\")
		{
                  legend(\"$opts{lp}\",legend=gcol,col=mycol[2:(length(gcol)+1)],pch=mypch[1:length(gcol)],pt.bg=paste(mycol[2:(length(gcol)+1)],\"$opts{aa}\",sep=\"\"))
	}

######3D-PCA#######
        
        }
	else if(td==\"T\")
	{
	  
	  #tiff(\"$opts{spe}.pca.si.3d.tiff\",w=$opts{w},h=$opts{h},pointsize=12)
	  pdf(\"$opts{spe}.pca.si.3d.pdf\")
  
          library(scatterplot3d)
	  layout(matrix(c(1,2),1,2),width=c(3,1))
          pca <-pca[,1:3]
	  par(mar=c(3,3,0,2))
          s3d <-scatterplot3d(pca,color=\"white\",pch=1,main=\"3D-PCA\",type=\"p\",angle=$opts{ag},scale.y=1.1,box=TRUE,xlab=paste(\"PC1  \",round(pc[1],2),\"%\",sep=\"\"),ylab=paste(\"PC2  \",round(pc[2],2),\"%\",sep=\"\"),zlab=paste(\"PC3  \",round(pc[3],2),\"%\",sep=\"\"))
          s3d\$points3d(pca,col=pcol,pch=ppch,type=\"p\",bg=paste(pcol,\"$opts{aa}\",sep=\"\"))
	  plot.new() 
          par(mar=c(0,1,3,0))         
          legend(\"$opts{lp}\",legend=gcol,col=mycol[2:(length(gcol)+1)],pch=mypch[1:length(gcol)],bty=\"n\",pt.bg=paste(mycol[2:(length(gcol)+1)],\"$opts{aa}\",sep=\"\"))          
   } 
  
   ###### display sample and species #######

	}
	else if(displ==12)
	{
       		#tiff(\"$opts{spe}.pca.si-sp.tiff\",w=$opts{w},h=$opts{h})
			pdf(\"pca.si-sp.pdf\")
      
	plot(pca,pch=ppch,col=pcol,bg=paste(pcol,\"$opts{aa}\",sep=\"\") ,xlim=c(min(pca[,1])-me,max(pca[,1])+me*3),xlab=paste(\"PC1  \",round(pc[1],2),\"%\",sep=\"\"),ylab=paste(\"PC2  \",round(pc[2],2),\"%\",sep=\"\"),main=\"PCA\",cex=$opts{pcex},las=1)
        
        #abline(h=0)
        #abline(v=0)

        points(pc.spe[,1:2],pch=spepch,col=specol,bg=paste(specol),cex=$opts{pcex})
        text(pc.spe[,1:2]*1.05,rownames(spe),col=specol,cex=$opts{cex})

        library(\"maptools\")
        if(gl==\"t\"||gl==\"b\"){
                  pointLabel(x=pca[,1],y=pca[,2],labels=rownames(pca),cex=$opts{cex},col=pcol)
          }
        if(gl==\"l\"||gl==\"b\"){
                  legend(\"$opts{lp}\",legend=gcol,col=mycol[2:(length(gcol)+1)],pch=mypch[1:length(gcol)],pt.bg=paste(mycol[2:(length(gcol)+1)],\"$opts{aa}\",sep=\"\"))
	}

	}



}else if(method ==\"ca\"){
#############  CA  ###################  
	    if(tur ==\"T\"){
	    	spe.ca <-cca(t(spe),scaling=\"$opts{scal}\")
	    }else {
    		spe.ca <-cca(spe,scaling=\"$opts{scal}\")
  		  }

    	pc.sites<-summary(spe.ca)\$sites
    	pc.spe<-summary(spe.ca)\$species
   	pc.cont<-summary(spe.ca)\$cont\$importance
    	write.table(pc.sites,\"pc.sites.xls\",sep=\"\t\",row.names=TRUE)
    	write.table(pc.spe,\"pc.spe.xls\",sep=\"\t\",row.names=TRUE)
    	write.table(pc.cont,\"pc.cont.xls\",sep=\"\t\",row.names=TRUE,col.names=NA)
     
    	pc <- pc.cont[2,]*100

    	axis<-read.table(file=\"pc.sites.xls\",header=T,check.names=FALSE,sep=\"\t\")       
        ca <-axis
                  
	ca <-ca[,1:2]
	me<-0.1*abs(max(ca)-min(ca))  

   ###### display samples ####### 
   	
	if(displ==2){
		#tiff(\"$opts{spe}.ca.si.tiff\",w=$opts{w},h=$opts{h})
        pdf(\"ca.si.pdf\")
  
	plot(ca,pch=ppch,col=pcol,bg=paste(pcol,\"$opts{aa}\",sep=\"\") ,xlim=c(min(ca[,1])-me,max(ca[,1])+me*3),xlab=paste(\"CA1  \",round(pc[1],2),\"%\",sep=\"\"),ylab=paste(\"CA2  \",round(pc[2],2),\"%\",sep=\"\"),main=\"CA\",cex=$opts{pcex},las=1)
        #abline(h=0)
        #abline(v=0)
        library(\"maptools\")
        if(gl==\"t\"||gl==\"b\"){
                  pointLabel(x=ca[,1],y=ca[,2],labels=rownames(ca),cex=$opts{cex},col=pcol)
          }
          if(gl==\"l\"||gl==\"b\")
		{
                  legend(\"$opts{lp}\",legend=gcol,col=mycol[2:(length(gcol)+1)],pch=mypch[1:length(gcol)],pt.bg=paste(mycol[2:(length(gcol)+1)],\"$opts{aa}\",sep=\"\"))
		}

  
   ###### display sample and species #######

	}
	else if(displ==12)
	{
       		#tiff(\"$opts{spe}.ca.si-sp.tiff\",w=$opts{w},h=$opts{h})  
            pdf(\"ca.si-sp.pdf\") 
	plot(ca,pch=ppch,col=pcol,bg=paste(pcol,\"$opts{aa}\",sep=\"\") ,xlim=c(min(ca[,1])-me,max(ca[,1])+me*3),xlab=paste(\"CA1  \",round(pc[1],2),\"%\",sep=\"\"),ylab=paste(\"CA2  \",round(pc[2],2),\"%\",sep=\"\"),main=\"CA\",cex=$opts{pcex},las=1)
        
        #abline(h=0)
        #abline(v=0)

        points(pc.spe[,1:2],pch=spepch,col=specol,bg=paste(specol),cex=$opts{pcex})
        text(pc.spe[,1:2]*1.05,rownames(spe),col=specol,cex=$opts{cex})

        library(\"maptools\")
        if(gl==\"t\"||gl==\"b\"){
                  pointLabel(x=ca[,1],y=ca[,2],labels=rownames(ca),cex=$opts{cex},col=pcol)
          }
          if(gl==\"l\"||gl==\"b\")
		{
                  legend(\"$opts{lp}\",legend=gcol,col=mycol[2:(length(gcol)+1)],pch=mypch[1:length(gcol)],pt.bg=paste(mycol[2:(length(gcol)+1)],\"$opts{aa}\",sep=\"\"))
		}

	}



}



dev.off()
";
system ('R --restore --no-save < cmd.r');

