#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use DBI;
use Getopt::Long;
use List::Util qw ( sum);
my %opts;
my $dbhost  = "localhost";    #
my $dbuser  = "blast2go";     #
my $dbpass  = "blast4it";     #
my $dbname  = "b2g";          #
my $dbport  = "3306";         #
my $dbh;
my @orf;
my %golist;

my $orfname;
my $goid;
GetOptions(
	\%opts,"i=s","o=s","l1=i","l2=i","l3=i","pie=s","bar=s","pw=f","ph=f","pc1=f","pc2=f","bw=f","bh=f","bc1=f","bc2=f","d=s","h=s","p=s","u=s","port=i","help!");

my $usage = <<"USAGE";
       Program : $0
       Usage :perl $0 [options]
                -i* 	[GO files from Interpro-extract]
				-o*      output prefix for the Rscript
				-l1      GO level,2
				-l2      GO level,3
				-l3      GO level,4
				-bar    STRING  [T/F]plot bar or not(default: T ) 
				-pie    STRING  [T/F]plot pie or not(default: T )
                
		pie options :
             	 -pw  	FLOAT    width      default: 20
             	 -ph    FLOAT    height     default: 15
             	 -pc1  FLOAT    text size of pies   default: 1
                 -pc2  FLOAT    text size of legend default: 1
		
		bar options :                        			       
                 -bw    FLOAT    width     default: 12
                 -bh    FLOAT    height    default: 20
				-bc1 	FLOAT    text size for each go name     default: 0.5
				-bc2 	FLOAT    text size for 3 go categaries  default: 1
		 
                -h		[Mysql database host:Defualt "localhost"]
                -u		[Mysql database username:Defualt "blast2go"]
                -p		[Mysql database passwod:Defualt "blast4it"]
                -d		[Mysql GO database name:Defualt "b2g"]
                -port	[Mysql database port:Defualt "3306"]
                -help	Display this usage information
                * 		must be given Argument 
       example:go_9pies.pl -i GO.list -l1 4 -l2 5 -l3 6 -pie T 
USAGE

die $usage
  if ( !( $opts{i} ) || $opts{help} );
$dbhost = $opts{h}    if ( $opts{h} );
$dbuser = $opts{u}    if ( $opts{u} );
$dbpass = $opts{p}    if ( $opts{p} );
$dbname = $opts{d}    if ( $opts{d} );
$dbport = $opts{port} if ( $opts{port} );

#define defaults
$opts{bar}=$opts{bar}?$opts{bar}:"T";
$opts{pie}=$opts{pie}?$opts{pie}:"T";
$opts{pw}=$opts{pw}?$opts{pw}:20;
$opts{ph}=$opts{ph}?$opts{ph}:16;
$opts{pc1}=$opts{pc1}?$opts{pc1}:1;
$opts{pc2}=$opts{pc2}?$opts{pc2}:1;
$opts{bw}=$opts{bw}?$opts{bw}:12;
$opts{bh}=$opts{bh}?$opts{bh}:20;
$opts{bc1}=$opts{bc1}?$opts{bc1}:0.5;
$opts{bc2}=$opts{bc2}?$opts{bc2}:1;
$opts{l1}=$opts{l1}?$opts{l1}:2;
$opts{l2}=$opts{l2}?$opts{l2}:3;
$opts{l3}=$opts{l3}?$opts{l3}:4;

open( FILE, "<$opts{i}" ) || die "Can't open $opts{i} or file not exsit!\n";
while (<FILE>) {
	chomp;
	@orf = split( /\s+/, $_ );
	if ( $#orf < 1 ) { next }
	if ( $orf[0] =~ /^GO:\d+$/ ) {
		die "$opts{i} file format error!\n";
	}
	else {
		for ( my $n = 1 ; $n <= $#orf ; $n++ ) {
			if($orf[$n]=~/\;/){
				my @b=split(/;\s*/,$orf[$n]);
				push( @{ $golist{ $orf[0] }{'GO'} }, @b);
			}
			#until ( $orf[$n] =~ /^GO:\d+$/ ) {
				#die "$opts{i} file format error!\n";
			#}
			else{
				push( @{ $golist{ $orf[0] }{'GO'} }, $orf[$n] );
			}
		}
	}
}
close FILE;

$dbh = DBI->connect( "DBI:mysql:database=$dbname;mysql_socket=/var/lib/mysql/mysql.sock;host=$dbhost;port=$dbport",
	$dbuser, $dbpass )
  || die "Could not connect to database: $DBI::errstr";


my %seqlist;

foreach $orfname ( keys(%golist) ) {
	foreach $goid ( @{ $golist{$orfname}{'GO'} } ) {
		my $sth = $dbh->prepare(
"SELECT DISTINCT ancestor.*, graph_path.term1_id AS ancestor_id FROM term INNER JOIN graph_path ON (term.id=graph_path.term2_id) INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)  WHERE term.acc='$goid';"
		);
		if ( !$sth ) {
			die "Error:" . $dbh->errstr . "\n";
		}
		if ( !$sth->execute ) {
			die "Error:" . $sth->errstr . "\n";
		}
		my $ref = $sth->fetchall_hashref('id');

		#print scalar(keys(%$ref))." \n";
		foreach my $id ( keys(%$ref) ) {

			#print("$ref->{$id}->{'name'} \t $ref->{$id}->{'acc'} \n");
			my $termname = $ref->{$id}->{'name'};
			#unless ( exists( $seqlist{$termname}{$orfname} ) ) {
				#$terms{$termname}++;
				#$golist{$orfname}{$termname} = $ref->{$id}->{'acc'};
				$seqlist{$termname}{$orfname}{$goid}=1;
			#}
		}
		$sth->finish;
	}
}

sub getDescendants() {
	my ($name) = @_;
	my %goterm;
	my $goname;
	my $sth1 = $dbh->prepare(
	"SELECT DISTINCT descendant.acc, descendant.name, descendant.term_type FROM  term  INNER JOIN graph_path ON (term.id=graph_path.term1_id)  INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id) WHERE distance = 1 and term.name=".$dbh->quote($name).";") or die "Error:" . $dbh->errstr . "\n";
	$sth1->execute();
	my $ref = $sth1->fetchall_hashref('acc');
	foreach my $acc ( keys(%$ref) ) {
		$goname = $ref->{$acc}->{'name'};
		$goterm{$goname} = $acc;
	}
	return %goterm;
	$sth1->finish;
}


my $sum =scalar keys(%golist);
open(OUT1,">$opts{i}.level$opts{l1}.txt");
print OUT1 "total proteins number with gos: $sum \n";

my (%bcm);
$bcm{'molecular_function'}{'1'}{'molecular_function'} = "GO:0003674";
$bcm{'biological_process'}{'1'}{'biological_process'} = "GO:0008150";
$bcm{'cellular_component'}{'1'}{'cellular_component'} = "GO:0005575";

### level1
        my $x =$opts{l1};
	#my $x =$opts{l};
	$x--;
	print OUT1 "term_type\tterm\tnumber\tpercent\tGO\n";
	foreach my $key ( keys(%bcm) ) {
		for ( my $i = 1 ; $i<= $x ; $i++ ) {
			foreach my $tname ( keys( %{ $bcm{$key}{$i} } ) ) {
				my %newterm = &getDescendants($tname);
				foreach my $m ( keys(%newterm) ) {
					$bcm{$key}{ $i + 1 }{$m} = $newterm{$m};
				}
			}
		}
		 my %zterms;
		 foreach my $lterm (keys( %{ $bcm{$key}{$opts{l1}} } )){
		   if (exists($seqlist{$lterm})){
		   	$zterms{$lterm}=scalar(keys %{$seqlist{$lterm}});
		   }		   
	     }
	    #if($opts{list}){
	    	open LIST1, ">> $opts{i}.level$opts{l1}.xls" or die "Error:Cannot open file $opts{i}.level$opts{l1} : $! \n";
			print LIST1 "Term_type\tTerm\tGO_id\tNumber_of_protein\tProtein_GO_list\n";
	    	foreach my $nterm (keys(%zterms)){
	    		my @seqs;
	    		foreach my $x (keys %{$seqlist{$nterm}}){
	    			push(@seqs,$x."(".join(",",keys %{$seqlist{$nterm}{$x}}).")");
	    		}
	    		print LIST1 "$key\t$nterm\t".$bcm{$key}{$opts{l1}}{$nterm}."\t".$zterms{$nterm}."\t".join(";",@seqs)."\n";
	    	}
	    	close LIST1;
	    #}
		
	     foreach my $nterm (keys(%zterms)){
	     	my $m = $zterms{$nterm};
	     	my $percent=$m/$sum;
	     	print OUT1 "$key\t$nterm\t$m\t$percent\t$bcm{$key}{$opts{l1}}{$nterm}\n";
	     }
	}

### level2
open(OUT2,">$opts{i}.level$opts{l2}.txt");
print OUT2 "total proteins number with gos: $sum \n";

        my $y =$opts{l2};
	#my $x =$opts{l};
	$y--;
	print OUT2 "term_type\tterm\tnumber\tpercent\tGO\n";
	foreach my $key ( keys(%bcm) ) {
		for ( my $i = 1 ; $i<= $y ; $i++ ) {
			foreach my $tname ( keys( %{ $bcm{$key}{$i} } ) ) {
				my %newterm = &getDescendants($tname);
				foreach my $m ( keys(%newterm) ) {
					$bcm{$key}{ $i + 1 }{$m} = $newterm{$m};
				}
			}
		}
		 my %zterms;
		 foreach my $lterm (keys( %{ $bcm{$key}{$opts{l2}} } )){
		   if (exists($seqlist{$lterm})){
		   	$zterms{$lterm}=scalar(keys %{$seqlist{$lterm}});
		   }		   
	     }
	    #if($opts{list}){
	    	open LIST2, ">> $opts{i}.level$opts{l2}.xls" or die "Error:Cannot open file $opts{i}.level$opts{l2} : $! \n";
			print LIST2 "Term_type\tTerm\tGO_id\tNumber_of_protein\tProtein_GO_list\n";
	    	foreach my $nterm (keys(%zterms)){
	    		my @seqs;
	    		foreach my $y (keys %{$seqlist{$nterm}}){
	    			push(@seqs,$y."(".join(",",keys %{$seqlist{$nterm}{$y}}).")");
	    		}
	    		print LIST2 "$key\t$nterm\t".$bcm{$key}{$opts{l2}}{$nterm}."\t".$zterms{$nterm}."\t".join(";",@seqs)."\n";
	    	}
	    	close LIST2;
	    #}
		
	     foreach my $nterm (keys(%zterms)){
	     	my $m = $zterms{$nterm};
	     	my $percent=$m/$sum;
	     	print OUT2 "$key\t$nterm\t$m\t$percent\t$bcm{$key}{$opts{l2}}{$nterm}\n";
	     }
	}

### level3
open(OUT3,">$opts{i}.level$opts{l3}.txt");
print OUT3 "total proteins number with gos: $sum \n";
        my $z =$opts{l3};
	#my $x =$opts{l};
	$z--;
	print OUT3 "term_type\tterm\tnumber\tpercent\tGO\n";
	foreach my $key ( keys(%bcm) ) {
		for ( my $i = 1 ; $i<= $z ; $i++ ) {
			foreach my $tname ( keys( %{ $bcm{$key}{$i} } ) ) {
				my %newterm = &getDescendants($tname);
				foreach my $m ( keys(%newterm) ) {
					$bcm{$key}{ $i + 1 }{$m} = $newterm{$m};
				}
			}
		}
		 my %zterms;
		 foreach my $lterm (keys( %{ $bcm{$key}{$opts{l3}} } )){
		   if (exists($seqlist{$lterm})){
		   	$zterms{$lterm}=scalar(keys %{$seqlist{$lterm}});
		   }		   
	     }
	    #if($opts{list}){
	    	open LIST3, ">> $opts{i}.level$opts{l3}.xls" or die "Error:Cannot open file $opts{i}.level$opts{l3} : $! \n";
			print LIST3 "Term_type\tTerm\tGO_id\tNumber_of_protein\tProtein_GO_list\n";
	    	foreach my $nterm (keys(%zterms)){
	    		my @seqs;
	    		foreach my $z (keys %{$seqlist{$nterm}}){
	    			push(@seqs,$z."(".join(",",keys %{$seqlist{$nterm}{$z}}).")");
	    		}
	    		print LIST3 "$key\t$nterm\t".$bcm{$key}{$opts{l3}}{$nterm}."\t".$zterms{$nterm}."\t".join(";",@seqs)."\n";
	    	}
	    	close LIST3;
	    #}
		
	     foreach my $nterm (keys(%zterms)){
	     	my $m = $zterms{$nterm};
	     	my $percent=$m/$sum;
	     	print OUT3 "$key\t$nterm\t$m\t$percent\t$bcm{$key}{$opts{l3}}{$nterm}\n";
	     }
	}	
	
$dbh->disconnect();

### plot pies
open RCMD, ">$opts{o}.r";
print RCMD "
options(warn=-1)
################ function of plot pies #################
plot_pies<-function(L1,L2,L3,l1,l2,l3,pdfname,width,height,pcex1,pcex2){
    #x<-c(0,0.2,0.19,0.32,0.31,0.51,0.50,0.63,0.62,0.82,0.81,1)
    x<-c(0,0.2,0.18,0.32,0.30,0.50,0.48,0.62,0.60,0.8,0.78,1)
    y<-c(0.2,0.4,0.41,0.61,0.62,0.82)
    
    #ramp_bp<-colorRamp(as.vector(c(\"turquoise4\",\"white\")))
    ramp_bp<-colorRamp(as.vector(c(\"#6B8E23\",\"white\")))
    heatcol_BP2<-paste(rgb(ramp_bp(seq(0,1,length=10)),max=255),\"E5\",sep=\"\")
    #ramp_cc<-colorRamp(as.vector(c(\"darkorange2\",\"white\")))
    ramp_cc<-colorRamp(as.vector(c(\"#7EC0EE\",\"white\")))
    heatcol_CC2<-paste(rgb(ramp_cc(seq(0,1,length=10)),max=255),\"E5\",sep=\"\")
    #ramp_mf<-colorRamp(as.vector(c(\"olivedrab4\",\"white\")))
    ramp_mf<-colorRamp(as.vector(c(\"#FF69B4\",\"white\")))
    heatcol_MF2<-paste(rgb(ramp_mf(seq(0,1,length=10)),max=255),\"E5\",sep=\"\")    
    ColorList<-brewer.pal(10,\"Paired\")   
  
    pdf(pdfname,width,height)
    layout(matrix(1:9,3,3,byrow=F),heights=rep(1,9))
    
    gostat_level2<-read.delim(l1,header=T,check.names=F,skip=1)
    gostat_level3<-read.delim(l2,header=T,check.names=F,skip=1)
    gostat_level4<-read.delim(l3,header=T,check.names=F,skip=1)
    #### LEVEL-2
    BP_2<-gostat_level2[which(gostat_level2[[1]]==\"biological_process\"),]
    BP_2_sort<-BP_2[order(BP_2[[3]],decreasing=T),][1:ifelse(nrow(BP_2)>10,10,nrow(BP_2)),2:3]   
    BP_2_splice<-as.numeric(as.character(BP_2_sort[[2]]))
    BP_2_label<-as.character(BP_2_sort[[2]])
    #BP_2_legend<-as.character(BP_2_sort[[1]])
    BP_2_legend<-paste(as.character(BP_2_sort[[1]]),\" (\",BP_2_label,\")\",sep=\"\")
    
    CC_2<-gostat_level2[which(gostat_level2[[1]]==\"cellular_component\"),]
    CC_2_sort<-CC_2[order(CC_2[[3]],decreasing=T),][1:ifelse(nrow(CC_2)>10,10,nrow(CC_2)),2:3]
    CC_2_splice<-as.numeric(as.character(CC_2_sort[[2]]))
    CC_2_label<-as.character(CC_2_sort[[2]])
    #CC_2_legend<-as.character(CC_2_sort[[1]])
    CC_2_legend<-paste(as.character(CC_2_sort[[1]]),\" (\",CC_2_label,\")\",sep=\"\")
    
    MF_2<-gostat_level2[which(gostat_level2[[1]]==\"molecular_function\"),]
    MF_2_sort<-MF_2[order(MF_2[[3]],decreasing=T),][1:ifelse(nrow(MF_2)>10,10,nrow(MF_2)),2:3]
    MF_2_splice<-as.numeric(as.character(MF_2_sort[[2]]))
    MF_2_label<-as.character(MF_2_sort[[2]])
    #MF_2_legend<-as.character(MF_2_sort[[1]])
    MF_2_legend<-paste(as.character(MF_2_sort[[1]]),\" (\",MF_2_label,\")\",sep=\"\")
    
    #### LEVEL-3
    BP_3<-gostat_level3[which(gostat_level3[[1]]==\"biological_process\"),]
    BP_3_sort<-BP_3[order(BP_3[[3]],decreasing=T),][1:ifelse(nrow(BP_3)>10,10,nrow(BP_3)),2:3]
    BP_3_splice<-as.numeric(as.character(BP_3_sort[[2]]))
    BP_3_label<-as.character(BP_3_sort[[2]])
    #BP_3_legend<-as.character(BP_3_sort[[1]])
    BP_3_legend<-paste(as.character(BP_3_sort[[1]]),\" (\",BP_3_label,\")\",sep=\"\")
    
    
    CC_3<-gostat_level3[which(gostat_level3[[1]]==\"cellular_component\"),]
    CC_3_sort<-CC_3[order(CC_3[[3]],decreasing=T),][1:ifelse(nrow(CC_3)>10,10,nrow(CC_3)),2:3]
    CC_3_splice<-as.numeric(as.character(CC_3_sort[[2]]))
    CC_3_label<-as.character(CC_3_sort[[2]])
    #CC_3_legend<-as.character(CC_3_sort[[1]])
    CC_3_legend<-paste(as.character(CC_3_sort[[1]]),\" (\",CC_3_label,\")\",sep=\"\")
    
    MF_3<-gostat_level3[which(gostat_level3[[1]]==\"molecular_function\"),]
    MF_3_sort<-MF_3[order(MF_3[[3]],decreasing=T),][1:ifelse(nrow(MF_3)>10,10,nrow(MF_3)),2:3]
    MF_3_splice<-as.numeric(as.character(MF_3_sort[[2]]))
    MF_3_label<-as.character(MF_3_sort[[2]])
    #MF_3_legend<-as.character(MF_3_sort[[1]])
    MF_3_legend<-paste(as.character(MF_3_sort[[1]]),\" (\",MF_3_label,\")\",sep=\"\")
    
    #### LEVEL-4
    BP_4<-gostat_level4[which(gostat_level4[[1]]==\"biological_process\"),]
    BP_4_sort<-BP_4[order(BP_4[[3]],decreasing=T),][1:ifelse(nrow(BP_4)>10,10,nrow(BP_4)),2:3]
    BP_4_splice<-as.numeric(as.character(BP_4_sort[[2]]))
    BP_4_label<-as.character(BP_4_sort[[2]])
    #BP_4_legend<-as.character(BP_4_sort[[1]])
    BP_4_legend<-paste(as.character(BP_4_sort[[1]]),\" (\",BP_4_label,\")\",sep=\"\")
    
    CC_4<-gostat_level4[which(gostat_level4[[1]]==\"cellular_component\"),]
    CC_4_sort<-CC_4[order(CC_4[[3]],decreasing=T),][1:ifelse(nrow(CC_4)>10,10,nrow(CC_4)),2:3]
    CC_4_splice<-as.numeric(as.character(CC_4_sort[[2]]))
    CC_4_label<-as.character(CC_4_sort[[2]])
    #CC_4_legend<-as.character(CC_4_sort[[1]])
    CC_4_legend<-paste(as.character(CC_4_sort[[1]]),\" (\",CC_4_label,\")\",sep=\"\")
    
    MF_4<-gostat_level4[which(gostat_level4[[1]]==\"molecular_function\"),]
    MF_4_sort<-MF_4[order(MF_4[[3]],decreasing=T),][1:ifelse(nrow(MF_4)>10,10,nrow(MF_4)),2:3]
    MF_4_splice<-as.numeric(as.character(MF_4_sort[[2]]))
    MF_4_label<-as.character(MF_4_sort[[2]])
    #MF_4_legend<-as.character(MF_4_sort[[1]])
    MF_4_legend<-paste(as.character(MF_4_sort[[1]]),\" (\",MF_4_label,\")\",sep=\"\")
    
    ## MF-level2
    par(mai=c(0,0,0,0))
    par(fig=c(x[1],x[2],y[1],y[2]),new=F)
	pie(MF_2_splice,labels=MF_2_label,col=heatcol_MF2,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    legend(\"topleft\",\"MF\",bty=\"n\",cex=2,text.col=\"#FF69B4\",text.font=2)
    par(fig=c(x[3],x[4],y[1],y[2]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=MF_2_legend,fill=heatcol_MF2,cex=pcex2,bty=\"n\")
    ## CC-level2
    par(fig=c(x[1],x[2],y[3],y[4]),new=T)
    pie(CC_2_splice,labels=CC_2_label,col=heatcol_CC2,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    legend(\"topleft\",\"CC\",bty=\"n\",cex=2,text.col=\"#7EC0EE\",text.font=2)
    par(fig=c(x[3],x[4],y[3],y[4]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=CC_2_legend,fill=heatcol_CC2,cex=pcex2,bty=\"n\")
    ## BP-level2
    par(fig=c(x[1],x[2],y[5],y[6]),new=T)
    pie(BP_2_splice,labels=BP_2_label,col=heatcol_BP2,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    legend(\"topleft\",\"BP\",bty=\"n\",cex=2,text.col=\"#6B8E23\",text.font=2)
    par(fig=c(x[3],x[4],y[5],y[6]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=BP_2_legend,fill=heatcol_BP2,cex=pcex2,bty=\"n\")
    legend(\"topleft\",paste(\"level\",L1,sep=\"\"),bty=\"n\",cex=2,text.font=2)
    
    ## MF-level3
    par(fig=c(x[5],x[6],y[1],y[2]),new=T)
    pie(MF_3_splice,labels=MF_3_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    #symbols(0,0,circles=0.1,add=T,bg=\"red\")
    par(fig=c(x[7],x[8],y[1],y[2]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=MF_3_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    ## CC-level3
    par(fig=c(x[5],x[6],y[3],y[4]),new=T)
    pie(CC_3_splice,labels=CC_3_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    par(fig=c(x[7],x[8],y[3],y[4]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=CC_3_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    ## BP-level3
    par(fig=c(x[5],x[6],y[5],y[6]),new=T)
    pie(BP_3_splice,labels=BP_3_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    par(fig=c(x[7],x[8],y[5],y[6]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=BP_3_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    legend(\"topleft\",paste(\"level\",L2,sep=\"\"),bty=\"n\",cex=2,text.font=2)
    
    ## MF-level4
    par(fig=c(x[9],x[10],y[1],y[2]),new=T)
    pie(MF_4_splice,labels=MF_4_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    par(fig=c(x[11],x[12],y[1],y[2]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=MF_4_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    ## CC-level4
    par(fig=c(x[9],x[10],y[3],y[4]),new=T)
    pie(CC_4_splice,labels=CC_4_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    par(fig=c(x[11],x[12],y[3],y[4]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=CC_4_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    ## BP-level4
    par(fig=c(x[9],x[10],y[5],y[6]),new=T)
    pie(BP_4_splice,labels=BP_4_label,col=ColorList,clockwise=T,cex=pcex1,radius=0.85,border=\"white\")
    par(fig=c(x[11],x[12],y[5],y[6]),new=T)
    plot(1:2,xlab=\"\",ylab=\"\",axes=F,type=\"n\")
    legend(\"left\",legend=BP_4_legend,fill=ColorList,cex=pcex2,bty=\"n\")
    legend(\"topleft\",paste(\"level\",L3,sep=\"\"),bty=\"n\",cex=2,text.font=2)
    dev.off()    
}

################# function of plot bar #################
plot_bars<-function(input,width,height,pcex1,pcex2){
   
}


################# ploting pies and bars ###################
plot_pie<-\"$opts{pie}\"
plot_bar<-\"$opts{bar}\"
Level_1<-paste(\"$opts{i}\",\".level\",$opts{l1},\".txt\",sep=\"\")
Level_2<-paste(\"$opts{i}\",\".level\",$opts{l2},\".txt\",sep=\"\")
Level_3<-paste(\"$opts{i}\",\".level\",$opts{l3},\".txt\",sep=\"\")
pdfname<-paste(\"$opts{i}\",\".Level\",\"$opts{l1}\",\"$opts{l2}\",\"$opts{l3}\",\".pdf\",sep=\"\")
if(plot_pie==\"T\"){
    library(\"RColorBrewer\")
    library(\"plotrix\")
    plot_pies(L1=$opts{l1},L2=$opts{l2},L3=$opts{l3},l1=Level_1,l2=Level_2,l3=Level_3,pdfname=pdfname,width=$opts{pw},height=$opts{ph},pcex1=$opts{pc1},pcex2=$opts{pc2})
         
}
if(plot_bar==\"T\"){
    plot_bars()    
}

";

system ("R --restore --no-save < $opts{o}.r");
system ("rm $opts{o}.r");








