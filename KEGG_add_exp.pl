#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Bio::SearchIO;
use Getopt::Long;
use DBI qw(:sql_types);
use SOAP::Lite;
#use autodie;
use Try::Tiny;
use LWP::Simple;
use LWP::UserAgent;
use HTML::TreeBuilder;
#use URI::Escape;
use Math::Round qw(:all);
#use HTML::Template;
use HTML::Manipulator;





my %opts;
my $VERSION="2.0";
GetOptions( \%opts,"i=s", "format=s","o=s","maxEvalue=f","minIdentity=i","log=f","org=s","fresh!","exp=s","exptype=s","rank=i","database=s","QminCoverage=i","HminCoverage=i","use_proxy!","proxy_server=s","server=i","parse_id!", "h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : liubinxu
       Lastest modify:2013-06-15
       Discription:parse blast to genes databse result and get kegg pathway info and map
       				please install perl model: bioperl DBI DBD::SQLite SOAP::Lite autodie Try::Tiny
       Usage :perl $0 [options]
                -i*		blastn.out		blast to genes database output,can use wildcard character in bash,like '*_blast.out',but must use '' when using  wildcard character     
                -format		blastformat		the format of blast output
							kobas	kobas2 anntation file		wise       Genewise -genesf format
				-exp		expressfile
				-exptype		edgeR|cufflink
				-maxEvalue	1e-5			Max E value,default:1e-6
				-minIdentity	75			Min percent(Positives) identity over the alignment,default:75
				-QminCoverage	0			Min Converage percent of Query ,Suggest set to >70 when Query seq is protien
				-HminCoverage	30			Min Converage percent of Hit,Defualt:30
				-rank	10		rank cutoff for valid hits from BLAST result, Defualt:10
				-parse_id			parse_id from query description not the Query ID
				-o		dir			output dir,defualt: pathways                 
                -org		organism		organism name of three letters ,list in http://www.genome.jp/kegg/catalog/org_list.html ,like hsa
                								defuat:ko
                								also can use:map
				-log		2 or 10, default: 2			
                -fresh		fresh database from network
                -database	database path		default/mnt/lustre/database/kegg/kegg.db
                -use_proxy	whether use http proxy
				-proxy_server		http proxy server address default:http://101.168.10.1:8888/

                -h					Display this usage information
                * 					must be given Argument
                exmaple:perl $0 -i 'unfinish_*.out' -format blastxml -minIdentity 70
USAGE

die $usage if ((!$opts{i})||$opts{h});
$opts{format}=$opts{format}?$opts{format}:"kobas";
$opts{o}=$opts{o}?$opts{o}:"./pathways";
$opts{org}=$opts{org}?$opts{org}:"ko";
$opts{log}=$opts{log}?$opts{log}:2;
$opts{exptype}=$opts{exptype}?$opts{exptype}:"edgeR";
#$opts{database}=$opts{database}?$opts{database}:"/home/db/kegg/kegg.db";
unless($opts{database}){
	if(-f "/state/partition1/kegg/kegg.db"){
		$opts{database}="/state/partition1/kegg/kegg.db";
	}else{
		$opts{database}="/mnt/lustre/database/kegg/kegg.db";
	}
}
$opts{proxy_server}=$opts{proxy_server}?$opts{proxy_server}:"http://101.168.10.1:8888/";

$opts{maxEvalue}=$opts{maxEvalue}?$opts{maxEvalue}:"1e-6";
$opts{minIdentity}=$opts{minIdentity}?$opts{minIdentity}:75;
$opts{HminCoverage}=$opts{HminCoverage}?$opts{HminCoverage}:30;
$opts{rank}=$opts{rank}?$opts{rank}:"10";
unless(-f $opts{database}){
	warn("Database not exists,Create new ...\n");
	$opts{fresh}=1;
}

my $ua =LWP::UserAgent->new();

$ua->timeout(60);
$ua->proxy('http',$opts{proxy_server}) if($opts{use_proxy});

my $dbh = DBI->connect("dbi:SQLite:dbname=$opts{database}","","",{AutoCommit => 1});
my $check=$dbh->prepare("select count(*) from sqlite_master where type='table' and name='pathway_".$opts{org}."'");

$check->execute();
my @row_ary  = $check->fetchrow_array;
if ($row_ary[0]<=0){
	$opts{fresh}=1;
	warn("Local database has no info of this organism,getting from kegg network ...\n");
}

if($opts{format} eq 'kobas' && $opts{org} ne 'ko'){
	$check=$dbh->prepare("select count(*) from sqlite_master where type='table' and name='gene_pathway_".$opts{org}."'");

	$check->execute();
	my @row_ary  = $check->fetchrow_array;
	if ($row_ary[0]<=0){
		$opts{fresh}=1;
		warn("Local database has no info of this organism,getting from kegg network ...\n");
	}	
}

&freshdatabase($opts{org}) if($opts{fresh});

my %expression;
my %exp_label;
my %ko_exps;

open(EXP,"< $opts{exp}") || die "can not open $opts{exp}";
my $head=1;
while(<EXP>){
      chomp;
      my @line = split (/\t/,$_);
      my $add = 0;
      if($opts{exptype} eq "cufflink"){
	      $add = 1;
      }
      if($head == 1){
	      $head=0;
	      $exp_label{id}=$line[0];
	      $exp_label{count1}=$line[1+$add];
	      $exp_label{count2}=$line[2+$add];	      
      }else{
	    $expression{$line[0]}{count1} = $line[1+$add];
	    $expression{$line[0]}{count2} = $line[2+$add];
      }      
}
close EXP;
my @file= glob $opts{i};
warn("Input blast result files:\n");
warn(join("\n",@file)."\n");
#my %hash;
my $pathway = &getpathways($opts{org});
my $kos =&getpathwaykos($opts{org});

my %seqkos;

mkdir("$opts{o}","493") or die "Can't create dir at $opts{o}\n" unless( -e $opts{o});
open(KEGG, "> $opts{o}/kegg_table.xls") || die "Can't open $opts{o}/kegg_table.xls\n";
unless($opts{format} eq 'kobas'){
	print KEGG "Queryname\tHitname\tHit_discription\tevalue\tScore\ttopHSP_strand\tMax_identity\tQuery_length\ttopHSP_Query_converage\tHit_length\tmaxHSP_Hit_coverage\tkos\tecs\tpathway\n";
	foreach my $f (@file){
		warn("Parsing blast result file $f ...\n");
		my $searchio= Bio::SearchIO->new(-format => $opts{format},
									 -file => $f,
									 -best => 1,
									);
		while(my $result = $searchio->next_result){
				my $algorithm=$result->algorithm();
				die "Only support blastp and blastx result!\n" unless($algorithm=~/blastx|blastp/i);
				my $query_name=$result->query_name;
				if($opts{parse_id}||$query_name=~/^Query_\d+$/){
					$query_name=$result->query_description;
					$query_name=~/^\s*(\S+)/;
					$query_name=$1;				
				}else{
					$query_name=$result->query_name;
					$query_name=~/^\s*(\S+)/;
					$query_name=$1;
				}
				my $query_length=$result->query_length;
	
				my @quiery_ko;
			while(my $hit = $result->next_hit){
				last if $hit->rank() > $opts{rank};
				my $hit_length=$hit->length();
				my $score=$hit->score();
				my @paths;
				my @kos;
				my @ecs;
				my $hsp= $hit->hsp; #Bio::Search::HSP::HSPI
				my ($query_hsp_length,$hit_hsp_length);
				
					$query_hsp_length=$hsp->length('qeury');
					$hit_hsp_length=$hsp->length('hit');
					#print "$query_name\t$b\n";
	
				my ($query_coverage,$hit_coverage);
				$query_coverage=$query_hsp_length/$query_length;
				$hit_coverage=$hit_hsp_length/$hit_length;
				if($opts{'QminCoverage'}){
					next if $query_coverage <$opts{'QminCoverage'}/100;
				}
				if($opts{'HminCoverage'}){
					next if $hit_coverage <$opts{'HminCoverage'}/100;
				}
				if($opts{'maxEvalue'}){
					last if $hsp->evalue > $opts{'maxEvalue'};
				}
				my $identity=nearest(.01, $hsp->frac_conserved('total')*100);
	
				
				if($opts{'minIdentity'}){
					last if $identity < $opts{'minIdentity'};
				}
				$identity=$identity."%";
				
				#$hash{$result->query_name}{des}=$hit->description;
				my $des=$hit->description;
				#$hash{$result->query_name}{evalue}=$hsp->evalue;
				my $evalue=$hsp->evalue;
				#$hash{$result->query_name}{hitname}=$hit->name;
				my $hitname=$hit->name;
				#$hash{$result->query_name}{strand}=$hit->strand("query")==1?"+":"-";
				my $strand;
				if($algorithm=~/blastx/i){
					$strand=$hit->strand("query")==1?"+":"-";
				}else{
					$strand=" ";
				}
				
				
				while($des =~ /\s+(K\d{4,6})\s+/g ){
					#$hash{$result->query_name}{ko}=$1;
					my $ko=$1;
					push(@kos,$ko);	
					#print "$ko\n";			
					if(exists($kos->{"ko:".$ko})){
						foreach my $p (keys(%{$kos->{"ko:".$ko}})){
							#print "$p\n";	
							push(@paths,$p);
							push(@{$pathway->{$p}{'kos'}},"ko:".$ko);
							push(@{$pathway->{$p}{'seqs'}},$query_name);					
						}
					}else{
						warn("warn:$ko is not in all pathways in this organism $opts{org},if your database is newest, this is ok ... \n");
					}
				}
				
				while($des =~ /[\[\(](EC:[\d\.\-\s\,]+)[\]\)]/g){
					#$hash{$result->query_name}{ec}=$1;
					my $ec=$1;
					push(@ecs,$ec);
				}
				my $paths_ref=&uniq(\@paths);
				my $kos_ref=&uniq(\@kos);
				my $ecs_ref=&uniq(\@ecs);
				push(@quiery_ko,@$kos_ref);
				#print $result->query_name."\t".$hit->name."\t".$hash{$result->query_name}{ko}."\t".$hash{$result->query_name}{ec}."\t".$hit->description."\t".$hash{$result->query_name}{strand}."\t".$hsp->strand('hit')."\t".$hsp->evalue."\n";
				print KEGG "$query_name\t$hitname\t$des\t$evalue\t$score\t$strand\t$identity\t$query_length\t".sprintf("%.2f",$query_coverage*100)."%"."\t$hit_length\t".sprintf("%.2f",$hit_coverage*100)."%"."\t".join(";",@$kos_ref)."\t".join(";",@$ecs_ref)."\t".join(";",@$paths_ref)."\n";
				#print "$query_name\t$hitname\t$des\t$evalue\t$strand\t$identity\t".join(";",@$kos_ref)."\t".join(";",@$ecs_ref)."\t".join(";",@$paths_ref)."\n";
						
			}
			$seqkos{$query_name}=&uniq(\@quiery_ko);
		}
	}
}else{
	print KEGG "Protein\tKo_id\tPaths\n";
	foreach my $f (@file){
		open ANNOT, "< $f" or die "Error:Cannot open file  $f : $! \n";
		my $line=<ANNOT>;
		while(<ANNOT>){
			chomp;
			next if(/^\s*#/);
			last if(/^\/\/\/\//);
			if(/^(\S*)\s+(.*)$/){
				my $protein=$1;
				my $ko=$2;
				next if $ko=~/^None/;	

				if(exists $expression{$protein}){
					if(exists $ko_exps{$ko}){
						$ko_exps{$ko}{gene}.= ";".$protein;
						$ko_exps{$ko}{exp1}+= $expression{$protein}{count1};
						$ko_exps{$ko}{exp2}+= $expression{$protein}{count2};					
					}else{
						$ko_exps{$ko}{gene} = $protein;
						$ko_exps{$ko}{exp1} = $expression{$protein}{count1};
						$ko_exps{$ko}{exp2} = $expression{$protein}{count2};
					}
				}
				
				
				#print $ko."\n";
				my @paths;
				my @kos;
				push(@kos,$ko);	
				$seqkos{$protein}=\@kos;
				my $koid;
				if($opts{format} eq 'kobas' && $opts{org} ne 'ko'){
					$koid=$ko;
				}else{
					$koid="ko:".$ko;
				}
				
				if(exists($kos->{$koid})){
					foreach my $p (keys(%{$kos->{$koid}})){
						#print "$p\n";	
						push(@paths,$p);
						push(@{$pathway->{$p}{'kos'}},$koid);
						push(@{$pathway->{$p}{'seqs'}},$protein);					
					}
				}else{
					warn("warn:$ko is not in all pathways in this organism $opts{org},if your database is newest, this is ok ... \n");
				}
				my $paths_ref=&uniq(\@paths);
				print KEGG "$protein\t$ko\t".join(";",@$paths_ref)."\n";
			}
		}
		close ANNOT;
	}
}
close KEGG;

#get ko color

open(KOEXP,"> ko_exp.xls") || die "can not open ko_exp.xls";

print  KOEXP "ko\t$exp_label{count1}\t$exp_label{count2}\tFC($exp_label{count1}/$exp_label{count2})\tProteins\n";
foreach(keys %ko_exps){
	my $ko = $_;
	#if ($ko_exps{$ko}{exp1} == 0){
	#	$ko_exps{$ko}{fc} = -9999;
	#}elsif($ko_exps{$ko}{exp2} == 0){
	#	$ko_exps{$ko}{fc} = 9999;
	#}else{
		$ko_exps{$ko}{fc} = ko_exps{$ko}{exp1}/$ko_exps{$ko}{exp2};
	#}
	
	print  KOEXP $ko."\t".$ko_exps{$ko}{exp1}."\t".$ko_exps{$ko}{exp2}."\t".$ko_exps{$ko}{fc}."\t".$ko_exps{$ko}{protein}."\n";
	
	
	my %back;
	my %font;
	$back{red}=0;
	$back{green}=255;
	$back{blue}=0;
	$font{red}=255;
	$font{green}=0;
	$font{blue}=0;  
	
	if($ko_exps{$ko}{fc}>=2 || $ko_exps{$ko}{fc}<=0.5){
		$font{red}=0;
		$font{green}=255;
	}
	
	if($ko_exps{$ko}{fc}>=1.3){
		$back{red}=int($ko_exps{$ko}{fc}*128+0.5);
		if ($back{red} >255){
			$back{red}=255;
		}
		$back{green}=255-$back{red};
	}else if($ko_exps{$ko}{fc}<=0.7) {
	        $back{blue}=int(-$ko_exps{$ko}{fc}*128+0.5);
		if ($back{blue} >255){
			$back{blue}=255;
		}
		$back{green}=255-$back{blue};
	}
	
	
	$back{red}=sprintf("%.2x",$back{red});
	$back{green}=sprintf("%.2x",$back{green});
	$back{blue}=sprintf("%.2x",$back{blue});
	$font{red}=sprintf("%.2x",$font{red});
	$font{green}=sprintf("%.2x",$font{green});
	$font{blue}=sprintf("%.2x",$font{blue});
	
	

	$ko_exps{$ko}{color} = "%09%23$back{red}$back{green}$back{blue},%23$font{red}$font{green}$font{blue}";
	#print KOEXP $ko_exps{$ko}{color}."\n";
	
}

close KOEXP;

open(BAR,">scal_plot.r") || die "Can't open scal_plot.r\n";
print BAR '
x=seq(-2,2,by=0.0078125)
cramp <-"blue-green-red"
colramp <-unlist(strsplit(cramp,"-"))
ramp <-colorRamp(as.vector(colramp))
heatcol <- paste(rgb(ramp(seq(0,1,length=511)),max=255),"E5",sep="")
png(filename="'.$opts{o}.'/color_scale.png",width=460,height=70)
par(mar = c(2,2,2,2))
image(x,0:1,z=matrix(x,513,1),col=heatcol,yaxt="n",xlab="",ylab="",axes=TRUE,add=FALSE,cex.lab=1.7)
title("'.$exp_label{count1}.'/'.$exp_label{count2}.'",cex=0.6)
dev.off()
';

close BAR;
system("R --no-save <  scal_plot.r");


open(PATHWATY,"> $opts{o}/pathway_table.xls" ) || die "Can't open $opts{o}/pathway_table.xls\n";
print PATHWATY "Pathway\tPathway_definition\tNumber_of_Proteins\tprotein_ko_list\tPathway_image_name\n";
warn("outputing Pathway table ...\n");

foreach my $p (keys(%$pathway)){
	next unless exists($pathway->{$p}{'kos'});
	my $kolist=&uniq(\@{$pathway->{$p}{'kos'}});
	my $pathfile=&MarkPathway($p,$kolist);
	#print $pathfile."\n";
	my $imgname=&getimgname($p);
	my $htmlfile=&getimgname1($p);	
	#$filepath="./".$filepath unless($filepath=~/^\//);
	warn("Geting pathway image from  $pathfile  ...\n");
	#getstore($pathfile,$filepath);
	&savekegg($pathfile,$imgname,$htmlfile);
	my $seqlist=&uniq(\@{$pathway->{$p}{'seqs'}});
	my $seq_ko_list;
	foreach my $n (@$seqlist){
		$seq_ko_list.=$n."(".join(",",@{$seqkos{$n}}).");";
	}
	print PATHWATY "$p\t".$pathway->{$p}{'definition'}."\t".scalar(@$seqlist)."\t".$seq_ko_list."\t".$imgname."\n"; 
}

close KOEXP;

close PATHWATY;
warn("All done!\n");
sub savekegg(){
	my $url=shift;
	my $imgfile=shift;
	my $htmlfile=shift;
	
	try{
	    my $response =$ua->get($url);
	    my $filepath=$opts{o}."/".$imgfile;
	    my $html;
	    if($response->is_success){
	      $html = $response->decoded_content;
	      $html=&formathtml($html);
 	   
	      my $imgurl;
	      my $datapage=HTML::TreeBuilder->new_from_content($html);
	      
	      my @data=$datapage->find_by_attribute("usemap","#mapdata"); 
	   		
	      $imgurl=$data[0]->attr("src") or return &savekegg($url,$imgfile,$htmlfile);
	      $datapage->delete();

	      getstore($imgurl,$filepath) or &savekegg($url,$imgfile,$htmlfile);
	   
	      $html=&formathtml1($html,$imgfile);
	   	 open (HTML,"> $opts{o}/$htmlfile") or die "Can't create file $opts{o}/$htmlfile\n";
 		 print HTML  $html;
 		 close HTML;
		 
	      my $r=$ua->mirror($imgurl,$filepath);
	      
              unless($r->is_success||$r->code eq '304'){
                       warn("Saveing image file error!:".$r->status_line." retrying get from $imgurl ... \n");
                       return &savekegg($url,$imgfile,$htmlfile);
              }
              if($r->code eq '304'){
                      # warn("It seems to image haven't being coloring,download uncolored image form $imgurl.This problem is form kegg server.\n");
              }

	    }else{
		warn($response->status_line." retrying ... \n");
		&savekegg($url,$imgfile,$htmlfile);
	    }
	}catch{
                warn("Server connection serious error:$_,Geting pathway image from  $url  ...\n");
                return &savekegg($url,$imgfile,$htmlfile);
        }	
}

sub formathtml1(){
	my $htm =shift;
	my $imgname=shift;	
	$htm =~ s/<img src\=\".*\" usemap\=\"#mapdata\" border\=\"0\" \/>/<img src\=\"$imgname\" usemap\=\"#mapdata\" border\=\"0\" \/>/g;
	return $htm;
}

sub formathtml(){
	my($htm)=@_;
	$htm =~ s/\"\//\"http\:\/\/www.kegg.jp\//g;
	$htm =~ s/\'\//\'http\:\/\/www.kegg.jp\//g;
	return $htm;
}

sub getimgname(){
	my $path=shift;
	my @a=split(":",$path);
	return $a[1].".png";
}

sub getimgname1(){
	my $path=shift;
	my @a=split(":",$path);
	return $a[1].".html";
}

sub MarkPathway(){
	my $pathway_id = shift;
	my $list=shift;		
	#$list = SOAP::Data->type();
	#my $fg_list  = ['#ff0000'];
	#my $bg_list  = ['#ffff00'];
	
	$pathway_id =~ s/path://g ;
	
	my $color_result = "http://www.kegg.jp/kegg-bin/show_pathway?";
	$color_result .= $pathway_id;	
	
	
	#my $result = "http://rest.kegg.jp/";	
	#$result .= "get/";	
	#$result .= $pathway_id;
	foreach(@$list){
	  my $ko= $_;
	  $ko =~ s/ko://;	  
	  $color_result .='/'.$_.$ko_exps{$ko}{color};    
	}
	if (length ($color_result) > 4096){
	   $color_result = "http://www.kegg.jp/kegg-bin/show_pathway?".$pathway_id;
	}
	   #$color_result .= "/kgml"
	
	return $color_result;
	
	
	#$result .= "/kgml";
	
	
	
	
	
 	
#        try{
#		if($opts{server}==1){
#			 $result= $serv->get_html_of_marked_pathway_by_objects($pathway_id,$list);
#		}else{
#			 $result= $serv->KeggGetHtmlOfMarkedPathwayByObjects($pathway_id,$list);
#		}
#
#		if($result){
#			return $result;
#		}else{
#			warn("Server return error,retrying getting $pathway_id ...\n");
#			undef $serv;
#			$service= &serverref();
#			&MarkPathway($pathway_id,$list,$service);
#		}
#	}catch{
#		warn("Server connection serious error:$_,retrying getting $pathway_id ...\n");
#		$service= &serverref();
#		&MarkPathway($pathway_id,$list,$service);
#	}	
}

sub getpathways(){
	my $org=shift;
	my $pw=$dbh->prepare(<<SQL
select class,definition from pathway_$org;
SQL
			  );
	$pw->execute();
	my $ref = $pw->fetchall_hashref('class');
	$pw->finish;
	return $ref;
}

sub getpathwaykos(){
	my $org=shift;
	my %kolist;
	my $mm;
	if($opts{format} eq 'kobas' && $opts{org} ne 'ko'){
		$mm=$dbh->prepare(<<SQL
	select * from gene_pathway_$org;
SQL
			  );
	}else{
		$mm=$dbh->prepare(<<SQL
	select * from ko_pathway_$org;
SQL
			  );
	}
	
	$mm->execute();
	my $ref = $mm->fetchall_hashref('id');
	foreach my $ids ( keys(%$ref) ) {
			#push(@{$list->{$ref->{$id}->{'ko'}}},$ref->{$id}->{'pathway'});
			#print $ref->{$ids}->{'ko'}."\t".$ref->{$ids}->{'pathway'}." $ids\n";
			my $k=$ref->{$ids}->{'ko'};
			my $p=$ref->{$ids}->{'pathway'};
			$kolist{$k}{$p}=1;
	}
	$mm->finish;
	return \%kolist;
}


sub freshdatabase(){
	my $org=shift;
	warn("Freshing database from kegg netwok,please wating ...\n");
	$dbh->do(<<SQL
	drop table if exists pathway_$org;
SQL
 		 );
	$dbh->do(<<SQL
		 CREATE TABLE  pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 class  varchar(50) NOT NULL,
			 definition	 varchar(10) NOT NULL
		 );
SQL
 );
 	$dbh->do(<<SQL
	CREATE INDEX IF NOT EXISTS i_pathway_$org\_class ON pathway_$org(class);
SQL
 		 );
 	my $insert;
    if($opts{format} eq 'kobas' && $opts{org} ne 'ko'){
    	$dbh->do(<<SQL
	drop table if exists gene_pathway_$org;
SQL
 		 );
 	$dbh->do(<<SQL
	 CREATE TABLE gene_pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 pathway  varchar(50) NOT NULL,
			 ko	varchar(10) NOT NULL
	);
SQL
 		 );	 
	 $dbh->do(<<SQL
	 CREATE INDEX IF NOT EXISTS i_gene_pathway_$org ON gene_pathway_$org(pathway);
SQL
 		 );	
    }else{
    	  	$dbh->do(<<SQL
	drop table if exists ko_pathway_$org;
SQL
 		 );
 	$dbh->do(<<SQL
	 CREATE TABLE ko_pathway_$org( 
 			 id  INTEGER PRIMARY KEY ASC, 
			 pathway  varchar(50) NOT NULL,
			 ko	varchar(10) NOT NULL
	);
SQL
 		 );	 
	 $dbh->do(<<SQL
	 CREATE INDEX IF NOT EXISTS i_ko_pathway_$org ON ko_pathway_$org(pathway);
SQL
 		 );		 
		
		 
 	#$dbh->commit;
 
    }
 		$insert = $dbh->prepare(<<SQL
INSERT INTO pathway_$org(class,definition) VALUES (?,?);
SQL
);
	warn("getting pathway list ....\n");
 	my $pathway=&listpathways($org);
 	foreach my $p (@$pathway){
 		my $n=$p->{'entry_id'};
 		my $m=$p->{'definition'};
 		#print "$n $m \n";
 		$insert->execute($n,$m);
 	}
	#$dbh->commit;
	if($opts{format} eq 'kobas' && $opts{org} ne 'ko'){
		$insert = $dbh->prepare(<<SQL
INSERT INTO gene_pathway_$org(pathway,ko) VALUES (?,?);
SQL
); 
	}else{
		$insert = $dbh->prepare(<<SQL
INSERT INTO ko_pathway_$org(pathway,ko) VALUES (?,?);
SQL
); 
	}
	
	warn("getting kos list for each pathway ....\n");
	
	foreach my $p (@$pathway){
		warn("getting kos list for $p->{'entry_id'} ....\n");		
 		my $kos=&listkos($p->{'entry_id'});
 		foreach my $x (@$kos){
 			$insert->execute($p->{'entry_id'},$x);	
 		}
 	}
 	#undef $service;
 	#$dbh->commit;
}


#########################################
#list all pathways for organism
sub listpathways(){
    my $organism = shift;
	try{
		my $response=$ua->get("http://rest.kegg.jp/list/pathway/$organism");
		if($response->is_success){
			my $result=$response->decoded_content;
			my @lines=split(/\n+/,$result);
			my @a;
			foreach(@lines){
				if(/^(.+)\t+(.*)$/){
					my %h=('entry_id'=>$1,'definition'=>$2);
					push(@a,\%h);
				}
			}
			return \@a;
		}else{
			warn "Server response error:".$response->status_line."\n";
			warn("Server return error,retrying getting $organism ...\n");
			return &listpathways($organism);
		}
	}catch{		
		warn("Server connection serious error:$_,retrying getting $organism ...\n");
		return &listpathways($organism);
	}	
}

#####
#get all kos for a pathway
sub listkos(){
	my $pathway_id=shift;

	my $response;
        try{
		if($opts{org} =~ /^ko$|^map$/i){
			$response=$ua->get("http://rest.kegg.jp/link/ko/$pathway_id");
		
		}else{
			$response=$ua->get("http://rest.kegg.jp/link/genes/$pathway_id");
		}
		
		if($response->is_success){
				my $result=$response->decoded_content;
				my @lines=split(/\n+/,$result);
				my @a;
				foreach(@lines){
					if(/^(.+)\t+(.*)$/){
						push(@a,$2);
					}
				}
				return \@a;
		}else{
				warn "Server response error:".$response->status_line."\n";
				warn("Server return error,retrying  getting $pathway_id...\n");
				return &listkos($pathway_id);
		}		
	}catch{
		warn("Server connection serious error:$_,retrying  getting $pathway_id...\n");
		return &listkos($pathway_id);
    }
}

sub uniq {
	my $array      = shift;
	my %hash       = map { $_ => 1 } @$array;
	my @uniq_array = sort( keys %hash );
	return \@uniq_array;
}

