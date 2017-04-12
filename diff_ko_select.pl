#!/usr/bin/perl -w
use strict;
use warnings;

my $usage = <<"USAGE";
		Program : $0
        Discription: search Ko id of  
        Usage:perl $0 -g protein_list -k ko_annotation
		ko_anonotation is output file of annotate.py
USAGE

if (!$ARGV[3]) {
	print $usage;
	exit;
}


my $protein_list;
my $pathway_annotation;

for (my $i = 0; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] =~ /-g/) {
		$protein_list = $ARGV[$i+1];
		$i++;
    }elsif($ARGV[$i] =~ /-k/){
		$pathway_annotation = $ARGV[$i+1];
		$i++;
    }
}

(open PROTEIN, "<$protein_list") || die $usage;
(open PATHWAY, "<$pathway_annotation") || die $usage;
(open OUT, ">$protein_list.ko_annot") || die $usage;

my %protein;
my %annotation;

while(<PROTEIN>){
	chomp;
	$annotation{$_} = 1;
}
close PROTEIN;

my $context = 0;

while(<PATHWAY>){
	chomp;
	
	if ($context == 0){
		print OUT $_."\n";
		if (/#QueryKO	ID|KO name|Hyperlink/){
			$context =1;
		}elsif(/#Query  Ko id|Ko name|hyperlink/){
			$context =1;
		}
	}else{
		if ($_ =~ /^([^\t]*)\t(K\d{5}).*$/){
			
			if ( exists $annotation{$1} ){
				print OUT $1."\t".$2."\n";
			}
		}elsif($_ =~ /^(.*)\t([a-z]{3}:[^|]*)|.*$/){
			if ( exists $annotation{$1} ){
				print OUT $1."\t".$2."\n";
			}
		}
	}

}
close PATHWAY;
close OUT;
