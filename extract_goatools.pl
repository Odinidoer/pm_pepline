#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Getopt::Long;
use List::Util qw ( sum);
my %opts;
my $VERSION = "1.0";

GetOptions(
	\%opts, "enrich=s", "diff=s", "lib=s", "help!"
);

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : quan.guo\@majorbio.com
       Usage :perl $0 [options]
                -enrich		enrichment result of goatools
				-diff		diff expression protein list
                -lib		listfile  listall genenames for all the GO terms of this level 
                -help		Display this usage information
                * 		must be given Argument 
       example:gene-ontology.pl -i go.txt -l 2 
USAGE

die $usage
  if ( !( $opts{enrich} && $opts{diff} ) || $opts{help} );
my $lib	= $opts{lib}  if ( $opts{lib} );


my %diff;
my %lib;
open( DIFF, "<$opts{diff}" ) || die "Can't open $opts{diff} or file not exsit!\n";
while (<DIFF>){
	chomp;
	$diff{$_} = 1;
}
close DIFF;

open( LIB, "<$opts{lib}") || die "Can't open $opts{lib} or file not exsit!\n";
my $go_name = "";
while (<LIB>){
	chomp;
	if ( $_ =~ /id: (GO:\d+)/ ){
		$go_name = $1;
	}elsif( $_ =~ /namespace: (.*)$/ ){
		$lib{$go_name} = $1;
	}
}
close LIB;

open( FILE, "<$opts{enrich}" ) || die "Can't open $opts{enrich} or file not exsit!\n";
open( OUT, ">$opts{diff}.enrichment.detail.xls" ) || die "Can't open $opts{diff}.enrichment.detail.xls or file not exsit!\n";
open( OUT2, ">$opts{diff}.enrichment.xls" ) || die "Can't open $opts{diff}.enrichment.xls or file not exsit!\n";

my %go_all;
my %go_diff;

while (<FILE>) {
	chomp;
	my @name = split(/\t/, $_);
	$name[0] =~ s/\s.*//g;
	if($name[1] eq "p"){
		next;
	}

	
	if ($#name == 1){
		if (exists $go_all{$name[1]}){
			$go_all{$name[1]} .= ";".$name[0]
		}else{
			$go_all{$name[1]} = $name[0];
		}
		if (exists $diff{$name[0]} ){
			if (exists $go_diff{$name[1]}){
				$go_diff{$name[1]} .= ";".$name[0]
			}else{
				$go_diff{$name[1]} = $name[0];
			}			
		}
	}elsif($name[0] eq "id"){
		print OUT2 $_."\n";
		my $title = $name[0]."\t".$name[1]."\t".$name[2]."\t".$name[3]."\t".$name[4]."\t".$name[5]."\t".$name[6];
		$title .= "\ttype\tdiff_proteins\n";
		print OUT $title;
	}else{
		print OUT2 $_."\n";
		#my $line = join ("\t",@name);
		if($name[6] > 1){
			$name[6] = 1;
		}
		if($name[1] eq "p" || $name[5] > 0.05){
			next;
		}
		my $line = $name[0]."\t".$name[1]."\t".$name[2]."\t".$name[3]."\t".$name[4]."\t".$name[5]."\t".$name[6];
		#print $name[0];
		my $type = $lib{$name[0]};
		my $proteins = $go_all{$name[0]};
		
		my @protein_list = split(/;/,$proteins);
		my %protein_h;
		@protein_h{@protein_list}=();
		my @protein_new = sort keys %protein_h;
		my $proteins_new = join(";",@protein_new);
		
		my $diff_proteins = $go_diff{$name[0]};
		my @diff_protein_list = split(/;/,$diff_proteins);
		my %diff_protein_h;
		@diff_protein_h{@diff_protein_list}=();
		my @diff_protein_new = sort keys %diff_protein_h;
		my $diff_proteins_new = join(";",@diff_protein_new);
		
		$line .= "\t".$type."\t".$diff_proteins_new."\n";
		print OUT $line;
	}

}
close FILE;
close OUT;
close OUT2;
