#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 protein.list full_uniprot_2_string.04_2015.tsv id2string.list
Note: protein.list should be Uniprot accession or entry
" unless @ARGV==3;

my ($list, $all, $out) = @ARGV;

open ALL, "<$all" || die $!;
my %idmapping;
my $last_id1 = "aaa";

my $header = <ALL>;
while (<ALL>) {
	chomp;
	my @line = split(/\t/, $_);
	my ($taxa, $id, $id2) = @line[0..2];
	
	my @id1 = split(/\|/, $id);
	my $id1 = $id1[0];
	
	next if ($id1 eq $last_id1);
	$idmapping{$id1} = $taxa.".".$id2;
	
	$last_id1 = $id1;
}
close ALL;

open OUT, ">$out" || die $!;
open LIST, "<$list" || die $!;
while (<LIST>) {
	chomp;
	if(exists $idmapping{$_}) {
		print OUT $_."\t".$idmapping{$_}."\n";
	}
}
close OUT;
close LIST;

