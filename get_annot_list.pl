#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 exp.list uniprot2go/ko/cog.list GO.list/pathway.txt/COG.list\n" unless @ARGV==3;
my ($list, $all, $out) = @ARGV;
open ALL, "<$all" || die $!;
my %idmapping;
while (<ALL>) {
	chomp;
	my ($id1, $id2) = split(/\t/, $_);
	$idmapping{$id1} = $id2;
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

