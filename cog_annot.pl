#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 <COG.list> <COG classification> <COG.annot.xls>\n" unless @ARGV==3;
my ($list, $all, $out) = @ARGV;
open ALL, "<$all" || die $!;
my %h;
while (<ALL>) {
	chomp;
	my @lines = split(/\t/, $_);
	my $cog = $lines[0];
	$h{$cog} = $_;
}
close ALL;

open OUT, ">$out" || die $!;
print OUT "Accession\tCOG\tFunction\tName\n";
open LIST, "<$list" || die $!;
while (<LIST>) {
	chomp;
	my @lines = split (/\t/, $_);
	my $id = $lines[0];
	my $cog = $lines[1];
	if(exists $h{$cog}) {
		print OUT $id."\t".$h{$cog}."\n";
	}
}
close OUT;
close LIST;

