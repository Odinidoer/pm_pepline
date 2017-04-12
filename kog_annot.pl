#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 <KOG.list> <KOG classification> <KOG.annot.xls>\n" unless @ARGV==3;
my ($list, $all, $out) = @ARGV;
open ALL, "<$all" || die $!;
my %h;
while (<ALL>) {
	chomp;
	my @lines = split(/\t/, $_);
	my $kog = $lines[0];
	$h{$kog} = $_;
}
close ALL;

open OUT, ">$out" || die $!;
print OUT "Accession\tKOG\tFunction\tName\n";
open LIST, "<$list" || die $!;
while (<LIST>) {
	chomp;
	my @lines = split (/\t/, $_);
	my $id = $lines[0];
	my $kog = $lines[1];
	if(exists $h{$kog}) {
		print OUT $id."\t".$h{$kog}."\n";
	}
}
close OUT;
close LIST;

