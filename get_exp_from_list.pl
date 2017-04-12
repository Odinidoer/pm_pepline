#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 DE.list diff.exp.xls DE.exp.xls\n" unless @ARGV==3;
my ($list, $all, $out) = @ARGV;
open ALL, "<$all" || die $!;
my %h;
while (<ALL>) {
	chomp;
	my @lines = split(/\t/, $_);
	my $id = $lines[0];
	$h{$id} = $_;
}
close ALL;

open OUT, ">$out" || die $!;
open LIST, "<$list" || die $!;
while (<LIST>) {
	chomp;
	if(exists $h{$_}) {
		print OUT $h{$_}."\n";
	}
}
close OUT;
close LIST;

