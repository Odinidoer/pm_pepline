#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","s=s","o=s","header!","h!");

my $usage = <<"USAGE";
	Program:	$0
	Version:	v1_20160427
	Contact:	ting.kuang\@majorbio.com
	Discription:	extract the special column (e.g. Accession 114/113 115/113 116/113 ...) from the database searching result (ProteinSummary.txt)
	Usage: perl $0
		-i	input file, ProteinSummary.txt
		-s	the special column, e.g. Accession,114:113,115:113,116:113
		-o	output file, the non-redundant expression file used to diffexp analysis (exp.xls) or all non-redundant proteins (protein.list)
		-header	add header for output file
		-h	Display this usage information
	Example:	perl $0 -i ProteinSummary.txt -s Accession,114/113,115/113,116/113 -o exp.txt -header
			perl $0 -i ProteinSummary.txt -s Accession -o protein.list
USAGE

die $usage if((! $opts{i}) || $opts{h} );

open OUT, ">$opts{o}" || die $!;
my @markers = split(/,/, $opts{s});
if ($opts{header}) {
	print OUT join("\t", @markers)."\n";
}

my %strs;
open SUM, "<$opts{i}" || die $!;
my $firstline = readline SUM;
chomp $firstline;
my @header = split (/\t/, $firstline);

foreach my $marker (@markers) {
	
	for (my $i = 0; $i < @header; $i ++) {
		if ($header[$i] eq $marker) {
			$strs{$marker} = $i;
		}
	}
}

my @list;
while (<SUM>) {
	chomp;
	my @lines = split(/\t/, $_);
	my @results;
	foreach my $marker (@markers) {
		if (defined $lines[$strs{$marker}]) {
			push @results, $lines[$strs{$marker}];
		} esle {
			next;
		}
	}

	if ($#results == $#markers) { #rm protein missing expression value
		my $accession = shift @results;
		if (grep { $_ eq $accession } @list) { #rm repeat protein
			next;
		} else {
			push @list, $accession;
			my $acc;
			my $name;
			next if ($accession =~ /^RRRRR/); #rm the re database result
			if ($accession =~ /.*\|(\S+)\|(\S+)/) {
				$acc = $1; #uniprot accession
				$name = $2; #uniprot name
				print OUT $acc."\t".join("\t", @results)."\n";
			} else {
				print OUT $accession."\t".join("\t", @results)."\n";
			}
		}
	}
}

close SUM;
close OUT;

