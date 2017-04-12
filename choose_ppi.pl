#!/usr/bin/perl

use strict;
use warnings;
@ARGV==3 || die "Usage:perl $0 <id2string.list> <links.txt> <score> > network.txt
###id2string.list  (no header)####
#C4G7Z7  626523.GCWU000342_00163
#C4G7Z8  626523.GCWU000342_00164
#C4G815  626523.GCWU000342_00186
###links.txt (with header)###
#protein1 protein2 combined_score
#626523.GCWU000342_00163	626523.GCWU000342_00164	613
#626523.GCWU000342_00121	626523.GCWU000342_00123	358
###score###
#700 or 900
";

my $score = $ARGV[2];


my (%idmapping, $taxaid);
open MAP,"<$ARGV[0]";

while (<MAP>) {
	chomp;
	my ($id1, $id2) = split(/\t/);
	$idmapping{$id2} = $id1;
	$taxaid = (split(/\./, $id2))[0];
}

close MAP;

open ACT,"<$ARGV[1]";
my $title= <ACT>;
print "Protein1\tProtein2\tCombined_score\n";
while(<ACT>){
	chomp;
	next if (!/^$taxaid\./);
	my @b=split /\t/;
	next if ($b[-1] < $score);
	
	if ((exists $idmapping{$b[0]}) and (exists $idmapping{$b[1]})){
		#print "$b[0]\t$b[1]\t$b[2]\n";
		print "$idmapping{$b[0]}\t$idmapping{$b[1]}\t$b[2]\n";
	}
}
close ACT;
