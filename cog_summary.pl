#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Bio::SearchIO;
use Getopt::Long;
use FindBin qw($RealBin);
use DBI qw(:sql_types);
my %opts;
my $VERSION="v1";
GetOptions( \%opts,"i=s","db=s","h!");
my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Discription:	COG summary and plot cog function classification bars
       Usage :perl $0 [options]
			-i*	input file, COG.list
			-db	COG database,defualt: /mnt/ilustre/app/rna/database/COG/cog.db
			-h	display help message
			
        exmaple:perl $0 -i COG.list
USAGE

die $usage if((! $opts{i}) || $opts{h} );

unless($opts{db}){
	$opts{db}="/mnt/ilustre/app/rna/database/COG/cog.db";
	$opts{db}="/mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/cog.db" unless (-e "/mnt/ilustre/app/rna/database/COG/cog.db");
}

my $dbh = DBI->connect("dbi:SQLite:dbname=$opts{db}","","",{AutoCommit => 1});

my @seqs;
my %cogs;
open IN, "<$opts{i}" || die $!;
while (<IN>) {
	chomp;
	my @id_cogs = split (/\t/, $_);
	my $id = $id_cogs[0];
	my $cog = $id_cogs[1];
	push @seqs, $id;
	$cogs{$id} = $cog;
}
close IN;

my %sumary;
my $general=&getGeneral();
foreach my $acc (keys(%cogs)){
	my $cat = &getCat($cogs{$acc});
	if($cat){
		foreach my $aa (@$cat){
			$sumary{$general->{$aa}{type}}{$aa}{"COG"}{$acc}=1;	
		}
	}
}

open(SUMARY, ">COG.class.catalog.xls") || die $!;
print SUMARY "#Type\tFunctional_categories\tCOG\n";
my @types=('INFORMATION STORAGE AND PROCESSING','CELLULAR PROCESSES AND SIGNALING','METABOLISM','POORLY CHARACTERIZED');
foreach my $t (@types){
	foreach my $c (sort keys %{$sumary{$t}}){
		print SUMARY "$t\t[$c] $general->{$c}{'name'}";			
		print SUMARY "\t".scalar(keys %{$sumary{$t}{$c}{"COG"}});
		print SUMARY "\n";
	}
}
close SUMARY;

warn("Generating image file ...\n");
open(PIC, ">pic_COG") || die $!;
foreach my $t (@types){
	foreach my $c (sort keys %{$sumary{$t}}){
		my $tile=$general->{$c}{'name'};
		$tile=~/^\s*(\S+.*\S+)\s*$/;
		print PIC "$1\t".scalar(keys(%{$sumary{$t}{$c}{COG}}))."\n";
	}
}
close PIC;


sub getCat(){
	my $cog=shift;
	my @catelog;
	my $categeries=$dbh->prepare(<<SQL
select `type` from categeries where id=?;
SQL
	);
	$categeries->execute($cog);
	while(my $res1=$categeries->fetch()){
		push(@catelog,$res1->[0]);
	}
	return \@catelog;
}

sub getGeneral(){
	my $query=$dbh->prepare(<<SQL
select `id`,`name`,`type` from `general`;
SQL
			  );
	$query->execute();
	my $l = $query->fetchall_hashref('id');	
	return $l;
}


$dbh->disconnect();


