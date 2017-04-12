#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions (\%opts,"f=s","l=s","o=s","w=i","h=i");

my $usage = <<"USAGE";
	Program : $0
	Discription: plot venn proportional to the number of corresponding set elements
	Usage:perl $0 [options]
		-f	files	a,b,c
		-l	string	labels x,y,z
		-o	string	out file prefix
		-w	int	image width
		-h	int	image height
USAGE
die $usage if (!($opts{f}&&$opts{l}&& $opts{o}));
$opts{w}=$opts{w}?$opts{w}:10;
$opts{h}=$opts{h}?$opts{h}:10;


my @files = split /,/, $opts{f};
my @ids = split /,/, $opts{l};
my %hash;
for(my $i = 0; $i < @files; $i ++){
	open FA, $files[$i] or die $!;
	while(<FA>){
		chomp;
		$hash{$_} = 0;
	}
	close FA;
}

my %res;
foreach(keys %hash){
	for(my $i = 0; $i < @files; $i ++){
		my $flag = 0;
		open FA, $files[$i] or die $!;
		while(my $line = <FA>){
			chomp $line;
			if($line eq $_){
				push @{$res{$_}}, 1;
				$flag = 1;
				last;
			}
		}
		close FA;
		
		if($flag eq 1){
			next;
		}else{
			push @{$res{$_}}, 0;
		}
	}
}

#my $out = join "-", @ids;
my $head = join "\t", @ids;
open OUT, "> $opts{o}.Venn.xls" or die $!;
print OUT "Accession\t$head\n";
foreach my $g(keys %res){
	my $string = join "\t", @{$res{$g}};
	print OUT "$g\t$string\n";
}
close OUT;

open RCMD, ">cmd.r";
print RCMD "

dat <- read.table(\"$opts{o}.Venn.xls\", header = T, sep = \"\\t\", row.names = 1, as.is = T)
library(\"Vennerable\")
probe <- rownames(dat)
#sample <- colnames(dat)
sample <- unlist(strsplit(\"$opts{l}\",\",\",fix=T))
data <- as.list(dat)
names(data) <- sample
for (i in 1:length(sample)){
	hit <- dat[[i]]
	out <- probe[hit>0]
	data[[i]] <- out
}
Vstem <- Venn(data)
pdf(\"$opts{o}.Venn.pdf\", w = $opts{w}, h = $opts{h})

#plot(Vstem, doWeights = T, show = list(Faces = T), doEuler = T, type = \"circles\")
V <- compute.Venn(Vstem, doWeights = T, type = \"circles\")  # less equal 3 sample
V <- compute.Venn(Vstem, doWeights = F, type=\"ellipses\")  # equal 4 sample
grid.newpage()
SetLabels <- VennGetSetLabels(V)
V <- VennSetSetLabels(V, SetLabels)
plot(V)

dev.off()

";

`R --restore --no-save < cmd.r`;
#`rm *.r`;

