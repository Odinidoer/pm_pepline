#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: $0 [pathway_table.xls]\n" unless @ARGV==1;
my $pathway_table = $ARGV[0];

my $kegg_brite = "KEGG_brite";
my $layer = "/mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/kegg_layer.txt";
open LAYER, "<$layer" || die $!;
my (%layer);
while (<LAYER>) {
	chomp;
	my @layers = split(/\t/, $_);
	my $pathway = $layers[2];
	$layer{$pathway} = $_;
}
close LAYER;

open PATHWAY, "<$pathway_table" || die $!;
my %pathway;
while (<PATHWAY>) {
	chomp;
	my @lines = split(/\t/, $_);
	$pathway{$lines[1]} = join("\t", ($lines[0], $lines[2], $lines[3]));
}
close PATHWAY;

#obtain all brite level
open BRITE, ">$kegg_brite" || die $!;
print BRITE "#Brite_A\tBrite_B\tPathway_Name\tPathway\tNumber_of_Protein\tProtein_KO_list\n";
foreach my $pathway (keys %pathway) {
	if (exists $layer{$pathway}) {
		print BRITE "$layer{$pathway}\t$pathway{$pathway}\n";
	}
	
}
close BRITE;

`sort -k 1,2 $kegg_brite > $kegg_brite.C.xls`;
`rm $kegg_brite`;


open BRITE, "<$kegg_brite.C.xls" || die $!;
my (%brite_b, %brite_a, %b_belongs_a);
while (<BRITE>) {
	next if /^#/;
	chomp;
	my @lines = split(/\t/, $_);
	if (exists $brite_b{$lines[1]}) {
		$brite_b{$lines[1]} .= $lines[-1].";";
	} else {
		$brite_b{$lines[1]} = $lines[-1];
	}
	
	if (exists $brite_a{$lines[0]}) {
		$brite_a{$lines[0]} .= $lines[-1].";";
	} else {
		$brite_a{$lines[0]} = $lines[-1];
	}
	
	$b_belongs_a{$lines[1]} = $lines[0];
}
close BRITE;

#obtain brite A
open A, ">$kegg_brite.A.xls" || die $!;
print A "#Brite_A\tNumber_of_Protein\tProtein_KO_list\n";
foreach my $brite_a ( sort {$a cmp $b} keys %brite_a) {
	my @proteins = split(/;{1,2}/, $brite_a{$brite_a});
	my %hash;
	my @uniq_proteins = grep { ++$hash{$_} < 2 } @proteins; #non-rebundance proteins
	my $number = @uniq_proteins;
	my $uniq_proteins = join(";", @uniq_proteins);
	print A "$brite_a\t$number\t$uniq_proteins\n";
}
close A;


#obtain brite B
open B, ">$kegg_brite.B.xls" || die $!;
print B "#Brite_A\tBrite_B\tNumber_of_Protein\tProtein_KO_list\n";
foreach my $brite_b (sort {$b_belongs_a{$a} cmp $b_belongs_a{$b}} keys %b_belongs_a) {
	my @proteins = split(/;{1,2}/, $brite_b{$brite_b});
	my %hash;
	my @uniq_proteins = grep { ++$hash{$_} < 2 } @proteins; #non-rebundance proteins
	my $number = @uniq_proteins;
	my $uniq_proteins = join(";", @uniq_proteins);
	print B "$b_belongs_a{$brite_b}\t$brite_b\t$number\t$uniq_proteins\n";
}
close B;

open RCMD, ">kegg_brite.r" || die $!;
print RCMD "
library(ggplot2)
dat <- read.table(\"KEGG_brite.B.xls\", head=T, sep = \"\\t\", comment.char = \"\")
colnames(dat) <- c(\"brite_a\", \"brite_b\", \"number\", \"proteins\")
dat <- dat[order(dat[,1],dat[,2]),] #å¯¹ç¬¬ä¸€åˆ—å’Œç¬¬äºŒåˆ—æŽ’åºï¼Œå³typeå’Œterm

pdf(\"KEGG_brite.pdf\", w=12, h=10)

dat\$brite_b <- factor(dat\$brite_b, levels=dat\$brite_b, ordered=T)  #ggplot²»»á×Ô¶¯ÅÅÐò

p <- ggplot(dat, aes(y=number, x=brite_b, fill=brite_a, group=factor(1)))
p <- p + geom_bar(stat=\"identity\", show.legend = T)
p <- p + theme(axis.text.x=element_text(angle=0, vjust=0.5, color=\"black\"))
p <- p + geom_text(aes(label=number, hjust=-0.2))
p <- p + labs(x=\"KEGG classification\", y=\"Number of Proteins\")
p <- p + ylim(min(dat\$number, 0)*1.05, max(dat\$number)*1.05)
p <- p + coord_flip() #x,yÖá»¥»»
p <- p + labs(fill=\"Classification\") + guides(fill=guide_legend(reverse=TRUE))

p

dev.off()

";
close RCMD;

`Rscript kegg_brite.r`;
