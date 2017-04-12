#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","up=s","down=s","h!");
my $usage = <<"USAGE";
Usage: perl plot_diffexp_scatter.pl [options]
Description: This program is used to plot scatter diagram for no duplicate study.
Version: v2
Options:
	-i*	STRING	input file, e.g. *_vs_*.diff.exp.xls
	-up		FLOAT	up regulate, default: 1.2
	-down	FLOAT	down regulate, default:0.8
	-h	usage	information
Usage: $0 -i *_vs_*.diff.exp.xls -up 1.2 -down 0.8
USAGE

$opts{up}=defined$opts{up}?$opts{up}:1.2;
$opts{down}=defined$opts{down}?$opts{down}:0.8;

die $usage if(!$opts{i}||$opts{h});
my @outs = split(/\./,$opts{i});
my $out = $outs[0];

my $xlim;
if ($opts{up} > 1.5) {
	$xlim = 2;
} elsif ($opts{up} > 1.2 && $opts{up} <= 1.5){
	$xlim = 1.5;
} else {
	$xlim = 1.2;
}

open RCMD, ">$opts{i}.r" || die $!;
print RCMD "
require(ggplot2)
library(scales)
##change theme##
old_theme <- theme_update(
  axis.ticks=element_line(colour=\"black\"),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  panel.border=element_rect(fill='transparent', color='black'),
  axis.line=element_line(size=0.5)
)

dat <- read.table(\"$opts{i}\", header = T, sep = \"\t\")
df <- data.frame(c(dat[,1]),c(dat[,5]))
colnames(df) <- c(\"Accession\", \"log2FC\")
up <- log($opts{up}, 2)
down <- log($opts{down}, 2)
df.G <- subset(df, (log2FC < down )) #define Green
df.G <- cbind(df.G, rep(1, nrow(df.G)))
colnames(df.G)[3] <- \"Color\"
df.B <- subset(df, (log2FC >= down & log2FC <= up )) #define gray
df.B <- cbind(df.B, rep(2, nrow(df.B)))
colnames(df.B)[3] <- \"Color\"
df.R <- subset(df, (log2FC > up )) #define Red
df.R <- cbind(df.R, rep(3, nrow(df.R)))
colnames(df.R)[3] <- \"Color\"
df.t <- rbind(df.G, df.B, df.R)
df.t\$Color <- as.factor(df.t\$Color)
##Construct the plot object
p <- ggplot(data = df.t, aes(x = log2FC, y = as.factor(Accession), color = Color )) + labs(title = \"Protein ration distribution\",x = \"log2FC\", y = \"\")  + scale_y_discrete(\"\", breaks=NULL) +
  geom_point(size = 1.75,position='jitter',shape=17) + theme( legend.position = \"none\")+
  xlim(c(-$xlim, $xlim))  + scale_color_manual(values = c(\"green\", \"gray\", \"red\"))  + geom_vline(xintercept=down,col=\"gray\") + geom_vline(xintercept=up,col=\"gray\")

ggsave(\"$out.Scatter.pdf\",p)

";

`R --no-save --restore < $opts{i}.r`;
