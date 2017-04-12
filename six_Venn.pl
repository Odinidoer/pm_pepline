#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"f=s","l=s","o=s","w=i","h=i","ls=f","ns=f","type=s");

my $usage = <<"USAGE";
        Program : $0
        Discription: plot venn for differently expressed genes(six groups ).
        Usage:perl $0 [options]
                -f      files     a,b,c
                -l      string    labels x,y,z
		-o	string	out file prefix
                -w      int       image width
                -h      int       image height
		-type    output format "Classic" or "Edwards";defult "Classic"

USAGE
die $usage if (!($opts{f}&&$opts{l}&& $opts{o}));
$opts{w}=$opts{w}?$opts{w}:500;
$opts{h}=$opts{h}?$opts{h}:415;
$opts{type}=$opts{type}?$opts{type}:"Classic";

my @files = split /,/, $opts{f};
my @ids = split /,/, $opts{l};

if ($opts{type} eq "Edwards"){
open CMD,">cmd.r" ;
print CMD "library(Vennerable)
sample<-c(\"$ids[0]\",\"$ids[1]\",\"$ids[2]\",\"$ids[3]\",\"$ids[4]\",\"[5]\")
a<-read.delim(\"$files[0]\",header=FALSE,check.names=FALSE,sep=\"\t\")
b<-read.delim(\"$files[1]\",header=FALSE,check.names=FALSE,sep=\"\t\")
c<-read.delim(\"$files[2]\",header=FALSE,check.names=FALSE,sep=\"\t\")
d<-read.delim(\"$files[3]\",header=FALSE,check.names=FALSE,sep=\"\t\")
e<-read.delim(\"$files[4]\",header=FALSE,check.names=FALSE,sep=\"\t\")
f<-read.delim(\"$files[5]\",header=FALSE,check.names=FALSE,sep=\"\t\")
data<-as.list(a)
data[[2]]<-b[[1]]
data[[3]]<-c[[1]]
data[[4]]<-d[[1]]
data[[5]]<-e[[1]]
data[[6]]<-f[[1]]
names(data)<-sample
Vstem <- Venn(data)
pdf(\"$opts{o}.pdf\",height=8,width=8)
plot(Vstem,doWeight=FALSE)
dev.off()\n";
`R --restore --no-save < cmd.r`;
}else
{
my %hash;
for(my $i = 0; $i < @files; $i ++)
{
	open FA, $files[$i] or die $!;
	while(<FA>)
	{
		chomp;
		$hash{$_} = 0;
	}
	close FA;
}

my %res;
foreach(keys %hash)
{
	for(my $i = 0; $i < @files; $i ++)
	{
		my $flag = 0;
		open FA, $files[$i] or die $!;
		while(my $line = <FA>)
		{
			chomp($line);
			if($line eq $_)
			{
				push @{$res{$_}}, 1;
				$flag = 1;
				last;
			}
		}
		close FA;
		
		if($flag eq 1)
		{
			next;
		}else{
			push @{$res{$_}}, 0;
		}
	}
}

#my $out = join "-", @ids;
my $head = join "\t", @ids;
open OUT, "> $opts{o}.xls" or die $!;
print OUT "gene_id\t$head\n";
foreach my $g(keys %res)
{
	my $string = join "\t", @{$res{$g}};
	print OUT "$g\t$string\n";
}
close OUT;

open FC,"$opts{o}.xls" or die $!;
my (%count,%all,%spe);
my $b = <FC>;
chop $b;
$all{all} = 0;
my @first = split /\t/,$b;
while (<FC>){
    chomp;
    my @line1 = split;
    my $sum = $line1[1]+$line1[2]+$line1[3]+$line1[4]+$line1[5]+$line1[6];
    if ($sum == 6){$all{all}++}
    for (my $t=1;$t<=6;$t++){
        if (($line1[$t] == 1) && ($sum == 1)){
        $spe{$first[$t]}++;
        }
    }
    for (my $j=1;$j<6;$j++){
        for (my $k=$j+1;$k<=6;$k++){
            if (($line1[$j] + $line1[$k] == 2)&&($sum == 2)){
            ($count{$first[$j]}{$first[$k]})++;
            }
        }
    }
}

close FC;

open FD,">$opts{o}.count.xls";
print FD "group\tnum\n";
print FD "all\t$all{all}\n";

##export specific num
for my $kk (1..6){
    print FD "$first[$kk]\t$spe{$first[$kk]}\n";
}

##if 2sample interation is null,def 0
for (my $j=0;$j<5;$j++){
        for (my $k=$j+1;$k<=5;$k++){
           if (exists $count{$ids[$j]}{$ids[$k]}){next;}
	   else {$count{$ids[$j]}{$ids[$k]} = 0;}
	}
}

##export 2sample interation
foreach my $key1(keys %count){
    foreach my $key2 (keys %{$count{$key1}}){
        print FD "$key1;$key2\t$count{$key1}{$key2}\n";
        }
}
close FD;

open SVG,">$opts{o}.svg"or die $!;
print SVG "
<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"$opts{w}\" height=\"$opts{h}\">
<defs></defs>
<g><path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(0,102,0)\" stroke=\"none\" d=\"M 0 11 L 254 160 L 174 235\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(90,155,212)\" stroke=\"none\" d=\"M 188 0 L 134 242 L 236 202\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(241,90,96)\" stroke=\"none\" d=\"M 338 52 L 135 123 L 191 242\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(207,207,27)\" stroke=\"none\" d=\"M 500 260 L 163 117 L 134 219\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(255,117,0)\" stroke=\"none\" d=\"M 250 415 L 133 150 L 203 67\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"none\"></path><g><path fill=\"rgb(192,152,83)\" stroke=\"none\" d=\"M 11 307 L 263 81 L 214 220\" fill-opacity=\"0.5\"></path></g>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 140 80 L 166 110\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 180 L 145 185\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 180 L 65 175\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 200 L 142 190\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 200 L 65 195\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 230 80 L 212 115\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 295 145 L 235 180\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 275 270 L 193 233\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 220 L 140 205\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 75 220 L 65 215\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<path fill=\"none\" stroke=\"rgb(0,0,0)\" d=\"M 150 270 L 183 230\" stroke-opacity=\"0.1\" stroke-miterlimit=\"10\"></path>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"99\" y=\"104\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[0]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"190\" y=\"64\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[1]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"281\" y=\"94\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[2]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"321\" y=\"219\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[3]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"213\" y=\"286\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[4]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"135\" y=\"74\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[0]}{$ids[1]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"145\" y=\"130\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[0]}{$ids[2]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"58\" y=\"171\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[0]}{$ids[3]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"143\" y=\"159\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[0]}{$ids[4]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"235\" y=\"74\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[1]}{$ids[2]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"230\" y=\"204\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[1]}{$ids[3]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"196\" y=\"96\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[1]}{$ids[4]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"305\" y=\"146\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[2]}{$ids[3]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"278\" y=\"282\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[2]}{$ids[4]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"217\" y=\"225\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[3]}{$ids[4]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"103\" y=\"254\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$spe{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"150\" y=\"282\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[4]}{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"56\" y=\"211\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[3]}{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"235\" y=\"127\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[2]}{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"150\" y=\"232\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[1]}{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"58\" y=\"191\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$count{$ids[0]}{$ids[5]}</text>
<text fill=\"rgb(0, 0, 0)\" fill-opacity=\"1\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"700\" text-decoration=\"none\" x=\"190\" y=\"184\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$all{all}</text>
<text fill=\"rgb(0, 102, 0)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"55\" y=\"24\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[0]</text>
<text fill=\"rgb(90, 155, 212)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"220\" y=\"19\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[1]</text>
<text fill=\"rgb(241, 90, 96)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"355\" y=\"74\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[2]</text>
<text fill=\"rgb(207, 207, 27)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"430\" y=\"214\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[3]</text>
<text fill=\"rgb(255, 117, 0)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"275\" y=\"399\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[4]</text>
<text fill=\"rgb(192, 152, 83)\" fill-opacity=\"0.5\" stroke=\"none\" font-family=\"Arial\" font-size=\"14px\" pt=\"\" font-style=\"normal\" font-weight=\"400\" text-decoration=\"none\" x=\"50\" y=\"314\" text-anchor=\"middle\" dominant-baseline=\"alphabetic\">$ids[5]</text>
</g></svg> \n";
close SVG;
`/usr/bin/convert $opts{o}.svg $opts{o}.png`;
}


