#!/usr/bin/perl -w

use strict;

die "usage:perl $0 <pathway table> \n" unless @ARGV == 1;

open A,"/mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/kegg_layer.txt";
open D,$ARGV[0];

my $base = $ARGV[0];
$base =~ s/.*\///g;
$base =~ s/\.*//g;
my %lay1;my %lay2;my %lay3;my %tar;my %tar2;my %kk;my %all;my %all2;my %st_gene;my %st_isof;

while(<A>){
    chomp;my @a=split /\t+/;$lay1{$a[2]}=[@a];
    $kk{$a[0]}=$.;
}
close A;

my %h;open OUT,">kegg_func.num";
while(<D>){
    chomp;next if(/^PathWay/);
    my @a=split /\t+/;my @b=split /;/,$a[3];my %h;
    for(my $i=0;$i < @b;$i++){
        my @c=split /[()]/,$b[$i];
        my @d=split /\_se/,$c[0];
        $h{$d[0]}=1;
        $all2{$lay1{$a[1]}[0]}{$d[0]}=1;
        $tar2{$lay1{$a[1]}[0]}{$lay1{$a[1]}[1]}{$d[0]}=1;
        $st_gene{$d[0]}++;$st_isof{$c[0]}++;

        $tar{$lay1{$a[1]}[0]}{$lay1{$a[1]}[1]}{$c[0]}++;
        $all{$lay1{$a[1]}[0]}{$c[0]}++;
    }
    my @site=sort keys %h;my $site=@site;
    print OUT "$site\t$a[2]\t$a[1]\n";
}
close D;close OUT;
my @num_gene=keys %st_gene;my $num_gene=$#num_gene;my @num_isof=keys %st_isof;my $num_isof=@num_isof;
open OUT,">to_figure_classify.txt";
open OUT2,">kegg_classification.txt";
foreach my $k1(sort {$kk{$a}<=>$kk{$b}} keys %kk){
    if(exists $all{$k1}){
        my @tmp=keys %{$all2{$k1}};my $num=@tmp;
        my @tmpp=keys %{$all{$k1}};my $numm=@tmpp;
        print OUT2 "$k1\t$numm($num)--$base(KO)\n";
    }
    else{
        print OUT2 "$k1\t0(0)\n";
    }
    foreach my $k2(sort keys %{$tar{$k1}}){
        print OUT2 " --$k2\t";
        if(exists $tar{$k1}{$k2}){
            my @tmp=keys %{$tar2{$k1}{$k2}};my $num=@tmp;
            my @tmpp=keys %{$tar{$k1}{$k2}};my $numm=@tmpp;
            print OUT2 "$numm($num)\n";
            print OUT "$k2";
            printf OUT "\t%d\t%d\t%.4f\t%d\t%.4f\n",$kk{$k1},$numm,$numm/$num_isof,$num,$num/$num_gene;
        }
        else{
            print OUT2 "0(0)\n";
        }
    }
}
close OUT;
`cp /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/kegg_classification_figure.m ./`;
`matlab < kegg_classification_figure.m`;
