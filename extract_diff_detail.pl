#!/usr/bin/perl -w

use Getopt::Long;
my %opts;
GetOptions(\%opts,"dup=s","iDiff=s","iProtein=s","iExp=s","iGO=s","iKO=s","iCOG=s","h!");
my $usage = <<"USAGE";
        Program:        $0
        Version:        v1_20160824
        Contact:        gan.li\@majorbio.com
        Discription:    add detail annot information to exp.diff.xls file
        Usage: perl $0
                -dup*   STRING  duplicate, T/F, yes:T or no:F
                -iDiff      input file: *_vs_*.exp.xls, e.g. A_vs_B.diff.exp.xls
                -iProtein   protein summary files  e.g.LQ-1_ProteinSummary.txt,LQ-2_ProteinSummary.txt
                -iExp       exp.txt
                -iGO        GO.list
                -iKO        pathway.txt
                -iCOG       COG.list or KOG.list
                -h          Display this usage information
        Example:            perl $0 -iDiff A_vs_B.diff.exp.xls -dup T -iProtein LQ-1_ProteinSummary.txt,LQ-2_ProteinSummary.txt -iExp exp.txt -iGO GO.list -iKO pathway.txt -iCOG KOG.list
USAGE

die $usage if((! $opts{iDiff}) || $opts{h} );
$strOUT = "";
@OUTheader = ();

@fileOUT = split(/\./,$opts{iDiff});
$fileOUT =$fileOUT[0];
$fileOUT.=".diff.exp.detail.xls";
open(FILEdiff,"<",$opts{iDiff})||die"cannot open the file: $!\n";
open(FILEexp,"<",$opts{iExp})||die"cannot open the file: $!\n";
open(FILEgo,"<",$opts{iGO})||die"cannot open the file: $!\n";
open(FILEko,"<",$opts{iKO})||die"cannot open the file: $!\n";
open(FILE_og,"<",$opts{iCOG})||die"cannot open the file: $!\n";
my $firstlineDiff = readline FILEdiff;
#my $firstlineExp = readline FILEexp;
#chomp $firstlineExp;
chomp $firstlineDiff;
#my @headerExp = split (/\t/, $firstlineExp);
my @headerDiff = split (/\t/, $firstlineDiff);
push @OUTheader,$headerDiff[0];
push @OUTheader,"Description";
#shift @headerExp;
shift @headerDiff;
my @headerDiffFC = @headerDiff;
splice @headerDiffFC,2;
my $n=0;
my @listN = ();
if ($opts{dup} eq "T") {
my $firstlineExp = readline FILEexp;
chomp $firstlineExp;
my @headerExp = split (/\t/, $firstlineExp);
shift @headerExp;
foreach $herderExpSs (@headerExp)
{
    my @herderExpS = split(/_/,$herderExpSs);
    my $headerExpHeader = $herderExpS[0];
    foreach (@headerDiffFC)
    {
        if ($headerExpHeader==$_)
        {
            push @OUTheader,$herderExpSs;
            push @listN,$n;
        }
    }
    $n++;
}
}
foreach (@headerDiff)
{
    push @OUTheader,$_;
}
push @OUTheader,"GO";
push @OUTheader,"KO";
@header_OG = split(/\./,$opts{iCOG});
$header_OG = $header_OG[0];
push @OUTheader,$header_OG;
while(<FILEdiff>)
{
    chomp;
    my @linesDiff = split(/\t/, $_);
    $keyDiff = $linesDiff[0];
    shift @linesDiff;
    $valueDiff = "";
    foreach (@linesDiff)
    {
        $valueDiff.= $_."\t";
    }
    push @keys,$keyDiff;
    $hashDiff{$keyDiff} = $valueDiff;
}
if ($opts{dup} eq "T") {
while(<FILEexp>)
{
    chomp;
    my @linesExp = split(/\t/, $_);
    $keyExp = $linesExp[0];
    shift @linesExp;
    $valueExp = "";
    foreach (@listN)
    {
        $valueExp.= $linesExp[$_]."\t";
    }
    $hashExp{$keyExp} = $valueExp;
}
}
while(<FILEgo>)
{
    chomp;
    my @linesGO = split(/\t/, $_);
    $keyGO = $linesGO[0];
    $valueGO = $linesGO[1]."\t";
    $hashGO{$keyGO} = $valueGO;
}
while(<FILEko>)
{
    chomp;
    my @linesKO = split(/\t/, $_);
    $keyKO = $linesKO[0];
    $valueKO = $linesKO[1]."\t";
    $hashKO{$keyKO} = $valueKO;
}
while(<FILE_og>)
{
    chomp;
    my @lines_OG = split(/\t/, $_);
    $key_OG = $lines_OG[0];
    $value_OG = $lines_OG[1]."\n";
    $hash_OG{$key_OG} = $value_OG;
}
@fileName = split(/,/, $opts{iProtein});
foreach (@fileName)
{
    open(FILEname,"<",$_)||die"cannot open the file: $!\n";
    my $firstlineName = readline FILEname;
    chomp ($firstlineName);
    my @headerName = split (/\t/, $firstlineName);
    @markersName = ("Accession","Description","Name");
    @markersNameExist = ();
    foreach (@markersName)
    {
        for (my $i = 0; $i < @headerName; $i++)
        {
            if ($headerName[$i] eq $_)
            {
                push @markersNameExist,$_;
                $strs{$_} = $i;
            }
        }
    }
    while (<FILEname>)
    {
    	chomp;
        my @linesName = split(/\t/, $_);
        $lineName = $linesName[$strs{$markersNameExist[0]}];
        if ($lineName =~ /\|/)
        {
            @keyName = split(/\|/,$lineName);
            $keyName = $keyName[1];
        }
        else
        {
            $keyName = $lineName;
        }
        $valueName = $linesName[$strs{$markersNameExist[1]}]."\t";
        $hashName{$keyName} = $valueName;
    }
}
foreach (@OUTheader)
{
$strOUT .= "$_\t";
}
chop($strOUT);
$strOUT .= "\n";
foreach (@keys)
{
    $strOUT .= "$_\t";
    $strOUT .= (exists $hashName{$_})?"$hashName{$_}":"-\t";
    if ($opts{dup} eq "T") {
    $strOUT .= (exists $hashExp{$_})?"$hashExp{$_}":"-\t";
    }
    $strOUT .= (exists $hashDiff{$_})?"$hashDiff{$_}":"-\t";
    $strOUT .= (exists $hashGO{$_})?"$hashGO{$_}":"-\t";
    $strOUT .= (exists $hashKO{$_})?"$hashKO{$_}":"-\t";
    $strOUT .= (exists $hash_OG{$_})?"$hash_OG{$_}":"-\n";
}
open (OUT, ">" ,$fileOUT) || die $!;
print OUT $strOUT;
