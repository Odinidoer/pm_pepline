#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$Bin/;
use Cwd;
use Getopt::Long;

my $VERSION = "v4.201701"; 
my $time = localtime();
$time=~s/[\s+:]/_/g;
my $cwd = getcwd;
my $bin = $Bin."/bin";

my %opts;
GetOptions(\%opts,"dup=s","type=s", "orgC=s","up=s","mb=i","ml=i","clt=s","score=i","h!");
my $usage = <<"USAGE";
	Usage: perl proteomics_pipe_v4.pl [options]
	Description: This pipeline is used for analysis the proteomics project
	Version: $VERSION
	Options:
		-dup*	STRING	duplicate, T/F, yes:T or no:F
		-type*	STRING	species type, Prokaryotes: pro or eukaryotes: euk, default: pro
		-orgC	STRING	organism class: all/ animals/ plants/ bacteria, default:all
		-up	FLOAT	up regulate, default: 1.2
		-score	FLOAT	string combined score cutoff, default: 400
		-mb	FLOAT	white space of bottom, default:13
		-ml	FLOAT	white space of left, default:8
		-clt	STRING	cluster type(both,row,colum or none), default:both
		-h	Display this usage information
	Usage: $0 -dup T -type pro -orgC bacteria -up 1.20 -mb 13 -ml 8 -clt both -score 400
		$0 -dup T -type euk -orgC animals -up 1.20 -mb 13 -ml 8 -clt both -score 400
USAGE

die $usage if(!($opts{dup} && $opts{type}) || $opts{h});

$opts{up}=defined$opts{up}?$opts{up}:1.20;
$opts{score}=defined$opts{score}?$opts{score}:400;
$opts{mb}=defined$opts{mb}?$opts{mb}:13;
$opts{ml}=defined$opts{ml}?$opts{ml}:8;
$opts{clt}=defined$opts{clt}?$opts{clt}:"both";
$opts{orgC}=defined$opts{orgC}?$opts{orgC}:"all";

my $down = 1/$opts{up};

my $sh = "Analysis_$opts{up}_$time.sh";
my $outdir = "Analysis_$opts{up}_$time";
my $logdir = "log";
open SH,">$sh";
print SH "
#pipeline version: proteomics_pipe:$VERSION 
#PBS -l nodes=1:ppn=2
#PBS -l mem=5G
#PBS -q sg
cd \$PBS_O_WORKDIR
cd $cwd
mkdir $logdir

##GO and KEGG annotation
mkdir $outdir
mkdir $outdir/1.annot
cd $outdir/1.annot
cp $cwd/data/exp.list ./

#GO annotation
$bin/get_annot_list.pl exp.list $bin/idmapping/uniprot2GO.list GO.list
$bin/gene-ontology.pl -i GO.list -l 2 -list GO.level2.list > level2.go.txt
$bin/go-bar.pl -i level2.go.txt >> $cwd/$logdir/go_annot.log 2>&1
$bin/go_9pies.new.pl -i GO.list -o GO.list >> $cwd/$logdir/go_annot.log 2>&1

#KEGG pathway annotation
$bin/get_annot_list.pl exp.list $bin/idmapping/uniprot2ko.list pathway.txt
$bin/getKEGGfromBlast.pl -i pathway.txt -format kobas -o pathways -org ko -use_proxy -database /mnt/ilustre/users/bingxu.liu/workspace/annotation/db/kegg.db >> $cwd/$logdir/kegg_annot.log 2>&1
python $bin/keggclass.py -i pathways -orgC $opts{orgC} -annot T -enrich F >> $cwd/$logdir/kegg_annot.log 2>&1
rm -rf pathways; mv pathways.filter pathways
#$bin/st_kegg_layer.pl  pathways/pathway_table.xls >> $cwd/$logdir/kegg_annot.log 2>&1
$bin/kegg_brite.pl pathways/pathway_table.xls
$bin/kegg_pathway_top20bars.pl -t pathways/pathway_table.xls -o pathway >> $cwd/$logdir/kegg_annot.log 2>&1
";

my $cog;
if ($opts{type} eq "pro") {
	$cog = "COG";
	print SH "
#COG annotation
$bin/get_annot_list.pl exp.list $bin/idmapping/uniprot2COG.list COG.list
$bin/cog_annot.pl COG.list $bin/COG-KOG/COG-list.txt COG.annot.xls
$bin/cog_summary.pl -i COG.list
$bin/cog-bar.pl -i pic_COG
python $bin/KOGinfo.py -type COG -a COG.annot.xls -c COG.class.catalog.xls -o COG.classification.xls
";
}elsif($opts{type} eq "euk") {
	$cog = "KOG";
	print SH "
#KOG annotation
$bin/get_annot_list.pl exp.list $bin/idmapping/uniprot2KOG.list KOG.list
$bin/kog_annot.pl KOG.list $bin/COG-KOG/KOG-list.txt KOG.annot.xls
$bin/kog_summary.pl -i KOG.list
$bin/kog-bar.pl -i pic_KOG
python $bin/KOGinfo.py -type KOG -a KOG.annot.xls -c KOG.class.catalog.xls -o KOG.classification.xls
";
}

print SH "cd $cwd

##different expression analysis
cd $outdir
mkdir $cwd/$outdir/2.diffexp
cd $cwd/$outdir/2.diffexp
$bin/run_diff_analysis.pl -i $cwd/data/exp.txt -g $cwd/data/sample.config -u $opts{up}

cd $cwd/$outdir/2.diffexp
";

if ($opts{dup} eq "T") {
	print SH "
for i in `ls *.diff.exp.xls | awk -F \".diff.exp.xls\" '{print \$1}'`
do
	less \$i.diff.exp.xls | awk '{if(\$7==\"yes\"&&(\$4>$opts{up}||\$4<$down)){print \$1}}' > \$i.DE.list
	less \$i.diff.exp.xls | awk '{if(\$7==\"yes\"&&\$4>$opts{up}){print \$1}}' > \$i.up.list
	less \$i.diff.exp.xls | awk '{if(\$7==\"yes\"&&\$4<$down){print \$1}}' > \$i.down.list
	less \$i.diff.exp.xls | awk -F \"\\t\" '{printf \$1\"\\t\"\$5\"\\t\"\$6\"\\t\"; if(NR==1){print \"sig\"}else if(\$6<0.05){if(\$4>$opts{up}){printf \"up\"; if(\$6<0.01){print \"-p-0.01\"}else{print \"-p-0.05\"}}else if(\$4<$down){printf \"down\"; if(\$6<0.01){print \"-p-0.01\"}else{print \"-p-0.05\"}}else{print \"nosig\"}}else{print \"nosig\"}}' > \$i.volcano
	less \$i.diff.exp.xls | awk -F \"\\t\" '{printf \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"; if(NR==1){print \"sig\"}else if(\$6<0.05){if(\$4>$opts{up}){printf \"up\"; if(\$6<0.01){print \"-p-0.01\"}else{print \"-p-0.05\"}}else if(\$4<$down){printf \"down\"; if(\$6<0.01){print \"-p-0.01\"}else{print \"-p-0.05\"}}else{print \"nosig\"}}else{print \"nosig\"}}' > \$i.scatter
	$bin/Highchart.pl -yAxis_log -yAxis_min 0.0001 -yAxis_max 1 -type scatter -t \$i.volcano -yAxis_reversed -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the \"'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'\" -scatter_size 2,2,1,2,2 -scatter_symbol \"'triangle-down','triangle-down','diamond','triangle','triangle'\" >> $cwd/$logdir/diffexp.log 2>&1
	$bin/Highchart.pl -type scatter -t \$i.scatter -xAxis_log -yAxis_log -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the \"'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'\" -scatter_size 2,2,2,2,2 -scatter_symbol \"'triangle-down','triangle-down','diamond','triangle','triangle'\" >> $cwd/$logdir/diffexp.log 2>&1
done
";
}elsif($opts{dup} eq "F") {
	print SH "
for i in `ls *.diff.exp.xls | awk -F \".diff.exp.xls\" '{print \$1}'`
do
	less \$i.diff.exp.xls | awk '{if(NR>1&&(\$4>$opts{up}||\$4<$down)){print \$1}}' > \$i.DE.list
	less \$i.diff.exp.xls | awk '{if(NR>1&&\$4>$opts{up}){print \$1}}' >\$i.up.list
	less \$i.diff.exp.xls | awk '{if(NR>1&&\$4<$down){print \$1}}' > \$i.down.list
	$bin/plot_diffexp_scatter.pl -i \$i.diff.exp.xls -up $opts{up} -down $down >> $cwd/$logdir/diffexp.log 2>&1
done
";
}

print SH "
$bin/get_diff_up_down.py

cd $cwd

##heatmap including all samples
mkdir $outdir/5.heatmap
cd $outdir/5.heatmap
cat $cwd/$outdir/2.diffexp/*.DE.list | sort | uniq | awk \'BEGIN {print \"Accession\"} {print \$_}\' > All.DE.list
$bin/get_exp_from_list.pl All.DE.list $cwd/data/exp.txt all.diffexp.txt
$bin/plot_heatmap_trendline.pl -i all.diffexp.txt -clt $opts{clt} >> $cwd/$logdir/cluster.log 2>&1 2>&1
cd $cwd

##cluster including all samples without repeates
mkdir $outdir/5.cluster
cd $outdir/5.cluster
cat $cwd/$outdir/2.diffexp/*.DE.list | sort | uniq | awk \'BEGIN {print \"Accession\"} {print \$_}\' > All.DE.list
python $bin/get_groupmean.py -e $cwd/data/exp.txt -s $cwd/data/sample.config -o $cwd/data/group.exp.txt
$bin/get_exp_from_list.pl All.DE.list $cwd/data/group.exp.txt group.diffexp.txt
$bin/plot_heatmap_trendline.pl -i group.diffexp.txt -clt $opts{clt} >> $cwd/$logdir/cluster.log 2>&1 2>&1
cd $cwd

##different expression GO and KEGG annotation
mkdir $outdir/3.diffannot
cd $outdir/3.diffannot
ln -s $cwd/$outdir/2.diffexp/*.diff.exp.xls ./
cp $cwd/$outdir/1.annot/GO.list $cwd/$outdir/1.annot/pathway.txt ./
for i in *.diff.exp.xls
";

if ($opts{dup} eq "T") {
	print SH "
do
	less \$i | awk '{if(NR==1){print \$0}else if(\$7==\"yes\"&&(\$4>$opts{up}||\$4<$down)){print \$0}}' > \$i.new
	rm \$i
	mv \$i.new \$i
	less \$i | awk '{if(\$7==\"yes\"&&\$4>$opts{up}){print \$1}}' > \$i.up.list
	less \$i | awk '{if(\$7==\"yes\"&&\$4<$down){print \$1}}' > \$i.down.list
	less \$i | awk '{if(\$7==\"yes\"&&(\$4<$down||\$4>$opts{up})){print \$1}}' > \$i.DE.list
	$bin/get_annot_list.pl \$i.DE.list pathway.txt \$i.DE.pathway.txt
	$bin/go-multi-bars.pl -i GO.list -u \$i.up.list -d \$i.down.list >> $cwd/$logdir/go_diffannot.log 2>&1
	#python $bin/get_descrip_GOannot.py -i1 \$i.up.list-\$i.down.list-gobars.xls -i2 $cwd/data/protein.xls -o \$i.up-down.gobars.detail.xls >> $cwd/$logdir/go_diffannot.log 2>&1
    
	$bin/KEGG_add_exp2.pl -i \$i.DE.pathway.txt -exp \$i -o \$i.path -up $opts{up} -down $down -format kobas -org ko -use_proxy -database /mnt/ilustre/users/bingxu.liu/workspace/annotation/db/kegg.db >> $cwd/$logdir/kegg_diffannot.log 2>&1
	python $bin/keggclass.py -i \$i.path -orgC $opts{orgC} -annot T -enrich F >> $cwd/$logdir/kegg_diffannot.log 2>&1
	rm -rf \$i.path; mv \$i.path.filter \$i.path
done
";
}elsif ($opts{dup} eq "F") {
	print SH "
do
	less \$i | awk '{if(NR==1){print \$0}else if(\$4>$opts{up}||\$4<$down){print \$0}}' > \$i.new
	rm \$i
	mv \$i.new \$i
	less \$i | awk '{if(NR>1&&\$4>$opts{up}){print \$1}}' > \$i.up.list
	less \$i | awk '{if(NR>1&&\$4<$down){print \$1}}' > \$i.down.list
	less \$i | awk '{if(NR>1&&(\$4<$down||\$4>$opts{up})){print \$1}}' > \$i.DE.list
	$bin/get_annot_list.pl \$i.DE.list pathway.txt \$i.DE.pathway.txt
	$bin/go-multi-bars.pl -i GO.list -u \$i.up.list -d \$i.down.list >> $cwd/$logdir/go_diffannot.log 2>&1
	python $bin/get_descrip_GOannot.py -i1 \$i.up.list-\$i.down.list-gobars.xls -i2 $cwd/data/protein.xls -o \$i.up-down.gobars.detail.xls >> $cwd/$logdir/go_diffannot.log 2>&1

	$bin/KEGG_add_exp2.pl -i \$i.DE.pathway.txt -exp \$i -up $opts{up} -down $down -o \$i.path -format kobas -org ko -use_proxy -database /mnt/ilustre/users/bingxu.liu/workspace/annotation/db/kegg.db >> $cwd/$logdir/kegg_diffannot.log 2>&1
	python $bin/keggclass.py -i \$i.path -orgC $opts{orgC} -annot T -enrich F >> $cwd/$logdir/kegg_diffannot.log 2>&1
	rm -rf \$i.path; mv \$i.path.filter \$i.path
done	
";
}

print SH "
cd $cwd

##different expression GO and KEGG enrichment
mkdir $outdir/4.enrich
cd $outdir/4.enrich
ln -s /mnt/ilustre/users/bingxu.liu/workspace/RNA_Pipeline/RNA_database/gene_ontology.1_2.obo ./
ln -s $cwd/$outdir/2.diffexp/*.DE.list ./
cp $cwd/data/exp.list $cwd/$outdir/1.annot/GO.list ./
cp $cwd/$outdir/1.annot/pathway.txt ./
cat $bin/kegg_enrich.head.txt pathway.txt > pathway.new.txt
rm pathway.txt
mv pathway.new.txt pathway.txt

for i in *.DE.list
do
	$bin/goatools-master/scripts/find_enrichment.py \$i exp.list GO.list --alpha 0.05 --fdr > \$i.go_enrichment
	$bin/extract_goatools.pl -enrich \$i.go_enrichment -diff \$i -lib gene_ontology.1_2.obo >> $cwd/$logdir/go_enrich.log 2>&1
	$bin/enrich_barplot.pl -i \$i.enrichment.detail.xls -t GO -mb $opts{mb} -ml $opts{ml} >> $cwd/$logdir/go_enrich.log 2>&1
	
	$bin/diff_ko_select.pl -g \$i -k pathway.txt >> $cwd/$logdir/kegg_enrich.log 2>&1
	$bin/identify.py -f \$i.ko_annot -n BH -b pathway.txt -o \$i.kegg_enrichment.xls >> $cwd/$logdir/kegg_enrich.log 2>&1
	python $bin/keggclass.py -i \$i.kegg_enrichment.xls -orgC $opts{orgC} -annot F -enrich T >> $cwd/$logdir/kegg_enrich.log 2>&1
	rm -rf \$i.kegg_enrichment.xls; mv \$i.kegg_enrichment.xls.filter \$i.kegg_enrichment.xls
	$bin/enrich_barplot.pl -i \$i.kegg_enrichment.xls -t PATHWAY -mb $opts{mb} -ml $opts{ml} >> $cwd/$logdir/kegg_enrich.log 2>&1
done
cd $cwd




##GO_DAG  ##need run local##
mkdir $outdir/6.go_dag
cd $cwd/$outdir/6.go_dag
cp $cwd/$outdir/4.enrich/*.enrichment.detail.xls ./
for i in *.enrichment.detail.xls
do
	$bin/go_orthology_relation.pl -enrich \$i
done
cd $cwd

##Venn
mkdir $outdir/7.venn
cd $cwd/$outdir/7.venn
ln -s $cwd/$outdir/2.diffexp/*.DE.list ./
# $bin/Venn_rnaseq.pl -f -l -o
cd $cwd

##network by string database
mkdir $outdir/8.network
cd $outdir/8.network
ln -s $cwd/$outdir/2.diffexp/*.DE.list ./
for i in `ls *_vs_*.DE.list | awk -F \".DE.list\" \'{print \$1}\'`
do
	$bin/get_id2string.pl \$i.DE.list /mnt/ilustre/users/ting.kuang/database/string/all/full_uniprot_2_string.04_2015.tsv \$i.id2string.list
	$bin/choose_ppi.pl \$i.id2string.list /mnt/ilustre/users/ting.kuang/database/string/all/protein.links.v10.nr400.txt $opts{score} > \$i.network.txt
	$bin/get_exp_from_list.pl \$i.DE.list $cwd/$outdir/2.diffexp/\$i.diff.exp.xls \$i.DE.exp.xls
";

if ($opts{dup} eq "T") {
	print SH "\t$bin/igraph_plot_network.pl -iLinks \$i.network.txt -iExp \$i.DE.exp.xls -o \$i -n T -s T\n";
} elsif ($opts{dup} eq "F") {
	print SH "\t$bin/igraph_plot_network_no_replicate.pl -iLinks \$i.network.txt -iExp \$i.DE.exp.xls -o \$i -n T -s T\n";
}

print SH "
done
cd $cwd

##Result_Files
mkdir $outdir/Result_Files
cd $outdir/Result_Files
mkdir -p 0.Data 1.QualityControl
mkdir -p 2.Annotation/{2.1.GO/,2.2.KEGG/,2.3.$cog/}
mkdir -p 3.DiffExpAnalysis/{3.1.Statistics/{Venn/,Volcano/},3.2.Cluster/{Heatmap/,Trendline/},3.3.GO/{Annotation/,Enrichment/},3.4.KEGG/{Annotation/,Enrichment/},3.5.Network}
#mkdir -p 3.DiffExpAnalysis/3.3.GO/GO-DAG
";
if ($opts{dup} eq "T") {
        print SH "cp -r $cwd/$outdir/5.heatmap/*Heatmap.pdf 3.DiffExpAnalysis/3.2.Cluster/Heatmap;cp -r $cwd/$outdir/5.cluster/{*trendlines*.pdf,subclusters_*} 3.DiffExpAnalysis/3.2.Cluster/Trendline";
} elsif ($opts{dup} eq "F") {
        print SH "cp -r $cwd/$outdir/5.heatmap/*Heatmap.pdf 3.DiffExpAnalysis/3.2.Cluster/Heatmap;cp -r $cwd/$outdir/5.heatmap/{*trendlines*.pdf,subclusters_*} 3.DiffExpAnalysis/3.2.Cluster/Trendline";
}

print SH "
cp -r $cwd/$outdir/1.annot/{GO.list,GO.list.*.xls,*.Level234.pdf,level2.go.txt.pdf} 2.Annotation/2.1.GO
cp -r $cwd/$outdir/1.annot/{pathway.txt,pathway.top20.pdf,KEGG_brite.*,pathways} 2.Annotation/2.2.KEGG
cp -r $cwd/$outdir/1.annot/{$cog.list,$cog.classification.xls,$cog.class.catalog.xls,$cog.class.catalog.pdf} 2.Annotation/2.3.$cog

cp -r $cwd/$outdir/2.diffexp/all_diff_up_down.xls 3.DiffExpAnalysis/3.1.Statistics
cp -r $cwd/$outdir/2.diffexp/{*.diff.exp.xls,*_vs_*.pdf,*.up.list,*.down.list,*.DE.list} 3.DiffExpAnalysis/3.1.Statistics/Volcano
cp -r $cwd/$outdir/3.diffannot/{*gobars.pdf,*gobars.xls} 3.DiffExpAnalysis/3.3.GO/Annotation
cp -r $cwd/$outdir/4.enrich/{*.enrichment.detail.xls,*.go.pdf} 3.DiffExpAnalysis/3.3.GO/Enrichment
#cp -r $cwd/$outdir/6.go_dag/{*.png,*.svg} 3.DiffExpAnalysis/3.3.GO/GO-DAG
cp -r $cwd/$outdir/3.diffannot/*.diff.exp.xls.path 3.DiffExpAnalysis/3.4.KEGG/Annotation
rm 3.DiffExpAnalysis/3.4.KEGG/Annotation/*.diff.exp.xls.path/pic.m
cp -r $cwd/$outdir/4.enrich/{*.pathway.pdf,*.pathway.xls} 3.DiffExpAnalysis/3.4.KEGG/Enrichment

cp -r $cwd/$outdir/8.network/{*.network.pdf,*.network.txt} 3.DiffExpAnalysis/3.5.Network

cd 3.DiffExpAnalysis/3.1.Statistics/Volcano
cp $cwd/data/protein.xls ./
cp $cwd/$outdir/Result_Files/2.Annotation/2.3.$cog/$cog.list ./
for i in *_vs_*.diff.exp.xls
do
	$bin/extract_diff_detail.pl -dup $opts{dup} -iDiff \$i -iProtein protein.xls -iExp $cwd/data/exp.txt -iGO $cwd/$outdir/Result_Files/2.Annotation/2.1.GO/GO.list -iKO $cwd/$outdir/Result_Files/2.Annotation/2.2.KEGG/pathway.txt -iCOG $cog.list
done
rm ./protein.xls ./$cog.list

##zip the Result_Files
#cd $outdir
#tar czf Result_Files.tar.gz Result_Files
";

close SH;

`qsub $sh`;
