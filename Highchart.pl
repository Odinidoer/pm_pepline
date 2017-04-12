#!/usr/bin/perl -w

use strict;
use FindBin qw($RealBin);
use Getopt::Long;
my %opts;
my $VERSION="2.0";
GetOptions( \%opts,"t=s","type=s","title=s", "legendx=i","legendy=i","width=i","height=i","nolegend!","marker=s","piemark=s","color_the=s","head=i","point_start=f","xAxis_lab=s","yAxis_lab=s","xAxis_lab_rota=s","yAxis_lab_rota=s","yAxis_reversed!","Axis_min=f","xAxis_log!","x_tickInterval=f","yAxis_log!","xAxis_min=f","xAxis_max=f","yAxis_min=f","yAxis_max=f","scatter_series=s","scatter_size=s","scatter_symbol=s","pdfcompress=s","jsdelay=f","h!");


my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact: liubinxu
       Lastest modify:2014-09-25
       Discription:
       SAMPLE  :perl $0 -t table.xls -type -type		area|arearange|areaspline|areasplinerange|bar|boxplot|bubble|column|columnrange|errorbar|funnel|gauge|heatmap|line|pie|pyramid|scatter|series|solidgauge|spline|waterfall
                
       Usage :	perl $0 [options]
				-t			table.xls 
				-type		area|arearange|areaspline|areasplinerange|bar|boxplot|bubble|column|columnrange|errorbar|funnel|gauge|heatmap|line|pie|pyramid|scatter|series|solidgauge|spline|waterfall
				-title		picure tittle
				-width		width of chart area
				-height		height of chart area
				
				-color_the	default|google|googlelight					
				-head		is table have header default 1 , 0 for no head title
				-xAxis_lab	x lable
				-yAxis_lab	y lable
				-yAxis_min  the min yAxis
				
				-xAxis_log 	logarithmic to xAxis
				-yAxis_log	logarithmic to yAxis
				-yAxis_reversed
				-xAxis_min	0
				-xAxis_max	10000
				-yAxis_min	0
				-yAxis_max	10000
				
				-xAxis_lab_rota		-45 
				-yAxis_lab_rota		0
				
				-x_tickInterval
				
				-point_start
				
				-legendx
				-legendy
				
				-marker true or false
				
				
				-scatter_series		the sample sequence				
				-scatter_size

				-scatter_symbol
				
				
				-pdfcompress no
				-jsdelay 20000
				-h	
                                
USAGE
if(!$opts{t} || !$opts{type} || $opts{h}){
	die $usage;
}

my $gcolor = "'#3366CC','#DC3912','#FF9900','#109618','#990099','#0099C6','#DD4477','#66AA00','#B82E2E','#316395','#994499','#22AA99','#AAAA11','#45AFE2','#FF3300','#FFCC00','#14C21D','#DF51FD','#15CBFF','#FF97D2','#97FB00','#DB6651','#518BC6','#BD6CBD','#35D7C2','#E9E91F'";
my $gcolor_light = "'#45AFE2','#FF3300','#FFCC00','#14C21D','#DF51FD','#15CBFF','#FF97D2','#97FB00','#DB6651','#518BC6','#BD6CBD','#35D7C2','#E9E91F','#3366CC','#DC3912','#FF9900','#109618','#990099','#0099C6','#DD4477','#66AA00','#B82E2E','#316395','#994499','#22AA99','#AAAA11'";
my $colors = {aliceblue=>"#f0f8ff", antiquewhite=>"#faebd7", aqua=>"#00ffff", aquamarine=>"#7fffd4", azure=>"#f0ffff", beige=>"#f5f5dc", bisque=>"#ffe4c4", black=>"#000000", blanchedalmond=>"#ffebcd", blue=>"#0000ff", blueviolet=>"#8a2be2", brown=>"#a52a2a", burlywood=>"#deb887", cadetblue=>"#5f9ea0", chartreuse=>"#7fff00", chocolate=>"#d2691e", 
coral=>"#ff7f50", cornflowerblue=>"#6495ed", cornsilk=>"#fff8dc", crimson=>"#dc143c", cyan=>"#00ffff", darkblue=>"#00008b", darkcyan=>"#008b8b", darkgoldenrod=>"#b8860b", darkgray=>"#a9a9a9", darkgreen=>"#006400", darkgrey=>"#a9a9a9", darkkhaki=>"#bdb76b", darkmagenta=>"#8b008b", darkolivegreen=>"#556b2f", darkorange=>"#ff8c00", 
darkorchid=>"#9932cc", darkred=>"#8b0000", darksalmon=>"#e9967a", darkseagreen=>"#8fbc8f", darkslateblue=>"#483d8b", darkslategray=>"#2f4f4f", darkslategrey=>"#2f4f4f", darkturquoise=>"#00ced1", darkviolet=>"#9400d3", deeppink=>"#ff1493", deepskyblue=>"#00bfff", dimgray=>"#696969", dimgrey=>"#696969", dodgerblue=>"#1e90ff", firebrick=>"#b22222",
 floralwhite=>"#fffaf0", forestgreen=>"#228b22", fuchsia=>"#ff00ff", gainsboro=>"#dcdcdc", ghostwhite=>"#f8f8ff", goldenrod=>"#daa520", gold=>"#ffd700", gray=>"#808080", green=>"#008000", greenyellow=>"#adff2f", grey=>"#808080", honeydew=>"#f0fff0", hotpink=>"#ff69b4", indianred=>"#cd5c5c", indigo=>"#4b0082", ivory=>"#fffff0", khaki=>"#f0e68c",
 lavenderblush=>"#fff0f5", lavender=>"#e6e6fa", lawngreen=>"#7cfc00", lemonchiffon=>"#fffacd", lightblue=>"#add8e6", lightcoral=>"#f08080", lightcyan=>"#e0ffff", lightgoldenrodyellow=>"#fafad2", lightgray=>"#d3d3d3", lightgreen=>"#90ee90", lightgrey=>"#d3d3d3", lightpink=>"#ffb6c1", lightsalmon=>"#ffa07a", lightseagreen=>"#20b2aa",
 lightskyblue=>"#87cefa", lightslategray=>"#778899", lightslategrey=>"#778899", lightsteelblue=>"#b0c4de", lightyellow=>"#ffffe0", lime=>"#00ff00", limegreen=>"#32cd32", linen=>"#faf0e6", magenta=>"#ff00ff", maroon=>"#800000", mediumaquamarine=>"#66cdaa", mediumblue=>"#0000cd", mediumorchid=>"#ba55d3", mediumpurple=>"#9370d8", mediumseagreen=>"#3cb371",
 mediumslateblue=>"#7b68ee", mediumspringgreen=>"#00fa9a", mediumturquoise=>"#48d1cc", mediumvioletred=>"#c71585", midnightblue=>"#191970", mintcream=>"#f5fffa", mistyrose=>"#ffe4e1", moccasin=>"#ffe4b5", navajowhite=>"#ffdead", navy=>"#000080", oldlace=>"#fdf5e6", olive=>"#808000", olivedrab=>"#6b8e23", orange=>"#ffa500", orangered=>"#ff4500",
 orchid=>"#da70d6", palegoldenrod=>"#eee8aa", palegreen=>"#98fb98", paleturquoise=>"#afeeee", palevioletred=>"#d87093", papayawhip=>"#ffefd5", peachpuff=>"#ffdab9", peru=>"#cd853f", pink=>"#ffc0cb", plum=>"#dda0dd", powderblue=>"#b0e0e6", purple=>"#800080", red=>"#ff0000", rosybrown=>"#bc8f8f", royalblue=>"#4169e1", saddlebrown=>"#8b4513",
 salmon=>"#fa8072", sandybrown=>"#f4a460", seagreen=>"#2e8b57", seashell=>"#fff5ee", sienna=>"#a0522d", silver=>"#c0c0c0", skyblue=>"#87ceeb", slateblue=>"#6a5acd", slategray=>"#708090", slategrey=>"#708090", snow=>"#fffafa", springgreen=>"#00ff7f", steelblue=>"#4682b4", tan=>"#d2b48c", teal=>"#008080", thistle=>"#d8bfd8", tomato=>"#ff6347",
 turquoise=>"#40e0d0", violet=>"#ee82ee", wheat=>"#f5deb3", white=>"#ffffff", whitesmoke=>"#f5f5f5", yellow=>"#ffff00", yellowgreen=>"#9acd32"};


my $type = $opts{type};
my $table = $opts{t};
my $point_start = (defined $opts{point_start})?"pointStart: ".$opts{point_start}.",":"";
my $x_tickInterval = (defined $opts{x_tickInterval})?"tickInterval: ".$opts{x_tickInterval}.",":"";

my $base = $table;
#$base =~ 's/\.*//g';

my $js_dir = $RealBin."/js/";

my $title = $opts{title}?$opts{title}:"$base";
$title =~ s/.*\///g;
my $width = $opts{width}?$opts{width}:900;
my $height = $opts{height}?$opts{height}:500;
my $marker = $opts{marker}?$opts{marker}:"true";
my $piemark = (defined $opts{piemark})?$opts{piemark}:"";

my $legendx =  $opts{legendx}?$opts{legendx}:"-100";
my $legendy =  (defined $opts{legendy})?$opts{legendy}:"200";

my $outfile = $base.$type.".html";
my $header = (defined $opts{head})?$opts{head}:1;
my $xAxis_lab = $opts{xAxis_lab}?$opts{xAxis_lab}:"";
my $yAxis_lab = $opts{yAxis_lab}?$opts{yAxis_lab}:"";
# my $yAxis_min =  $opts{yAxis_min}?"min: ".$opts{yAxis_min}.",":"";

my $color_the = $opts{color_the}?$opts{color_the}:"default";

my $alignTicks ="true";

my $jsdelay = $opts{jsdelay}?$opts{jsdelay}:20000;


my @col_name;
my @row_name;

my %datas;
my $xAxis_cate ="";


open IN,"$opts{t}" || die "can not open $opts{t}";
if($header == 1){
	my $head = <IN>;
	chomp($head);
	@col_name = split(/\t/,$head);
	my $typer = shift @col_name;
	if(! $xAxis_lab){
		$xAxis_lab = $typer;
	}
	if(! $yAxis_lab){
		$yAxis_lab = join(" ",@col_name);
	}
}

while(<IN>){
	chomp;
	my @line =  split(/\t/,$_);
	push(@row_name,"'".$line[0]."'");
	
	if($header == 0){
		push(@col_name,(1..$#line));
		$header =1;
	}
	
	foreach(my $i=1;$i<=$#line;$i++){
		#print $line[0]."\n".$col_name[$i-1]."\n".$line[$i]."\n";
		$datas{"'".$line[0]."'"}{$col_name[$i-1]} = $line[$i];
	}
}


$xAxis_cate = "categories: [".join(",",@row_name)."],";
if($type eq "scatter"){
	$xAxis_lab = $col_name[0];
	$xAxis_cate = "";
	$yAxis_lab = $col_name[1];
	$alignTicks = "false";
}

close IN;
my $s_num = $#col_name;
if($type eq "pie"){
	$s_num = $#row_name;
}
#print $s_num;
$legendy = 200 - 10*$s_num;
if ($legendy <20){$legendy = 20;}

my $set_col_scheme = "";
if($color_the eq "default"){

}elsif($color_the eq "google"){
$set_col_scheme = "Highcharts.setOptions({
				colors: [$gcolor]
			});";
}elsif($color_the eq "googlelight"){
$set_col_scheme = "Highcharts.setOptions({
				colors: [$gcolor_light]
			});";
}else{
$set_col_scheme = "Highcharts.setOptions({
				colors: [$color_the]
			});";
}

my $legend_ht = "";
if($type eq "pie" ||$type eq "area" || $type eq "scatter" ){
	$legend_ht = "
		legend: {
			align: 'right',
			verticalAlign: 'center',
			layout: 'vertical',
			symbolHeight: 16,
			itemMarginBottom: 4,
			
			x: $legendx,
			y: $legendy,
		},";
}elsif($type eq "boxplot" || $opts{nolegend}){
	$legend_ht = "
		legend: {
			enabled: false
		},";
}

my $chart_ht = "
		chart: {
			type: '$type',
			alignTicks: $alignTicks
		},";
my $tytle_ht = "
		title: {
			text: '$title',
			style: {
				color: '#000000',
				fontWeight: 'bold',
				fontSize: 20,
				fontFamily: 'Times New Roman'
			}
		},";
my $tooltip_ht = "
";

my $grid_y =  "";
my $line_other ="";

if($type eq "scatter" || $type eq "area" ){
	if($type eq "scatter"){
		$xAxis_lab = $col_name[0];
		$yAxis_lab = $col_name[1];
	}
	$alignTicks = "false";
	$grid_y = "gridLineWidth: 0,";
	$line_other = "lineColor: '#000000',
			lineWidth: 1,";
}

my $x_other = $line_other;
my $y_other = $line_other;
if(defined $opts{xAxis_min}){
	$x_other .= "
			min: $opts{xAxis_min},
	";
}
if($opts{xAxis_max}){
	$x_other .= "
			max: $opts{xAxis_max},
	";
}

#print "\n****".$opts{yAxis_min}."****";
if(defined $opts{yAxis_min}){
	$y_other .= "
			min: $opts{yAxis_min},
	";
}
if(defined $opts{yAxis_max}){
	$y_other .= "
			max: $opts{yAxis_max},
	";
}

if(defined $opts{yAxis_reversed}){
	$y_other .= "
			reversed: true,
	";
}
			# type: 'logarithmic',

my $xAxis_log_ht ="";
my $yAxis_log_ht ="";

if ($opts{xAxis_log}){
	$xAxis_log_ht = "type: 'logarithmic',";	
}
if ($opts{yAxis_log}){
	$yAxis_log_ht = "type: 'logarithmic',";	
}


my $xAxis_ht = "
		xAxis: {
			$x_tickInterval
			$xAxis_cate
			$x_other
			$xAxis_log_ht
			title: {
				text: '$xAxis_lab',
				style: {
					color: '#000000',
				}
			}
		},";
my $yAxis_ht = "
		yAxis: {
			labels: {
				format: '{value}'
			},
			$grid_y
			$y_other
			$yAxis_log_ht
			title: {
				text: '$yAxis_lab',
				style: {
					color: '#000000',
				}
			}

		},";


if($type eq "pie"){
	$xAxis_ht = "";
	$yAxis_ht = "";
}
		
my $plotOptions_ht = &get_option_ht($type);

my $series_ht = &get_series_ht($type);

my $html = <<"HTML";
<html>
	<head>
		<meta http-equiv="$table.$type" content="text/html; charset=utf-8">
		<title>Highcharts Example</title>		
		<script type="text/javascript" src="${js_dir}jquery.min.js"></script>
		<style type="text/css">
		\$\{demo.css\}
		</style>
		<script type="text/javascript">
		\$(function () {
		$set_col_scheme			
		\$('\#container').highcharts(\{
$chart_ht
$tytle_ht
$tooltip_ht        
$plotOptions_ht 
$legend_ht
$xAxis_ht
$yAxis_ht
$series_ht
		\});
		});
		</script>
	</head>
	<body>

	<script src="$js_dir/highcharts.js"></script>
	<script src="$js_dir/highcharts-3d.js"></script>
	<script src="$js_dir/highcharts-more.js"></script>
	<div id="container" style="height: ${height}px; margin: auto; width: ${width}px"></div>
	</body>

</html>                                
HTML

open OUT,"> $opts{t}.html" || die "can not open $opts{t}.html";
print OUT $html;
close OUT;

my $pdf_width = int($width/3);
my $pdf_height = int($height/3);

my $draw = "wkhtmltopdf --no-stop-slow-scripts --page-width $pdf_width --page-height $pdf_height";
if(defined $opts{pdfcompress} && $opts{pdfcompress} eq "no"){
	$draw .= " --no-pdf-compression";
}
if($jsdelay){
	$draw .= " --javascript-delay $jsdelay";
}
$draw .=  " $opts{t}.html $opts{t}.pdf";

#print $draw."\n";
system($draw);

sub get_series_ht(){
	my $types =shift; 
	my $data= "";
	my @series;
	if($types eq "column" || $types eq "area" || $types eq "areaspline" || $types eq "bar" || $types eq "line"){
		foreach(@col_name){
			my $col = $_;
			my @col_data;
			foreach(@row_name){
				#print $_."\n".$col."\n".$datas{$_}{$col};
				push(@col_data, $datas{$_}{$col});
			}
			$data = "{
				name: '$col',
				data: [".join(",",@col_data)."],
				marker:{
					enabled: ".$marker."
				}
			}";
			push(@series,$data);
		}
	}elsif($types eq "pie" ){
		if($#col_name >0){
			die "$types can only have 1 cols";
		}
		my @col_data;
		foreach(@row_name){
			if($piemark ne "" &&  "'".$piemark."'" eq $_){
				my $mark_pie = "{
					name: $_,
					y: $datas{$_}{$col_name[0]},
					sliced: true,
					selected: true
				}";
				push (@col_data, $mark_pie);
			}else{
				push (@col_data, "[$_, $datas{$_}{$col_name[0]}]");
			}
		}
		$data = "{
				type: 'pie',
				name: '$col_name[0]',
				data: [
				".join(",\n\t\t\t\t",@col_data)."
				]
		}";
		push(@series,$data);			
	}elsif($types eq "scatter"){
		if($#col_name != 2){
			die "$types can only have 3 cols";
		}
		my @series;		

		my @scatter_series = split(/,/, $opts{scatter_series});
		
		my @scatter_size = split(/,/, $opts{scatter_size});
		my @scatter_symbol = split(/,/, $opts{scatter_symbol});
		# print STDERR @scatter_size;
		# print STDERR @scatter_symbol;
		my $i=0;
		#print STDERR @scatter_series;
		
		foreach(@scatter_series){
			
			my $sc_type = $_;
			my @row_data;
			foreach(@row_name){				
				
				if($datas{$_}{$col_name[2]} eq $sc_type){
					if( $opts{xAxis_min} && $datas{$_}{$col_name[0]} < $opts{xAxis_min} ){
						$datas{$_}{$col_name[0]} = $opts{xAxis_min};
					}
					if( $opts{yAxis_min} && $datas{$_}{$col_name[1]} < $opts{yAxis_min} ){
						$datas{$_}{$col_name[1]} = $opts{yAxis_min};
					}
					push (@row_data, "[$datas{$_}{$col_name[0]}, $datas{$_}{$col_name[1]}]");
				}
			}
			$data = "{
				name: '".$sc_type."',
				data:[
					".join(",",@row_data)."
				],
				marker:{
					radius: ".$scatter_size[$i].",
					symbol: ".$scatter_symbol[$i].",
				}
				
			}";	
			push(@series,$data);
			$i ++;
		}
		my $ser_result = "\t\tseries: [".join(",",@series)."]";	
		return $ser_result;		
		
	
	}elsif($types eq "arearange"){
		if($#col_name != 2){
			die "$types can only have 3 cols";
		}
		
		my @med_data;
		my @range_data;
		foreach(@row_name){
			push (@med_data, "['$_', $datas{$_}{$col_name[0]}]");
			push (@range_data, "['$_', $datas{$_}{$col_name[1]}, $datas{$_}{$col_name[2]}]");
		}
		
		$data = "
				{
				name: '$col_name[0]',
				data: [
				".join(",\n\t\t\t\t",@med_data)."
				];
				},{
				name: '$col_name[1]-$col_name[2]',
				data: [
				".join(",\n\t\t\t\t",@range_data)."
				];
				}";				

		push(@series,$data);
		my $ser_result = "\t\tseries: [".join(",",@series)."]";	
		return $ser_result;
		
	}elsif($types eq "boxplot"){
		my @all_data;
		foreach(@row_name){
			my $row = $_;
			my @row_data;
			foreach(@col_name){
				push(@row_data, $datas{$row}{$_});
			}
			push (@all_data, "[".join(",",@row_data)."]");		
		}	
		$data = "{
				name: '".join("_",@col_name)."',
				data:[
					".join(",\n\t\t\t\t\t",@all_data)."
				],
			}";			
		push(@series,$data);
	}else{
		die "can not draw $types";
	}
	my $ser_result = "\t\tseries: [".join(",",@series)."]";	
	return $ser_result;
}

sub get_option_ht(){
	my $types =shift;
	my $result = "";
	if($types eq "pie"){
		$result = "
		plotOptions: {
			pie: {
				allowPointSelect: true,
				cursor: 'pointer',
				dataLabels: {
					enabled: true,
					format: \"{y}\"
				},
				showInLegend: true
			}
		},";
	}elsif($types eq "area"){
		$result = "
		plotOptions: {
			series: {
                lineWidth: 2,
				fillOpacity: 0.5
            },
			".$type.": {
				".$point_start."
				marker: {
					enabled: false,
					symbol: 'circle',
					radius: 2,
					states: {
						hover: {
							enabled: true
						}
					}
				}
			}
		},
		"
	}else{
	}
	return $result;
}

