#!/usr/bin/perl -w
use warnings;
use strict;

use Cwd;
use Getopt::Long;
use Data::Dumper;

my %opts;
GetOptions (\%opts,"i=s","o=s","g=s","u=s","d=s","h!");

my $usage = <<"USAGE";
		Program:	$0
		Contact:	ting.kuang\@majorbio.com
        Discription:	This script is used for iTRAQ/TMT different expression anlysis 
        Usage:	perl $0 [options]
			-i*	expression matrix of iTRAQ
				############################
				#Accession	DMSO1	DMSO2	DM1	DM2	BPA1	BPA2
				#B9EKP5	1	0.966	1.036	0.998	0.987	0.977
				############################
			-o	output dir, default: diffexp
			-g	group info, e.g. DM/DMSO, BPA/DMSO, BPA/DM (if no this parameter, all pairwise comparisons: DMSO2/DMSO1, DM1/DMSO1, DM2/DMSO1, ...)
				############################
				#DMSO	DMSO1
				#DMSO	DMSO2
				#DM		DM1
				#DM		DM2
				#BPA	BPA1
				#BPA	BPA2
				############################
			-u	up-regulate, default: 1.2
			-h	show this information
        Example:	perl $0 -i MJ_expression.xls -g group.config -u 1.2
USAGE

die $usage if ( !$opts{i} || $opts{h} );

$opts{u}=defined$opts{u}?$opts{u}:1.20;
$opts{o}=defined$opts{o}?$opts{o}:"diffexp";
my $down = 1/$opts{u};

my @sample_names;
my @group_names;

my %sample_name_to_column = &get_sample_name_to_column_index($opts{i});

my %samples;
if ($opts{g}) {
	unless ($opts{g} =~ /^\//) {
		$opts{g} = cwd() . "/$opts{g}";
	}
	%samples = &parse_sample_info($opts{g});
} else {
	# no replicates, so assign each sample to itself as a single replicate
	foreach my $sample_name (keys %sample_name_to_column) {
		$samples{$sample_name} = [$sample_name];
	}
}

print Dumper(\%samples);
        
if ($opts{i} !~ /^\//) {
	## make full path
	$opts{i} = cwd() . "/$opts{i}";
}
    
mkdir $opts{o} or die "Error, cannot mkdir $opts{o}";
chdir $opts{o} or die "Error, cannot cd to $opts{o}";
    
if($opts{g}){
	@sample_names = @group_names;
}

my @A_samples = @sample_names;
my @B_samples = @sample_names;	


# all pairwise comparisons
for (my $i = 0; $i < $#A_samples; $i++) {
	my $sample_i = $A_samples[$i];
	for (my $j = $i + 1; $j <= $#B_samples; $j++) {
		my $sample_j = $B_samples[$j];
		my ($sample_a, $sample_b) = ($sample_i, $sample_j);
		
		&run_ttest_sample_pair($opts{i}, \%samples, \%sample_name_to_column, $sample_a, $sample_b, $opts{u}, $down);
		
	}
}



#####################################################################
sub parse_sample_info {
    my ($sample_file) = @_;

    my %samples;

    open SF, "<$sample_file" || die $!;
    while (<SF>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } # allow comments
        chomp;
        s/^\s+//; # trim any leading ws
        my @x = split(/\s+/); # now ws instead of just tabs
        if (scalar @x < 2) { next; }
        my ($sample_name, $replicate_name, @rest) = @x;
        
        #$sample_name =~ s/^\s|\s+$//g;
        #$replicate_name =~ s/^\s|\s+$//g;
        if (exists $samples{$sample_name}){

		}else{
			push(@group_names, $sample_name);
		}
		
        push (@{$samples{$sample_name}}, $replicate_name);

    }
    close SF;

    return(%samples);
}

#####################################################################
sub get_sample_name_to_column_index {
    my ($matrix_file) = @_;

    my %column_index;

    open MF, "<$matrix_file" || die $!;
    my $header_line = <MF>;
	chomp($header_line);
    $header_line =~ s/^\#//; # remove comment field.
    #$header_line =~ s/^\s+|\s+$//g;
    my @samples = split(/\t/, $header_line);
	shift @samples;
	@sample_names = @samples;


    { # check for disconnect between header line and data lines
        my $next_line = <MF>;
        my @x = split(/\t/, $next_line);
        print STDERR "Got " . scalar(@samples) . " samples, and got: " . scalar(@x) . " data fields.\n";
        print STDERR "Header: $header_line\nNext: $next_line\n";
        
        if (scalar(@x) == scalar(@samples)) {
            # problem... shift headers over, no need for gene column heading
            shift @samples;
            print STDERR "-shifting sample indices over.\n";
        }
    }
    close MF;
            
    
    my $counter = 0;
    foreach my $sample (@samples) {
        $counter++;
        
        $sample =~ s/\.(isoforms|genes)\.results$//; 
        
        $column_index{$sample} = $counter;
    }

    use Data::Dumper;
    print STDERR Dumper(\%column_index);
    

    return(%column_index);
    
}

####################################################################################

sub run_ttest_sample_pair {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B, $up, $down) = @_;
         
    my $output_prefix = join("_vs_", ($sample_A, $sample_B));
        
    my $Rscript_name = "$output_prefix.t-test.Rscript";
    # print $samples_href."****\n";
    my @reps_A = @{$samples_href->{$sample_A}};
    my @reps_B = @{$samples_href->{$sample_B}};

    my $num_rep_A = scalar(@reps_A);
    my $num_rep_B = scalar(@reps_B);
    
    my @rep_column_indices_a;
    foreach my $rep_name (@reps_A) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices_a, $column_index);
    }
	
	my @rep_column_indices_b;
	foreach my $rep_name (@reps_B) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices_b, $column_index);
    }
        

    ## write R-script to run t test
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
	print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, sep=\"\\t\", com='')\n";
	print $ofh "col_ordering_a = c(" . join(",", @rep_column_indices_a) . ")\n";
	print $ofh "col_ordering_b = c(" . join(",", @rep_column_indices_b) . ")\n";
	
	# no replicate sample
	if ($num_rep_A > 1 && $num_rep_B > 1) {
		print $ofh "
		
outfile <- matrix(nrow=nrow(data),ncol=8)
for(i in 1:nrow(data)){
	acc <- as.character(rownames(data[i,]))
	exp_a <- as.numeric(data[i,col_ordering_a])
	exp_b <- as.numeric(data[i,col_ordering_b])
	avg_a <- mean(exp_a, na.rm=T)
	avg_b <- mean(exp_b, na.rm=T)
	fc <- avg_b/avg_a
	logfc <- log(fc,2)
	
	if (length(col_ordering_a) == 2) {
		if(exp_a[1]==exp_a[2]){
			exp_a[2] = exp_a[2]+(1e-10)
		}
	}
	
	if (length(col_ordering_b) == 2) {
		if(exp_b[1]==exp_b[2]){
			exp_b[2] = exp_b[2]+(1e-10)
		}
	}
	
	t <- t.test(exp_a, exp_b, var.equal = T)
	p <- t\$p.value
	
	if (p <= 0.05) {
		sig <- \"yes\"
	} else {
		sig <- \"no\"
	}

	if (fc >= $up) {
		reg <- \"up\"
	} else if (fc <= $down) {
		reg <- \"down\"
	} else {
		reg <- \"no change\"
	}
	
	outfile[i,] <- c(acc,avg_a,avg_b,fc,logfc,p,sig,reg)

}

colnames(outfile)=c(\"Accession\", \"$sample_A\", \"$sample_B\", \"FC($sample_B/$sample_A)\", \"log2FC($sample_B/$sample_A)\",\"Pvalue($sample_B/$sample_A)\",\"significant\",\"regulate\")
write.table(outfile,file=\"$output_prefix.diff.exp.xls\",col.names=T,row.names=F,quote=F,sep=\"\\t\")

	";
	} else {
		print $ofh "
		
outfile <- matrix(nrow=nrow(data),ncol=6)
for(i in 1:nrow(data)){
	acc <- as.character(rownames(data[i,]))
	exp_a <- as.numeric(data[i,col_ordering_a])
	exp_b <- as.numeric(data[i,col_ordering_b])
	fc <- exp_b/exp_a
	logfc <- log(fc,2)
	
	if (fc >= $up) {
		reg <- \"up\"
	} else if (fc <= $down) {
		reg <- \"down\"
	} else {
		reg <- \"no change\"
	}
	
	outfile[i,] <- c(acc,exp_a,exp_b,fc,logfc,reg)
}
colnames(outfile)=c(\"Accession\", \"$sample_A\", \"$sample_B\", \"FC($sample_B/$sample_A)\", \"log2FC($sample_B/$sample_A)\",\"regulate\")
write.table(outfile,file=\"$output_prefix.diff.exp.xls\",col.names=T,row.names=F,quote=F,sep=\"\\t\")
";
	}	

    close $ofh;

    ## Run R-script
    my $cmd = "Rscript $Rscript_name";


    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "$@\n\n";
        print STDERR "\n\nWARNING: This t test comparison failed...\n\n";
        ## if this is due to data paucity, such as in small sample data sets, then ignore for now.
    }
    

    return;
}

#################################
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret ($ret) ";
    }

    return;
}

