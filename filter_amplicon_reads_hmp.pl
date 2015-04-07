#!/usr/bin/perl -w
#filename:filter_amplicon_reads_hmp.pl

use strict;
use warnings;
use Getopt::Long;

############### FILTERING AIMS ###############

#based on paper: Consortium, The Human Microbiome Project (2012), 'A framework for human microbiome research', Nature, 486 (7402), 215-21

#using the high stringency pipeline filtering

#remove: sequences with an ambiguous base call
#remove: sequences with homopolymer longer than 8nt
#remove: short sequences 

#calc average quality score within 50bp window moving along sequence, when ave quality drops below 35 trim sequence


#remove chimeras - ChimeraSlayer trained to Gold database aligned to the SILVA reference alignment


### other filtering ###
#remove: sequences longer than amplicon


### low stringency ###
#trim at position where cumulative average quality drops below 35





############### INPUTS ####################
my $input_format = "$0 <read_file> <chimera_slayer.CPC file> Options: -filter/-f <high>\n";

my $read_file = shift @ARGV or die "$input_format\n";
my $chimera_slayer_file = shift @ARGV or die "$input_format\n";
my $high_stringency = 0;

GetOptions
    ('f|filter=s'      => \$high_stringency);

my $total_reads = 0;
my $chimera_count = 0;
my $passed_reads = 0;
my $failed_reads = 0;

my $output_pass = "pass_filter_$read_file";
my $output_fail = "failed_filter_$read_file";
my $output_chimera = "failed_chimera_filter_$read_file";

############### PUT READS IN HASH TABLE ####################
my %read_hash;

use Bio::SeqIO;
my $inseq = Bio::SeqIO->new('-file' => "<$read_file",'-format' => 'fastq');
while (my $seq_obj = $inseq->next_seq)
{   
    $total_reads ++;
    
    my $id = $seq_obj->id;	
    $read_hash{$id} = $seq_obj;
}

############### CHIMERA SLAYER ####################

open (SLAYER, "$chimera_slayer_file") or die "missing output file from ChimeraSlayer\n";
while (<SLAYER>)
{
    my $file_line = $_;
    my @split_file_line = split ("\t",$file_line);

# parameters following mother chimera slayer default parameters (the implementation used by HMP) 
if (($split_file_line[4]>=1.007 && $split_file_line[5]>=90 && $split_file_line[6]>=90) || ($split_file_line[7]>=1.007 && $split_file_line[8]>=90 && $split_file_line[9]>=90))
    {
	my $id = $split_file_line[1];
	my $seq_obj = $read_hash{$id};
	if ($seq_obj) 
	{
	    $chimera_count ++;
	    $failed_reads ++;
	    
	    my $read = $seq_obj->seq;
	    chomp $read;
	    my $qual_ref = $seq_obj->qual;
	    my $qual = join ("", map {chr($_+33)} @{$qual_ref});
	    output($id,$read,$qual,$output_chimera);
	    
	    delete ($read_hash{$split_file_line[1]});
	}
    }
    
}

close SLAYER; 

############### PROCESS READ FILE ####################

foreach my $id (keys %read_hash)
{
    my $seq_obj = $read_hash{$id};
    my $read = $seq_obj->seq;
    chomp $read;
    my $qual_ref = $seq_obj->qual; 
    
    
    if (filter_long_seq ($read) and filter_short_seq ($read))#filter long reads (likely missed chimeras) before qual trim shortens them #remove short reads without further processing
    {
	my $qual;
	if ($high_stringency)
	{
	($read, $qual) = trim_by_quality_high_stringency ($read,$qual_ref);
	}
	else 
	{
	($read, $qual) = trim_by_quality_low_stringency ($read,$qual_ref);
	}	
	
	chomp $qual;	

	#remove reads qual trimmed below threshold 
	if (filter_short_seq ($read) and 
	    filter_N ($read) and
	    filter_homopolymer ($read))
	{
	    $passed_reads ++;
	    output($id,$read,$qual,$output_pass);
	}
	else
	{
	    $failed_reads ++;
	    output($id,$read,$qual,$output_fail); 
	}		
    }
    else
    {
	my $qual= join ("", map {chr($_+33)} @{$qual_ref});
	$failed_reads ++;
	output($id,$read,$qual,$output_fail);
    }
}
############### LOG FILE ####################
my $percent_passed = calc_percent($total_reads,$passed_reads);
my $percent_failed = calc_percent ($total_reads,$failed_reads);
my $percent_chimera = calc_percent ($total_reads,$chimera_count);

open OUTFILE, ">log_filter_$read_file";
if ($high_stringency)
{
print OUTFILE "Filtered with HIGH stringency options\n\n";
}
else
{
print OUTFILE "Filtered with LOW stringency options\n\n";
}
print OUTFILE "Total Reads\t$total_reads\nTotal Passed\t$passed_reads\t$percent_passed%\nTotal Failed\t$failed_reads\t$percent_failed%\nRemoved Chimeras\t$chimera_count\t$percent_chimera%";
close OUTFILE;




############### SUBROUTINES ##########################################################

sub trim_by_quality_high_stringency #calc average quality score within 50bp window moving along sequence, when ave quality drops below 35 trim sequence
{
    use List::Util "sum";
    my ($read,$read_qual_ref) = @_;
    
    my @qual_array = @{$read_qual_ref};

    my $current_window = (sum(@qual_array[0 .. 49]))/50;
    my $count = 0;

    while (($current_window >= 35) && ((50+$count) < (length $read)))
    {
	$count ++;
	$current_window = (sum(@qual_array[(0+$count) .. (49+$count)]))/50;
    }

    my $trimmed_read_quality;
    if ($current_window >= 35)
    {
	$read = substr ($read, 0, (50+$count));
	$trimmed_read_quality = join ("", map {chr($_+33)} @{$read_qual_ref}[0..((length $read)-1)]); 
    }
    else
    {
	$read = substr ($read, 0, (49+$count));
	$trimmed_read_quality = join ("", map {chr($_+33)} @{$read_qual_ref}[0..((length $read)-1)]); 
    }
    
    return ($read, $trimmed_read_quality);
}

sub trim_by_quality_low_stringency #trim at position where cumulative average quality drops below 35
{
use List::Util "sum";
    my ($read,$read_qual_ref) = @_;
    
    my @qual_array = @{$read_qual_ref};
my $count = 0;

my $cum_ave = $qual_array[0];

while (($cum_ave >= 35) && (($count+1) < (length $read)))
    {
	$count ++;
	$cum_ave = (sum(@qual_array[0 .. $count]))/($count+1);
    }

my $trimmed_read_quality;
    if ($cum_ave >= 35)
    {
	$read = substr ($read, 0, ($count+1));
	$trimmed_read_quality = join ("", map {chr($_+33)} @{$read_qual_ref}[0..((length $read)-1)]); 
    }
    else
    {
	$read = substr ($read, 0, ($count));
	$trimmed_read_quality = join ("", map {chr($_+33)} @{$read_qual_ref}[0..((length $read)-1)]); 
    }
    
    return ($read, $trimmed_read_quality);
}

sub filter_long_seq 
{
    my $read = $_[0];
    my $pass = 0;
    if ((length $read) <= 579)
    {
	$pass = 1;
    }
    return $pass;
}

sub filter_short_seq
{
    my $read = $_[0];
    my $pass = 0;
    if ((length $read) >= 200)
    {
	$pass = 1;
    }
    return $pass;
}

sub filter_N 
{
    my $read = $_[0];
    my $pass = 0;
    if ($read !~ m/[^AGCT]+/)
    {
	$pass = 1;
    }
    return $pass;
}

sub filter_homopolymer 
{
    my $read = $_[0];
    my $pass = 0;
    if ($read !~ m/.*A{9,}.*|.*G{9,}.*|.*C{9,}.*|.*T{9,}.*/)
    {
	$pass = 1;
    }
    return $pass;
}

sub output
{
    my ($id,$read,$trimmed_read_quality,$output_file) = @_;
    open OUTFILE, ">>$output_file";
    print OUTFILE "\@$id\n$read\n+\n$trimmed_read_quality\n";
    close OUTFILE;
}

sub calc_percent
{
    my ($total,$amount) = @_;
    my $percent = sprintf("%.1f",(($amount/$total)*100)); 
    return $percent;
}

