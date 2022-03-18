#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts('ECMSA:', \%options);

#####
#
#-E = Print alignment error scores for each tag
#-C = Print CIGAR strings for each tag
#-M = Print MD tag for each tag
#-S = Print Start/Stop for oligo alignment
#-A = Alignment cutoff to use for barcode assignment (default = 0.05 = 5% error)
#####

my $ERR_flag;
if(exists($options{E}))
	{
	print STDERR "Appending alignment scores\n";
	$ERR_flag = 1;
	}
else {$ERR_flag = 0;}

my $CIGAR_flag;
if(exists($options{C}))
	{
	print STDERR "Appending CIGAR strings\n";
	$CIGAR_flag = 1;
	}
else {$CIGAR_flag = 0;}

my $MD_flag;
if(exists($options{M}))
	{
	print STDERR "Appending cs Tags\n";
	$MD_flag = 1;
	}
else {$MD_flag = 0;}

my $POS_flag;
if(exists($options{S}))
	{
	print STDERR "Appending Alignment start/stop\n";
	$POS_flag = 1;
	}
else {$POS_flag = 0;}


my $aln_cutoff = $options{A} || 0.05;
print STDERR "Using $aln_cutoff error rate for alignment cutoff\n";


my $list = $ARGV[0]; #File with ID and filename
my $out = $ARGV[1];

#open (MATCH, "$match_file") or die("ERROR: can not read file ($match_file): $!\n");

my @inline;
my %file_list;
my @ordered_list;

open (LIST, "$list") or die("ERROR: can not read file ($list): $!\n");
	while (<LIST>)
		{
		chomp;
		@inline = split(/\s+/);
		$file_list{$inline[0]}=$inline[1];
		push(@ordered_list,$inline[0]);
		}

my %sample;
my %sample_stats;
my %sample_stats_B;
my %counts;
my %cigar;
my %md;
my %aln;
my %pos;

my %oligo_id;
my $sample_ID;
my $key;

my $barcode;
my $bc_ct;
my $flag_B;
my $oligo;
my $bc_flag;
my $bc_aln;
my $bc_cigar;
my $bc_md;
my $bc_pos;

my $cur_file;

foreach $sample_ID (@ordered_list)
	{
	$cur_file=$file_list{$sample_ID};

	print STDERR "Reading $sample_ID\n";
	print "Stats for $sample_ID\n";

	open (COUNTS, "$cur_file") or die("ERROR: can not read file ($cur_file): $!\n");
	while (<COUNTS>)
		{

		chomp;
		@inline = split("\t");
		$barcode=$inline[0];
		$bc_ct=$inline[1];
		$flag_B=$inline[2];
		$oligo=$inline[3];
		$bc_flag=$inline[4];
		$bc_aln=$inline[6];
		$bc_cigar=$inline[7];
		$bc_md=$inline[8];
		$bc_pos=$inline[9];

		$sample_stats{$sample_ID}{$bc_flag}{"ct"}++;
		$sample_stats{$sample_ID}{$bc_flag}{"sum"}+=$bc_ct;

		$sample_stats_B{$sample_ID}{$flag_B}{"ct"}++;
		$sample_stats_B{$sample_ID}{$flag_B}{"sum"}+=$bc_ct;

		if(($bc_flag eq 0 || $bc_flag eq 2) && $oligo ne "*")
			{
			if($aln_cutoff >= $bc_aln)
				{
				die "Barcode & Sample combination seen twice\n" if(exists $counts{$barcode}{$sample_ID});
				$counts{$barcode}{$sample_ID}=$bc_ct;
				$sample_stats{$sample_ID}{"counted"}{"ct"}++;
				$sample_stats{$sample_ID}{"counted"}{"sum"}+=$bc_ct;
				if(exists $oligo_id{$barcode})
					{
					die "Barcodes seen with different oligo IDs\n$barcode\n$oligo_id{$barcode}\n$cur_file\n$oligo\n" if($oligo_id{$barcode} ne $oligo);
					die "Barcodes seen with different flag IDs\n$barcode\n$aln{$barcode}\n$cur_file\n$bc_aln\n" if($aln{$barcode} ne $bc_aln);
					die "Barcodes seen with different cigar IDs\n$barcode\n$cigar{$barcode}\n$cur_file\n$bc_cigar\n" if($cigar{$barcode} ne $bc_cigar);
					die "Barcodes seen with different cs tag IDs\n$barcode\n$md{$barcode}\n$cur_file\n$bc_md\n" if($md{$barcode} ne $bc_md);
					die "Barcodes seen with different start/stop positions IDs\n$barcode\n$pos{$barcode}\n$cur_file\n$bc_pos\n" if($pos{$barcode} ne $bc_pos);

					}
				else
					{
					$oligo_id{$barcode}=$oligo;
					$aln{$barcode}=$bc_aln;
					$cigar{$barcode}=$bc_cigar;
					$md{$barcode}=$bc_md;
					$pos{$barcode}=$bc_pos;
					}
				}
			elsif($aln_cutoff == 0.05 && $bc_flag eq 0)
				{
				die "$barcode alignment score is under pipeline cutoff of 5% but flagged??\n$cur_file\n";
				}
			elsif($aln_cutoff == 0.05 && $bc_flag eq 2 && $bc_aln <= 0.05)
				{
				die "$barcode alignment score is under pipeline cutoff of 5% but flagged??\n$cur_file\n";
				}
			}
		}

	print "Flag\tBC Count\tRead Sum\n";
	foreach my $key (sort { $a cmp $b} keys %{$sample_stats{$sample_ID}})
		{
    	print join("\t",$key,$sample_stats{$sample_ID}{$key}{"ct"},$sample_stats{$sample_ID}{$key}{"sum"})."\n";
		}
	#print "Flag B\tBC Count\tRead Sum\n";
	#foreach my $key (sort { $a cmp $b} keys %{$sample_stats_B{$sample_ID}})
	#	{
    #	print join("\t",$key,$sample_stats_B{$sample_ID}{$key}{"ct"},$sample_stats_B{$sample_ID}{$key}{"sum"})."\n";
	#	}
	close COUNTS;
	}

print STDERR "Writing file\n";

open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");
print OUT "Barcode\tOligo\t";
print OUT "Error\t" if($ERR_flag == 1);
print OUT "CIGAR\t" if($CIGAR_flag == 1);
print OUT "cs\t" if($MD_flag == 1);
print OUT "Aln_Start:Stop\t" if($POS_flag == 1);
print OUT join ("\t",@ordered_list)."\n";

my $cur_bc;
my $cur_sample;

foreach $cur_bc (keys %counts)
	{
	print OUT "$cur_bc\t$oligo_id{$cur_bc}";
	print OUT "\t$aln{$cur_bc}" if($ERR_flag == 1);
	print OUT "\t$cigar{$cur_bc}" if($CIGAR_flag == 1);
	print OUT "\t$md{$cur_bc}" if($MD_flag == 1);
	print OUT "\t$pos{$cur_bc}" if($POS_flag == 1);

	foreach $cur_sample (@ordered_list)
		{
		if(exists $counts{$cur_bc}{$cur_sample})
			{
			print OUT "\t$counts{$cur_bc}{$cur_sample}";
			}
		else
			{
			print OUT "\t0";
			}
		}
		print OUT "\n"

	}
close OUT
