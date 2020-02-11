#!/usr/bin/perl
#
#by Ryan Tewhey
use strict;
use warnings;
use Getopt::Std;

##################################
#
#4/27/2014 - Changed "RC" reporter to after SNP ID
##################################

my %options=();
getopts('CB', \%options);

#####
#
#-C = Use new cigar format =/M
#-B = Filter for 0 bitflag
#####

my $CIGAR_flag;

if(exists($options{C}))
	{
	print STDERR "Using update CIGAR format =/M\n";
	$CIGAR_flag = 1;
	}
else {$CIGAR_flag = 0;}

my $BIT_flag;

if(exists($options{B}))
	{
	print STDERR "Using only fwd mapping reads (BITFLAG=0)\n";
	$BIT_flag = 1;
	}
else {$BIT_flag = 0;}

my $sam = $ARGV[0];
my $out = $ARGV[1];

open (FASTA, "$sam") or die("ERROR: can not read file ($sam): $!\n");
open (MATCH, ">$out") or die("ERROR: can not create $out .matched: $!\n");

my $score_cutoff = 0.05;

my %chr_size;
my $cigar;
my $cigar_unchanged;
my $chr;
my $pos;
my $start_ins;
my $read_pos;
my $seq;
my @bitflag;
my $match;
my $mismatch;
my $size;
my $cigar_part;
my @id;
my $seq_correct_ori;
my $updated_chr;
my @split_ID;
my $score;
my $parse_col;
my $md_col;
my $md_col_part;
my $md_col_unchanged;
my $mismatch_all;
my $score_all;
my $cigar_substitution;
my $aln_length;
my $unaln_length;
my $aln_info;

while (<FASTA>)
	{
	$_ =~ s/[\n\r]//g;
	my @line = split(/\t/);

	if($line[0] =~ /^@/)
		{
		$line[1] =~ s/SN://;
		$line[2] =~ s/LN://;
		$chr_size{$line[1]} = $line[2];
		}
	else
		{
		 $cigar = $line[5];
		 $cigar_unchanged = $line[5];
		 $chr = $line[2];
		 $pos = $line[3]-1;
		 $start_ins = 0;
		 $read_pos = 0;
		 $seq = $line[9];

		 @id = split(/#/,$line[0]);

		 $size = 0;
		 $size = $chr_size{$chr} if(exists $chr_size{$chr});

		$md_col = "*";
		###Parse MD for mismatches
		foreach $parse_col (@line)
			{
			$md_col = $parse_col if($parse_col =~ m/^MD:Z:/);
			$md_col =~ s/^MD:Z:// if($parse_col =~ m/^MD:Z:/);
			}

		 #ParseFlag $bitflag[7] eq strand (1=neg, 0=pos)
		 @bitflag = split(//,sprintf('%012b',$line[1]));

		 $seq_correct_ori = $seq;
		 $seq_correct_ori = reverse($seq) if($bitflag[7] == 1);
		 $seq_correct_ori =~ tr/ACGTNacgtn/TGCANtgcan/if($bitflag[7] == 1);

		 $match = 0;
		 $mismatch = 0;
		 $cigar_substitution=0;
		 $aln_length=0;

		 while ($cigar !~ /^$/)
			 {
			 if ($cigar =~ /^([0-9]+[MIDSH=X])/)
				 {
				 $cigar_part = $1;
				 if ($cigar_part =~ /(\d+)M/)
					{
					$match += $1;
					$aln_length += $1;
					}
				 elsif ($cigar_part =~ /(\d+)=/)
					{
					$match += $1;
					$aln_length += $1;
					}
				 elsif ($cigar_part =~ /(\d+)I/)
					{
					$mismatch += $1;
					}
				 elsif ($cigar_part =~ /(\d+)D/)
					{
					$mismatch += $1;
					$aln_length += $1;
					}
				 elsif ($cigar_part =~ /(\d+)S/)
					{
					$mismatch += $1;
					}
				 elsif ($cigar_part =~ /(\d+)H/)
                    {
                    $mismatch += $1;
                    }
				 elsif ($cigar_part =~ /(\d+)X/)
                    {
                    $cigar_substitution += $1;
					$aln_length += $1;
                    }
				 else
				 	{
				 	die "Unexpected cigar: $cigar_unchanged\n";
				 	}
				 $cigar =~ s/$cigar_part//;
				 }
			elsif ($cigar eq "*")
				{
				#Unmapped Read
				$cigar =~ s/\*//
				}
			else
				{
				die "Unexpected cigar: $cigar_unchanged\n";
				}
			 }



			 $mismatch_all=0;

			 $md_col_unchanged = $md_col;
			 while ($md_col !~ /^$/)
			 	{
				if ($md_col =~ /^(\d+)/)
					{
					$md_col_part = $1;
					$md_col =~ s/^$md_col_part//;
					}
				elsif ($md_col =~ /^\^([a-zA-Z]+)/)
					{
					$md_col_part = $1;
					$md_col =~ s/^\^$md_col_part//;
					}
				elsif ($md_col =~ /^([a-zA-Z]+)/)
					{
					$mismatch_all += length($1);
					$md_col_part = $1;
					$md_col =~ s/^$md_col_part//;
					}
				elsif ($md_col eq "*")
					{
					# Unmapped Read
					$md_col =~ s/\*//
					}
				else
					{
					die "Unexpected MD: $md_col $md_col_unchanged\n";
					}
				}


			   if($mismatch_all != $cigar_substitution && $CIGAR_flag==1)
			   	{
			   	print STDERR "Substitutions from CIGAR and MD tag are discordant:\n\t$line[0] $cigar_unchanged $md_col_unchanged\n\t$cigar_substitution $mismatch_all\n";
			   	}

			 	$updated_chr = $line[2];
				@split_ID=split(/_/,$line[2]) if($bitflag[7] == 1);
				$updated_chr = $split_ID[0]."_RC_".join("_",splice(@split_ID,1)) if($bitflag[7] == 1);

				$score = "NA";
				$score_all = "NA";
				$aln_info = "NA";

			if($cigar ne "*" && $size > 0 && ($line[1]==0 || $BIT_flag==0))
				{

				$unaln_length=$size-$aln_length;
				$score = sprintf("%.3f", $mismatch/$size);
				$score_all = sprintf("%.3f", ($mismatch+$mismatch_all+$unaln_length)/$size) if($CIGAR_flag==0);
				$score_all = sprintf("%.3f", ($mismatch+$cigar_substitution+$unaln_length)/$size) if($CIGAR_flag==1);

				$aln_info=$line[3].":".$aln_length;

				#print join("\t",$id[0],$id[1],$bitflag[7],$line[2],$line[4],$size,$line[5],$score,$seq_correct_ori,"Secondary")."\n"  if($bitflag[3] == 1); #Secondary Alignments

				if($score_all <= $score_cutoff)
					{
					print MATCH join("\t",$id[0],$id[1],$bitflag[7],$updated_chr,$line[2],$line[4],$size,$line[5],$score_all,$seq_correct_ori,"PASS",$score,$md_col_unchanged,$aln_info)."\n"  if($bitflag[3] == 0);
					}
				else
					{
					print MATCH join("\t",$id[0],$id[1],$bitflag[7],$updated_chr,$line[2],$line[4],$size,$line[5],$score_all,$seq_correct_ori,"FAIL",$score,$md_col_unchanged,$aln_info)."\n"  if($bitflag[3] == 0);
					}
				}
			else
				{
				print MATCH join("\t",$id[0],$id[1],$bitflag[7],$updated_chr,$line[2],$line[4],$size,$line[5],"-",$seq_correct_ori,"FAIL",$score,$md_col_unchanged,$aln_info)."\n";
				print $seq_correct_ori."\n";
				print $score."\n";
				print $md_col_unchanged."\n";
				print $aln_info."\n";
				}
		}
	}
close MATCH;
