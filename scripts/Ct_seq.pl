#!/usr/bin/perl
#
#by Ryan Tewhey
use strict;
use warnings;
##################################
#
##################################

my $mapped = $ARGV[0];
my $CT_COL = $ARGV[1]-1;
my $CMP_COL = $ARGV[2]-1;

open (MAPPED, "$mapped") or die("ERROR: can not read file (mapped): $!\n");

my $first = 0;
my $tmp_line;
my @last;

my $next_line;
my @next;

my $cur_barcode;
my %cur_hits;

my %cur_pass_flag;
my %cur_hits_score;
my %cur_cigar;
my %cur_mdtag;
my %cur_pos;

my @print_score;
my @tmp_sort_idx;
my @tmp_sort;
my @print_keys;
my @print_values;
my @print_passflag;
my @print_cigar;
my @print_md;
my @print_pos;

my $key;
my $sum;
my $ct_pass_flag; # 0 = ok, 1 = collision, 2 = other failure (mapping)
my $cmp_pass_flag; #0 = ok, 2 = mapping

my $aln_score;
my $cigar;
my $mdtag;
my $pos;

while (<MAPPED>)
	{
	chomp;
	my @line = split(/\t/);

	if($first == 0)
		{
		@last = @line;
		$first = 1;
		next;
		}

	$ct_pass_flag = 2;
	$cmp_pass_flag = 2;
	$ct_pass_flag = 0 if($last[10] eq "PASS");
	$cmp_pass_flag = 0 if($last[10] eq "PASS");

	$cur_barcode = $last[$CT_COL];
	$cur_hits{$last[$CMP_COL]}++;
	$aln_score=$last[8];
	$aln_score=1 if($last[8] eq "-");
	push(@{$cur_hits_score{$last[$CMP_COL]}},$aln_score);
	push(@{$cur_pass_flag{$last[$CMP_COL]}},$cmp_pass_flag);
	$cigar = $last[7];
	$mdtag = $last[12];
	$pos = $last[13];
	push(@{$cur_cigar{$last[$CMP_COL]}},$cigar);
	push(@{$cur_mdtag{$last[$CMP_COL]}},$mdtag);
	push(@{$cur_pos{$last[$CMP_COL]}},$pos);

	while($last[$CT_COL] eq $line[$CT_COL] && !(eof))
		{
		$cmp_pass_flag = 2;
		$cmp_pass_flag = 0 if($line[10] eq "PASS");
		$aln_score=$line[8];
		$aln_score=1 if($line[8] eq "-");

		$cur_hits{$line[$CMP_COL]}++;
		push(@{$cur_hits_score{$line[$CMP_COL]}},$aln_score);
		push(@{$cur_pass_flag{$line[$CMP_COL]}},$cmp_pass_flag);
		$ct_pass_flag = 0 if($line[10] eq "PASS");

		$cigar = $line[7];
		$mdtag = $line[12];
		$pos = $line[13];
		push(@{$cur_cigar{$line[$CMP_COL]}},$cigar);
		push(@{$cur_mdtag{$line[$CMP_COL]}},$mdtag);
		push(@{$cur_pos{$line[$CMP_COL]}},$pos);

		$tmp_line=<MAPPED>;
		chomp $tmp_line;
		@line = split(/\t/,$tmp_line);
		}

	if(eof) ##If end of file process last two lines
		{
		if($line[$CT_COL] eq $last[$CT_COL])
			{
			$cur_hits{$line[$CMP_COL]}++;
			push(@{$cur_hits_score{$line[$CMP_COL]}},$aln_score);
			push(@{$cur_pass_flag{$line[$CMP_COL]}},$cmp_pass_flag);
			push(@{$cur_cigar{$line[$CMP_COL]}},$cigar);
			push(@{$cur_mdtag{$line[$CMP_COL]}},$mdtag);
			push(@{$cur_pos{$line[$CMP_COL]}},$pos);

			$ct_pass_flag = 0 if($line[10] eq "PASS");
			}
		else  ## Run normal loop on last line if unique
			{
				$sum = 0;
				@print_keys=();
				@print_values=();
				@print_score=();
				@print_passflag=();
				@print_cigar=();
				@print_md=();
				@print_pos=();

				for $key (keys %cur_hits)
					{
					push(@print_keys,$key);
					push(@print_values, $cur_hits{$key});
  					$sum += $cur_hits{$key};

  					#Pick the oligo with the best alignment score.
  					@tmp_sort_idx = 0..$#{$cur_hits_score{$key}};
  					@tmp_sort_idx = sort { ${$cur_hits_score{$key}}[$a] <=> ${$cur_hits_score{$key}}[$b] } 0..$#{$cur_hits_score{$key}} if(scalar @{$cur_hits_score{$key}} > 1);

  					#print STDERR join (",",$key,@{$cur_pass_flag{$key}})."\n";
  					@tmp_sort = @{$cur_pass_flag{$key}}[@tmp_sort_idx];
  					#@tmp_sort = sort {$a <=> $b} @{$cur_pass_flag{$key}} if(scalar @{$cur_pass_flag{$key}} > 1);
  					push(@print_passflag, $tmp_sort[0]);

  					#print STDERR join (",",$key,@{$cur_hits_score{$key}})."\n";
  					@tmp_sort = @{$cur_hits_score{$key}}[@tmp_sort_idx];
  					#@tmp_sort = sort {$a <=> $b} @{$cur_hits_score{$key}} if(scalar @{$cur_hits_score{$key}} > 1);
  					push(@print_score,sprintf("%.3f",$tmp_sort[0]));

  					@tmp_sort = @{$cur_cigar{$key}}[@tmp_sort_idx];
  					push(@print_cigar, $tmp_sort[0]);

  					@tmp_sort = @{$cur_mdtag{$key}}[@tmp_sort_idx];
  					push(@print_md, $tmp_sort[0]);

  					@tmp_sort = @{$cur_pos{$key}}[@tmp_sort_idx];
  					push(@print_pos, $tmp_sort[0]);
					}
				print "$cur_barcode\t";
				print join(",",@print_keys)."\t";
				print join(",",@print_values)."\t";
				print $sum."\t";

				$ct_pass_flag = 1 if(scalar(keys %cur_hits) > 1);
				print $ct_pass_flag."\t";
				print join(",",@print_passflag)."\t";
				print join(",",@print_score)."\t";
				print join(",",@print_cigar)."\t";
				print join(",",@print_md)."\t";
				print join(",",@print_pos)."\n";

				%cur_hits = ();
				%cur_hits_score = ();
				%cur_pass_flag = ();
				%cur_cigar = ();
				%cur_mdtag = ();
				%cur_pos = ();

				@last = @line;
				$ct_pass_flag = 2;
				$ct_pass_flag = 0 if($last[10] eq "PASS");
				$cmp_pass_flag = 2;
				$cmp_pass_flag = 0 if($last[10] eq "PASS");

				$cur_barcode = $last[$CT_COL];
				$cur_hits{$last[$CMP_COL]}++;
				$aln_score=$last[8];
				$aln_score=1 if($last[8] eq "-");
				push(@{$cur_hits_score{$last[$CMP_COL]}},$aln_score);
				push(@{$cur_pass_flag{$last[$CMP_COL]}},$cmp_pass_flag);

				$cigar = $last[7];
				$mdtag = $last[12];
				$pos = $last[13];
				push(@{$cur_cigar{$last[$CMP_COL]}},$cigar);
				push(@{$cur_mdtag{$last[$CMP_COL]}},$mdtag);
				push(@{$cur_pos{$last[$CMP_COL]}},$pos);

			}

		}



	$sum = 0;
	@print_keys=();
	@print_values=();
	@print_score=();
	@print_passflag=();
	@print_cigar=();
	@print_md=();
	@print_pos=();

	for $key (keys %cur_hits)
		{
		push(@print_keys,$key);
		push(@print_values, $cur_hits{$key});
  		$sum += $cur_hits{$key};

  		#Pick the oligo with the best alignment score.
  		@tmp_sort_idx = 0..$#{$cur_hits_score{$key}};
  		@tmp_sort_idx = sort { ${$cur_hits_score{$key}}[$a] <=> ${$cur_hits_score{$key}}[$b] } 0..$#{$cur_hits_score{$key}} if(scalar @{$cur_hits_score{$key}} > 1);

  		#print STDERR join (",",$key,@{$cur_pass_flag{$key}})."\n";
  		@tmp_sort = @{$cur_pass_flag{$key}}[@tmp_sort_idx];
  		#@tmp_sort = sort {$a <=> $b} @{$cur_pass_flag{$key}} if(scalar @{$cur_pass_flag{$key}} > 1);
  		push(@print_passflag, $tmp_sort[0]);

  		#print STDERR join (",",$key,@{$cur_hits_score{$key}})."\n";
  		@tmp_sort = @{$cur_hits_score{$key}}[@tmp_sort_idx];
  		#@tmp_sort = sort {$a <=> $b} @{$cur_hits_score{$key}} if(scalar @{$cur_hits_score{$key}} > 1);
  		push(@print_score,sprintf("%.3f",$tmp_sort[0]));

  		@tmp_sort = @{$cur_cigar{$key}}[@tmp_sort_idx];
  		push(@print_cigar, $tmp_sort[0]);

  		@tmp_sort = @{$cur_mdtag{$key}}[@tmp_sort_idx];
  		push(@print_md, $tmp_sort[0]);

		@tmp_sort = @{$cur_pos{$key}}[@tmp_sort_idx];
		push(@print_pos, $tmp_sort[0]);
		}
	print "$cur_barcode\t";
	print join(",",@print_keys)."\t";
	print join(",",@print_values)."\t";
	print $sum."\t";

	$ct_pass_flag = 1 if(scalar(keys %cur_hits) > 1);
	print $ct_pass_flag."\t";
	print join(",",@print_passflag)."\t";
	print join(",",@print_score)."\t";
	print join(",",@print_cigar)."\t";
	print join(",",@print_md)."\t";
	print join(",",@print_pos)."\n";



	%cur_hits = ();
	%cur_hits_score = ();
	%cur_pass_flag = ();
	%cur_cigar = ();
	%cur_mdtag = ();
	%cur_pos = ();
	@last = @line;
	}
