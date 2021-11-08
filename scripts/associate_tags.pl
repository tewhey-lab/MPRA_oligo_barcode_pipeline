#!/usr/bin/perl

##################
#
#052314 - Updated to account for multimapping oligos
#
####
#	FLAGS
##################

use strict;
use warnings;


my $tags = $ARGV[0]; #matched tag file
my $enhancers = $ARGV[1];
my $out = $ARGV[2];
my $read = $ARGV[3];

if ($enhancers =~ /.gz$/) {
open(ENHANCERS, "gunzip -c $enhancers |") or die("ERROR: can not open pipe to file ($enhancers): $!\n");
}
else {
open(ENHANCERS, "$enhancers") or die("ERROR: can not read file ($enhancers): $!\n");
}
open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");

my @inline;
my %tags; #[ct,orientation,location,tagseq(rc)]
my $revcomp;


open (TAGS, "$tags") or die("ERROR: can not read file ($tags): $!\n");
while (<TAGS>){

	chomp;
	@inline = split("\t");

  ${$tags{$inline[1]}}[0]++;
	${$tags{$inline[1]}}[1] = -9 ;
	${$tags{$inline[1]}}[2] = "-" ;
	${$tags{$inline[1]}}[3] = "-" ;
	${$tags{$inline[1]}}[4] = $inline[1];
	${$tags{$inline[1]}}[5] = "NA"; #mapping score
	${$tags{$inline[1]}}[6] = "NA"; #CIGAR
	${$tags{$inline[1]}}[7] = "NA"; #MD tag
	${$tags{$inline[1]}}[8] = "NA"; #start/stop pos

}
close TAGS;


my $cur_tag;
my $cur_loc;
my $cur_flag;
my $cur_m_flag;
my $cur_m_aln;
my $cur_m_cigar;
my $cur_m_md;
my $cur_m_pos;

my @tmp_id;
my @tmp_passflg;
my @tmp_aln;
my @tmp_cigar;
my @tmp_md;
my @tmp_pos;
my @tmp_ori;
my $collision_ok;
my $p;
my $i;

while (<ENHANCERS>){
	chomp;
	@inline = split("\t");

	$cur_tag = $inline[0];
	$cur_loc = $inline[1];
	$cur_flag = $inline[4];
	$cur_m_flag = $inline[5];
	$cur_m_aln = $inline[6];
	$cur_m_cigar = $inline[7];
	$cur_m_md = $inline[8];
	$cur_m_pos = $inline[9];

	if($cur_flag > 0)
		{
		if(exists($tags{$cur_tag}))
			{
			if($cur_flag == 1)
				{
				$collision_ok=0;
				@tmp_id = split(/,/,$cur_loc);
				@tmp_passflg = split(/,/,$cur_m_flag);
				@tmp_aln = split(/,/,$cur_m_aln);
				@tmp_cigar = split(/,/,$cur_m_cigar);
				@tmp_md = split(/,/,$cur_m_md);
				@tmp_pos = split(/,/,$cur_m_pos);


				if($collision_ok == 0){
					${$tags{$cur_tag}}[1] = -4;
					${$tags{$cur_tag}}[2] = $cur_loc;
					${$tags{$cur_tag}}[3] = $cur_flag;
					${$tags{$cur_tag}}[5] = $cur_m_aln;
					${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
					${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
					${$tags{$cur_tag}}[8] = $cur_m_pos; #Alignment start/stop

				}
				if($collision_ok == 1){
					for($i=0;$i<scalar(@tmp_id);$i++){
						${$tags{$cur_tag}}[1] = -1;
						${$tags{$cur_tag}}[2] = $tmp_id[$i];
						${$tags{$cur_tag}}[3] = $tmp_passflg[$i];
						${$tags{$cur_tag}}[5] = $tmp_aln[$i];
						${$tags{$cur_tag}}[6] = $tmp_cigar[$i]; #CIGAR
						${$tags{$cur_tag}}[7] = $tmp_md[$i]; #MD tag
						${$tags{$cur_tag}}[8] = $tmp_pos[$i]; #Alignment start/stop

					}
				}
			}
			elsif($cur_flag == 2){
				${$tags{$cur_tag}}[1] = -5;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = $cur_flag;
				${$tags{$cur_tag}}[5] = $cur_m_aln;
				${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
				${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
				${$tags{$cur_tag}}[8] = $cur_m_pos; #Alignment start/stop

			}

			elsif($cur_flag > 2){
				${$tags{$cur_tag}}[1] = -6;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = $cur_flag;
				${$tags{$cur_tag}}[5] = $cur_m_aln;
				${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
				${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
				${$tags{$cur_tag}}[8] = $cur_m_pos; #Alignment start/stop

			}
		}
	}
	else{
		if(exists($tags{$cur_tag})){
			${$tags{$cur_tag}}[1] = 0;
			${$tags{$cur_tag}}[2] = $cur_loc;
			${$tags{$cur_tag}}[3] = $cur_flag;
			${$tags{$cur_tag}}[5] = $cur_m_aln;
			${$tags{$cur_tag}}[6] = $cur_m_cigar;
			${$tags{$cur_tag}}[7] = $cur_m_md;
			${$tags{$cur_tag}}[8] = $cur_m_pos;

		}
	}
}
close ENHANCERS;


my $key;

foreach $key (keys %tags)
	{
	print OUT join("\t",$key,@{$tags{$key}})."\n";
	}
