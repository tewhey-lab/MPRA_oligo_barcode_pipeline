#!/usr/bin/perl

## The FASTA input should be your reference FASTA used in the MPRAMatch.wdl pipeline
## The FASTA should be broken into groups for each "project" i.e. negative controls, positive controls, experimental
## Any oligos which belong to more than one group should processed separately, and the two projects should be separated by a comma

use strict;
use warnings;

my $oligo_seqs = $ARGV[0];
my $project_name = $ARGV[1];

open (FASTA, "$oligo_seqs") or die("ERROR: can not read file ($oligo_seqs): $!\n");
open (OUT, ">$project_name".".proj_list") or die("ERROR: can not open file ($project_name.proj_list): $!\n");

my @header;
my @array;
my $id;

while (<FASTA>) {
  if(/^>/){
    @header = split(/\s/, $_);
    @array = split(/\|/, $header[0]);
    $id = $array[0];
    $id =~ s/^@/>/;
    $id =~ s/^>//;
    $id =~ s/\/1$//;

    print OUT $id . "\t" . $project_name . "\n";
  }
}
