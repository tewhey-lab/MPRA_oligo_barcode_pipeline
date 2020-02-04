### filter out blacklisted barcodes from count data

#!/usr/bin/perl

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

    print OUT $id . "\t" . $project_name . "\n";
  }
}
