#!/usr/bin/perl

## Create attributes file based on count data, and a list of controls
## Oligo example: chr:pos:ref:alt_allele_window_strand
## output columns: id  snp  chromosome  snp_pos ref_allele  alt_allele  allele  window  strand  start stop  project

use strict;
use warnings;

my $oligo_project=$ARGV[0];
my $out=$ARGV[1];

open (PROJ, "$oligo_project") or die("ERROR: can not read file ($oligo_project): $!\n");

open (OUT, ">$out".".attributes") or die("ERROR: can not open file ($out.attributes): $!\n");

my %oligo_proj;

my @line;
my @attribute;
my @more;
my @loc;

my $oligo;
my $id;
my $snp;
my $chr;
my $allele;
my $window;
my $strand;
my $loc_info;
my $ref_allele;
my $alt_allele;
my $start;
my $stop;
my $snp_pos;
my $position;

while (<PROJ>){
  chomp;
  @line = split(/\t/, $_);
  $oligo_proj{$line[0]} = $line[1];
}
close PROJ;

foreach $oligo (keys %oligo_proj){
  @attribute = split(/_/, $oligo);
  $snp = $attribute[0];
  $allele = 'NA';
  $window = 'NA';
  $strand = 'fwd';
  $chr = 'NA';
  $ref_allele = 'NA';
  $alt_allele = 'NA';
  $start = 'NA';
  $stop = 'NA';
  $snp_pos = 'NA';
  if(scalar(@attribute)==3){
    $allele = 'ref' if($attribute[1] eq 'A');
    $allele = 'alt' if($attribute[1] eq 'B');
    $window = 'center' if($attribute[2] eq 'wC');
    $strand = 'rev' if($attribute[3] eq 'RC');
  }
  if(scalar(@attribute)==2){
    $allele = 'ref' if($attribute[1] eq 'A');
    $allele = 'alt' if($attribute[1] eq 'B');
    $window = 'center' if($attribute[2] eq 'wC');
  }
  if(scalar(@attribute)==1){
    $allele = 'ref' if($attribute[1] eq 'A');
    $allele = 'alt' if($attribute[1] eq 'B');
    $strand = 'rev' if($attribute[1] eq 'RC');
  }
  @more = split(/:/, $snp);
  if(scalar(@more) ==3){
    $chr = $more[0];
    $loc_info = $more[1];
    $ref_allele = $more[2];
    $alt_allele = $more[3];
  }
  else{
    $chr = $more[0];
    $loc_info = $more[1];
  }
  @loc = split(/-/, $loc_info);
  if(scalar(@loc)==1){
    $start = $loc[0];
    $stop = $loc[1];
  }
  if(scalar(@loc)==0){
    $snp_pos = $loc[0];
  }
  print OUT join("\t", $oligo, $snp, $chr, $snp_pos, $ref_allele, $alt_allele, $allele,
                  $window, $strand, $start, $stop, $oligo_proj{$oligo} . "\n");
}
