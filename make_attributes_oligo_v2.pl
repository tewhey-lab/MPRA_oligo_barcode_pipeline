#!/usr/bin/perl

## Create attributes file based on count data, and a list of controls
## Oligo example: chr:pos:ref:alt:allele:window(:haplotype)
## output columns: id snp chromosome snp_pos  ref_allele  alt_allele allele window  strand  project haplotype

use strict;
use warnings;

my $oligo_project=$ARGV[0];
my $out=$ARGV[1];

open (PROJ, "$oligo_project") or die("ERROR: can not read file ($oligo_project): $!\n");

open (OUT, ">$out".".attributes") or die("ERROR: can not open file ($out.attributes): $!\n");


my %oligo_proj;

my @line;
my @oligo_info;
my @attribute;
my @more;
my @loc;

my $oligo;
my $snp;
my $chr;
my $allele;
my $window;
my $strand;
my $haplotype;
my $ref_allele;
my $alt_allele;
my $snp_pos;

while (<PROJ>){
  chomp;
  @line = split(/\t/, $_);
  $oligo_proj{$line[0]} = $line[1];
}
close PROJ;

foreach $oligo (keys %oligo_proj)
  @oligo_info = split(/_/, $oligo)
  @attribute = split(/:/, $oligo_info[0]);
  @more = split(/-/, $oligo_info[1]);
  $snp = join(":", $attribute[0], $attribute[1], $attribute[2], $attribute[3]);
  $allele = 'ref' if($attribute[4] eq 'R');
  $allele = 'ref' if($attribute[4] eq 'A');
  $window = 'left' if($attribute[5] eq 'wL');
  $window = 'center' if($attribute[5] eq 'wC');
  $window = 'right' if($attribute[5] eq 'wR');
  $strand = 'fwd';
  $chr = $attribute[0];
  $ref_allele = $attribute[2];
  $alt_allele = $attribute[3];
  $snp_pos = $attribute[1];
  $haplotype = 'ref'
  if(scalar(@attribute)==6){
    $haplotype = 'alt' if($attribute[6] eq 'Alt')
  }

  print OUT join("\t", $oligo, $snp, $chr, $snp_pos, $ref_allele, $alt_allele, $allele,
                  $window, $strand, $oligo_proj{$oligo}, $haplotype . "\n");
}
