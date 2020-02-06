#!/usr/bin/perl

## Create attributes file based on count data, and a list of controls
## Oligo example: chr:pos:ref:alt_allele_window_strand
## output columns: id  snp allele  window  strand  haplotype project

use strict;
use warnings;

my $oligo_project=$ARGV[0];
my $out=$ARGV[1];

open (PROJ, "$oligo_project") or die("ERROR: can not read file ($oligo_project): $!\n");

open (OUT, ">$out".".attributes") or die("ERROR: can not open file ($out.attributes): $!\n");

my %oligo_proj;

my @line;
my @more

my $snp;
my $id
my $oligo_num;
my $oligo_name;
my $allele;
my $strand = 'fwd';
my $window = 'center';
my $haplotype = 'ref';

while (<PROJ>){
  chomp;
  @line = split(/\t/, $_);
  $oligo_proj{$line[0]} = $line[1];
}
close PROJ;

foreach $id (keys %oligo_proj)
  @more = split(/_/, $id);
  if(scalar(@more)==1){
    $snp = $more[0];
    $allele = 'ref' if($more[1] eq 'A');
    $allele = 'alt' if($more[1] eq 'B');
  }
  if(scalar(@more)==2){
    if($more[2] =~ /(\d+)/){
      $oligo_num =~ s/(\d+)//;
      $oligo_name =~ s/\D//g;
      $allele = 'ref';
      if($oligo_num % 2 ==0){
        $oligo_num += -1;
        $allele= 'alt'
      }
      $snp = $more[0]."_".$more[1]."_".$oligo_name.$oligo_num;
    }
    else{
      $snp = $more[1];
      $allele = $more[2];
    }
  }
  print OUT join("\t", $id, $snp, $allele, $window, $strand, $haplotype, $oligo_project{$id} . "\n");
}
