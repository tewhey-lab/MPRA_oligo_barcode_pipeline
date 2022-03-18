#!/usr/bin/perl

## Create attributes file based on a project file containing the oligo ID and the project(s) to which it belongs
## Oligo example: chr:pos:ref:alt:allele:window(:haplotype)
## output columns: id snp chromosome snp_pos  ref_allele  alt_allele allele window  strand  project haplotype

use strict;
use warnings;

# The oligo project should be the combined project lists from make_project_list.pl also available in the scripts directory
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
my $len;
my $i;
my $field;

while (<PROJ>){
  chomp;
  @line = split(/\t/, $_);
  $oligo_proj{$line[0]} = $line[1];
}
close PROJ;

print OUT join("\t","ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","project","haplotype")."\n";

foreach $oligo (keys %oligo_proj) {

  @attribute = split(/:/, $oligo);

  $len=scalar(@attribute);
  die "Need atleast 4 attributes in ID. Only $len found in $oligo\n" if($len < 4);

  $chr = $attribute[0];
  $snp_pos = $attribute[1];
  $ref_allele = $attribute[2];
  $alt_allele = $attribute[3];

  if($len >= 5) {
    $allele = $attribute[4];
    $allele = 'ref' if($attribute[4] eq 'R');
    $allele = 'alt' if($attribute[4] eq 'A');
    }
  else {$allele = "NA";}
  print STDERR "Allele should be R or A, set as '$allele' in $oligo\n" unless($allele eq 'ref' || $allele eq 'alt');


  if($len >= 6) {
    $window = $attribute[5];
    $window = 'left' if($attribute[5] eq 'wL');
    $window = 'center' if($attribute[5] eq 'wC');
    $window = 'right' if($attribute[5] eq 'wR');
    }
  else {$window = "NA";}
  print STDERR "Window should be wL, wC or wR, set as '$window' in $oligo\n" unless($window eq 'left' || $window eq 'center' || $window eq 'right');

  $snp = join(":", $attribute[0], $attribute[1], $attribute[2], $attribute[3]);
  $strand = 'fwd';
  $haplotype = 'ref';

  ###Process additional attributes
  if(scalar(@attribute) > 5) {
    for($i=5;$i<scalar(@attribute);$i++) {
      if($attribute[$i] =~ m/^Alt/)
      	{
    	###Process alt haplotypes
    	$haplotype = 'alt'
    	}
    }
  }
  print OUT join("\t", $oligo, $snp, $chr, $snp_pos, $ref_allele, $alt_allele, $allele,
                  $window, $strand, $oligo_proj{$oligo}, $haplotype . "\n");
}
