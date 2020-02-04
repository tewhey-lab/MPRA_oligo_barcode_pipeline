#!/usr/bin/perl

use strict;
use warnings;

my $fasta = $ARGV[0];
my $read = $ARGV[1];
my $out = $ARGV[2];
my $link_A_bc = $ARGV[3];
my $link_A_oligo = $ARGV[4];
my $end_A_oligo = $ARGV[5];
my $MIN_SEQ_SIZE = $ARGV[6];

open (FASTA, "$fasta") or die("ERROR: can not read file ($fasta): $!\n");
open (MATCH, ">$out".".match") or die("ERROR: can not create $out .matched: $!\n");
open (REJECT, ">$out".".reject") or die("ERROR: can not create $out .rejected: $!\n");

# my $link_A_bc = "TCTAGA";
# my $MIN_SEQ_SIZE = 100;
my $barcode_seq;
my $barcode_start;
my $id;
my $r1;
my $revcomp;
my $revcomp_barcode;
my $link_index;
my $oligo_start;
my $oligo_end;
my $oligo_length;
my $oligo_seq;
my $revcomp_oligo;


while (<FASTA>){
# Extract Sequence ID
  chomp;
  $id = $_;
  $id =~ s/^@/>/;
  $id =~ s/^>//;
  $id =~ s/\/1$//;

# Extract the sequence
  $r1 = <FASTA>;
  chomp $r1;
  if(length($r1) < $MIN_SEQ_SIZE){
    next;
  }

# If checking a Read 1 file take the reverse complement
  if($read == 1){
    $revcomp = reverse($r1);
  	$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
    $r1 = $revcomp;
  }
# Check for presence of linker A
  if(index(substr($r1, 18, 10), $link_A_bc) != -1){
    $link_index = index(substr($r1, 18, 10), $link_A_bc);
    $link_index += 18;
  }
  if(index(substr($r1, 18, 10), $link_A_bc) == -1){
    print REJECT "$id\n";
  }

  if($link_index - 20 <= 0){
    $barcode_start = 0;
    $barcode_seq = substr($r1, $barcode_start, 20);
  }
  if($link_index - 20 > 0){
    $barcode_start = $link_index - 20;
    $barcode_seq = substr($r1, $barcode_start, 20);
  }

  $oligo_start = index(substr($r1, $link_index+30, 12), $link_A_oligo);
  if($oligo_start == -1){
    $oligo_start = 0;
  }
  $oligo_start += 34+$link_index;
# Find the end of the oligo
  $oligo_end = index(substr($r1, -18, 6), $end_A_oligo);
  $oligo_end += -18;

# Define the substring that is the oligo
  $oligo_length = length($r1) + $oligo_end - $oligo_start;
  $oligo_seq = substr($r1, $oligo_start, $oligo_length);
  $revcomp_oligo = reverse($oligo_seq);
  $revcomp_oligo =~ tr/ACGTNacgtn/TGCANtgcan/;
  $oligo_seq = $revcomp_oligo;


  if($id ne "+"){
    print MATCH join("\t", $id, $barcode_seq, $oligo_seq, $oligo_length."\n");
  }
}
