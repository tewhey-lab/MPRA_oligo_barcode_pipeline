#!/usr/bin/perl
#
#by Ryan Tewhey
use strict;
use warnings;
use Getopt::Std;


my %options=();
getopts('SA:', \%options);


#####
#
#-S = Run in saturation mutagenesis mode
#-A = Atributes file
#####

my $sat_mut_mode;
my $attributes = "NA";

$attributes = $options{A} || "NA";

if(exists($options{S}))
	{
	print STDERR "Running in saturation mutagenesis mode\nWith $attributes as attribute file\n";
	$sat_mut_mode = "T";
	}
else {$sat_mut_mode = "F";}


my $mapped = $ARGV[0];

open (ATTRIBUTES, "$attributes") or die("ERROR: can not read file (attributes): $!\n") if($sat_mut_mode eq "T");
open (MAPPED, "$mapped") or die("ERROR: can not read file (mapped): $!\n");

my @line;

my @ids;
my @cov;
my @pass;
my @aln;
my @cigar;
my @md;
my @pos;

my $max;
my $i;
my $keep;
my $max_idx;

my $head = 1;
my %ref_hash;
my $sat_ref_col;
my $ID_col;
my @id_arry;
my @max_id_arry;

my $z;
my $y;

if($sat_mut_mode eq "T")
	{
	while (<ATTRIBUTES> )
		{
				
		chomp;
		my @line = split(/\t/);
		
		if($head == 1)
			{
			for($i=0;$i<scalar(@line);$i++)
				{
				$sat_ref_col = $i if($line[$i] eq "sat_ref_parent");
				$ID_col = $i if($line[$i] eq "ID");
				}
			
			print STDERR "ID Col: $ID_col\nSat Ref Col: $sat_ref_col\n";
			$head++;
			next;
			}

		$ref_hash{$line[$ID_col]}=$line[$sat_ref_col];	
		}
	}
	
while (<MAPPED>)
	{
	chomp;
	my @line = split(/\t/);
	
	my @ids = split(/,/,$line[1]);
	my @cov = split(/,/,$line[2]);
	my @pass = split(/,/,$line[5]);
	my @aln = split(/,/,$line[6]);
	my @cigar = split(/,/,$line[7]);
	my @md = split(/,/,$line[8]);
	my @pos = split(/,/,$line[9]);

	if(scalar(@ids) > 1)
		{
		$max = 0;
		$keep = "0"; # 0=Modify to use indiv hit; 1=keep original flag
		
		for (my $i = 0; $i < scalar(@ids); $i++) 
			{
  			if ($cov[$i] > $max) 
  				{
  				$max = $cov[$i];
   				$max_idx = $i;
  				}
			}
			
		for($i=0;$i<scalar(@cov);$i++)
			{
			if($cov[$i] == $max && $i ne $max_idx)
				{
				$keep = 1;
				}
			elsif($cov[$i] == $max && $max_idx eq $i) #legacy block keeping as a sanity check if needed
				{
				#$max_idx = $i; 
				}
			elsif($sat_mut_mode eq "T")
				{

				if($cov[$i]/$max > 0.1 && $aln[$i]*1 == 0 && $aln[$max_idx]*1 == 0) #High quality mapping - toss
					{
					$keep = 1; 
					}
				elsif($cov[$i]/$max > 0) #non-perfect and on same reference - keep
					{
					@id_arry=split_ID($ids[$i]);
					@max_id_arry=split_ID($ids[$max_idx]);
					
					for($z=0;$z<scalar(@id_arry);$z++)
						{

						for($y=0;$y<scalar(@max_id_arry);$y++)
							{
							if($id_arry[$z] ne "*" && $max_id_arry[$y] ne "*" && $ref_hash{$id_arry[$z]} eq $ref_hash{$max_id_arry[$y]})
								{
								if($cov[$i]/$max > 0.5)
									{
									$keep = 1; #same reference mapping. Loosen cutoffs. 
									}
								}
							elsif($cov[$i]/$max > 0.1)  #non-perfect and different reference - toss
								{
								$keep = 1
								}
							else
								{
								#print join("\t",$pass[$max_idx],$i,$max_idx,$aln[$i],$aln[$max_idx],$cov[$i],$max,$keep,$line[0],join(",",@ids),join(",",@cov),join(",",@aln))."\n";
								}
							}
						}
					}
				}
			elsif($sat_mut_mode eq "F" && $cov[$i]/$max > .1)
				{
				$keep = 1;
				}
			else
				{
				#others
				}
			}	
		print join("\t",$line[0],$ids[$max_idx],$cov[$max_idx],$line[3],$pass[$max_idx],$pass[$max_idx],$aln[$max_idx],$cigar[$max_idx],$md[$max_idx],$pos[$max_idx])."\n" if($keep==0);	
		print join("\t",@line)."\n" if($keep==1);	

		}
	else
		{
		print join("\t",@line)."\n";	
		}
	}
	
	
sub split_ID
{
my ($full_id) = @_;
$full_id =~ s/^\(//;
$full_id =~ s/\)$//;
my @split_ids=split(/;/,$full_id);
return(@split_ids)
}


