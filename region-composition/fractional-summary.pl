#!/usr/bin/perl

# 27/7/2017

use strict;
use warnings;

# Libraries.
use Statistics::Basic qw(mean median);
use Statistics::Robust::Scale qw(MAD);
use Text::CSV_XS;

# Parameters.
my $regionOfInterest = 'disorder';

# Command line parameters.
my $bfile = ''; # background list of genes
my $lfile = ''; # list of genes to summarize
my $pfile = ''; # file with complete pfam annotations

if ($#ARGV==0){
	print "After running pfam_format.pl, this script will take the _complete output\n";
  print "along with a list of genes and a motif name, and output stats for that list.\n\n";
	print "\nusage:\n $0\n";
  print "-b [background list with gene name in first column and uniprot ID in second]\n";
	print "-l [gene list with gene name in first column and uniprot ID in second]\n";
  print "-m [motif of interest [default: 'disorder']]\n";
  print "-p [file with Pfam motifs]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$regionOfInterest = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-p'){
			$i++;
			$pfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

sub isInBackground {
  my $gene = $_[0];
  my $hasBackground = $_[1];
  my %background = %{$_[2]};

  if (
    !$hasBackground
    || ($gene && exists $background{$gene})
  ) {
    return 1;
  }

  return 0;
}

sub printWilcox {
  my $outfile = $_[0];
  my @arr = @{$_[1]};
  open my $FH, '>', $outfile;
  foreach my $value (@arr) {
    print $FH "$value\n";
  }
  close $FH;
}

my %tsvParams = (
	sep_char => "\t",
);

# Read background list.
my %backgroundGenes;
if ($bfile) {
  print STDERR "Reading background list: $bfile\n";
  open my $backgroundFH, '<', $bfile or die "Could not open $bfile: $!";
  my $backgroundTSV = Text::CSV_XS->new(\%tsvParams);
  while(my $row = $backgroundTSV->getline($backgroundFH)) {
    $backgroundGenes{@{$row}[1]} = 1
  }
  close $backgroundFH;
}

# Read gene list.
print STDERR "Reading gene list: $lfile\n";
open my $listFH, '<', $lfile or die "Could not open $lfile: $!";
my $listTSV = Text::CSV_XS->new(\%tsvParams);
my @queryGenes;
while(my $row = $listTSV->getline($listFH)) {
  push @queryGenes, @{$row}[1];
}
close $listFH;

# Read motifs.
print STDERR "Parsing Pfam result\n";
open my $pfamFH, '<', $pfile or die "Could not open $pfile: $!";
my $pfamTSV = Text::CSV_XS->new(\%tsvParams);
## Skip header
$pfamTSV->getline($pfamFH);
my $lastGene = '';
my $motifTotalLength;
my %pfam;
while(my $row = $pfamTSV->getline($pfamFH)) {
  my $currGene = @{$row}[1];
  my $currMotif = @{$row}[2];
  my $motifLength = @{$row}[5];
  my $proteinLength = @{$row}[6];
  if (isInBackground($currGene, $bfile, \%backgroundGenes)) {
    if ($currGene ne $lastGene) {
      $lastGene = $currGene;
      $motifTotalLength = 0;
      $pfam{$currGene}{'fraction'} = 0;
      $pfam{$currGene}{'longest'} = 0;
      $pfam{$currGene}{'total'} = 0;
    }
    if ($currMotif eq $regionOfInterest) {
      $motifTotalLength += $motifLength;
      $pfam{$currGene}{'fraction'} = sprintf "%.3f", $motifTotalLength / $proteinLength;
      if (
        !(exists $pfam{$currGene}{'longest'})
        || $motifLength > $pfam{$currGene}{'longest'}
      ) {
        $pfam{$currGene}{'longest'} = $motifLength;
      }
      $pfam{$currGene}{'total'} = $motifTotalLength;
    }
  }
}
close $pfamFH;

# Summarize results.
print STDERR "Summarizing results\n";
my %summary;
foreach my $gene (keys %pfam) {
  push @{$summary{'all'}{'fraction'}}, $pfam{$gene}{'fraction'};
  push @{$summary{'all'}{'longest'}}, $pfam{$gene}{'longest'};
  push @{$summary{'all'}{'total'}}, $pfam{$gene}{'total'};
  if (grep(/^$gene$/, @queryGenes)) {
    push @{$summary{'query'}{'fraction'}}, $pfam{$gene}{'fraction'};
    push @{$summary{'query'}{'longest'}}, $pfam{$gene}{'longest'};
    push @{$summary{'query'}{'total'}}, $pfam{$gene}{'total'};
  }
}

# Create histogram.
my $binmin = 0;
my $binwidth = 0.05;
my %histo;
foreach my $value (@{$summary{'all'}{'fraction'}}) {
	my $bin = int( ( $value - $binmin ) / $binwidth );
	$histo{'all'}[ $bin ]++;
}
foreach my $value (@{$summary{'query'}{'fraction'}}) {
	my $bin = int( ( $value - $binmin ) / $binwidth );
	$histo{'query'}[ $bin ]++;
}
open my $histoFH, '>', 'output/histograms.txt';
for(my $i = 0; $i <= 20; $i++) {
	my $valueAll = 0;
	my $valueQuery = 0;
	if ($histo{'all'}[$i]) {
		$valueAll = $histo{'all'}[$i];
	}
	if ($histo{'query'}[$i]) {
		$valueQuery = $histo{'query'}[$i];
	}
	print $histoFH "$i\t$valueAll\t$valueQuery\n";
}

# Calculate stats.
my %stats;
$stats{'all'}{'mean'}{'fraction'} = mean @{$summary{'all'}{'fraction'}};
$stats{'all'}{'mean'}{'longest'} = mean @{$summary{'all'}{'longest'}};
$stats{'all'}{'mean'}{'total'} = mean @{$summary{'all'}{'total'}};
$stats{'all'}{'median'}{'fraction'} = median @{$summary{'all'}{'fraction'}};
$stats{'all'}{'median'}{'longest'} = median @{$summary{'all'}{'longest'}};
$stats{'all'}{'median'}{'total'} = median @{$summary{'all'}{'total'}};
$stats{'all'}{'mad'}{'fraction'} = MAD(\@{$summary{'all'}{'fraction'}});
$stats{'all'}{'mad'}{'longest'} = MAD(\@{$summary{'all'}{'longest'}});
$stats{'all'}{'mad'}{'total'} = MAD(\@{$summary{'all'}{'total'}});
$stats{'query'}{'mean'}{'fraction'} = mean @{$summary{'query'}{'fraction'}};
$stats{'query'}{'mean'}{'longest'} = mean @{$summary{'query'}{'longest'}};
$stats{'query'}{'mean'}{'total'} = mean @{$summary{'query'}{'total'}};
$stats{'query'}{'median'}{'fraction'} = median @{$summary{'query'}{'fraction'}};
$stats{'query'}{'median'}{'longest'} = median @{$summary{'query'}{'longest'}};
$stats{'query'}{'median'}{'total'} = median @{$summary{'query'}{'total'}};
$stats{'query'}{'mad'}{'fraction'} = MAD(\@{$summary{'query'}{'fraction'}});
$stats{'query'}{'mad'}{'longest'} = MAD(\@{$summary{'query'}{'longest'}});
$stats{'query'}{'mad'}{'total'} = MAD(\@{$summary{'query'}{'total'}});

# print arrays for wilcox test
printWilcox('output/wilcox_all_fraction.txt', \@{$summary{'all'}{'fraction'}});
printWilcox('output/wilcox_query_fraction.txt', \@{$summary{'query'}{'fraction'}});
printWilcox('output/wilcox_all_longest.txt', \@{$summary{'all'}{'longest'}});
printWilcox('output/wilcox_query_longest.txt', \@{$summary{'query'}{'longest'}});

# Output summary.
open my $summaryFH, '>', 'output/' . $regionOfInterest . '-summary.txt';

## Mean.
print $summaryFH "mean\tfractional protein length for region\tlongest region\ttotal region AA\tmad fraction\tmad longest\n";
print $summaryFH "all\t$stats{'all'}{'mean'}{'fraction'}\t$stats{'all'}{'mean'}{'longest'}\t$stats{'all'}{'mean'}{'total'}\t$stats{'all'}{'mad'}{'fraction'}\t$stats{'all'}{'mad'}{'longest'}\n";
print $summaryFH "query\t$stats{'query'}{'mean'}{'fraction'}\t$stats{'query'}{'mean'}{'longest'}\t$stats{'query'}{'mean'}{'total'}\t$stats{'query'}{'mad'}{'fraction'}\t$stats{'query'}{'mad'}{'longest'}\n\n";

## Median.
print $summaryFH "median\tfractional protein length for region\tlongest region\ttotal region AA\n";
print $summaryFH "all\t$stats{'all'}{'median'}{'fraction'}\t$stats{'all'}{'median'}{'longest'}\t$stats{'all'}{'median'}{'total'}\n";
print $summaryFH "query\t$stats{'query'}{'median'}{'fraction'}\t$stats{'query'}{'median'}{'longest'}\t$stats{'query'}{'median'}{'total'}\n";
close $summaryFH;
