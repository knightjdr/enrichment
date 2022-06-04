#!/usr/bin/perl

# 7/7/2017

use strict;
use warnings;

# modules
use FindBin;
use lib "$FindBin::RealBin/../lib"; 

use List::MoreUtils qw(uniq);
use Stats::BenjaminiHochberg qw(bhCorrection);
use Stats::Fishers qw(fishers);
use String::Util qw(trim);
use Text::CSV_XS;

# parameters
my $allowMissing = 0; # should genes not found in database be allowed to affect enrichment
my $fdr = 0.01;
my $outputPrefix = 'enrichment';

# command line parameters
my $bfile = ''; # background file
my $gfile = ''; # list of genes to test for enrichment
my $tfile = '';	# file with list of enrichment terms for each gene

if ($#ARGV == 0) {
	print "\nTakes a file with enrichment terms for each gene, and a list of genes\n";
	print "and calculates term enrichment (whatever the enrichment terms are).\n";
	print "\nusage:\n $0\n";
	print "-a [allow missing genes (absent from database) to contribute to enrichment score [1: true, 0: false (default)]]\n";
	print "-b [background list of genes (optional)]\n";
  print "-g [list of genes to test for enrichment]\n";
  print "-o [out file prefix, default 'enrichment']\n";
	print "-t [list of genes with enrichment terms]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-a'){
			$i++;
			$allowMissing = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-o'){
			$i++;
			$outputPrefix = $ARGV[$i];
		} elsif($ARGV[$i] eq '-t'){
			$i++;
			$tfile = $ARGV[$i];
		} else{
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

sub joinGenes {
  my @genes = @{$_[0]};
  my %idMap = %{$_[1]};
  my @withID;
  foreach my $gene (@genes) {
    push @withID, "$gene ($idMap{$gene})";
  }
  return join ', ', @withID;
}

my %tsvParams = (
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
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

# Read list of genes.
print STDERR "Reading gene file: $gfile\n";
open my $geneFH, '<', $gfile or die "Could not open $gfile: $!";
my $geneTSV = Text::CSV_XS->new(\%tsvParams);
my @geneList;
my %geneExists; # This hash is for checking for duplicate entries. Additional entries are ignored.
while(my $row = $geneTSV->getline($geneFH)) {
  my $gene = @{$row}[1];
  if (!(exists $geneExists{$gene})) {
    push @geneList, @{$row}[1];
    $geneExists{$gene} = 1;
  }
}
close $geneFH;

# Read domain information.
print STDERR "Reading term information: $tfile \n";
open my $domainFH, '<', $tfile or die "Could not open $tfile: $!";
my $termTSV = Text::CSV_XS->new(\%tsvParams);
$termTSV->getline($domainFH); # discard header
my %domainsPerGene;
my %idMap;
while(my $row = $termTSV->getline($domainFH)) {
  my $gene = @{$row}[0];
  my $id = @{$row}[1];
  my @domains = split ';', @{$row}[2];
  if (isInBackground($id, $bfile, \%backgroundGenes)) {
	  @{$domainsPerGene{$id}} = @domains;
    $idMap{$id} = $gene;
  }
}
close $domainFH;

# Find all genes per domain.
my %genesPerDomain;
for my $id (keys %domainsPerGene) {
	for my $domain (@{$domainsPerGene{$id}}) {
		push @{$genesPerDomain{$domain}}, $id;
	}
}

# Remove genes with no term information.
if (!$allowMissing) {
  for (my $i = scalar @geneList - 1; $i >= 0; $i--) {
    my $gene = $geneList[$i];
    if (!(exists $domainsPerGene{$gene})) {
      splice @geneList, $i, 1;
	  }
  }
}

# Count genes for each term.
my %domainCount;
my %genesWithTerm;
for my $gene (@geneList) {
  for my $geneDomain (@{$domainsPerGene{$gene}}) {
		if (exists $domainCount{$geneDomain}) {
			$domainCount{$geneDomain}++;
		} else {
			$domainCount{$geneDomain} = 1;
		}
		push @{$genesWithTerm{$geneDomain}}, $gene;
	}
}

# Calculate enrichment.
my $backgroundSize = scalar keys %domainsPerGene;
my %correctionArray;
my %queryOutput;
foreach my $term (keys %domainCount) {
  $queryOutput{$term}{'genes'} = joinGenes(\@{$genesWithTerm{$term}}, \%idMap);
  $queryOutput{$term}{'matched'} = $domainCount{$term};
  $queryOutput{$term}{'backgroundWithDomain'} = scalar @{$genesPerDomain{$term}};
  $queryOutput{$term}{'querySize'} = scalar @geneList;
  $queryOutput{$term}{'fEnrichment'} = sprintf '%.3f', ($queryOutput{$term}{'matched'} / $queryOutput{$term}{'querySize'}) / ($queryOutput{$term}{'backgroundWithDomain'} / $backgroundSize);
  $queryOutput{$term}{'pValue'} = fishers($queryOutput{$term}{'matched'}, $queryOutput{$term}{'querySize'}, $queryOutput{$term}{'backgroundWithDomain'}, $backgroundSize);
  $correctionArray{$term} = $queryOutput{$term}{'pValue'};
}
  
# Perform correction.
my ($adjustedPValueRef, $correctedFDRRef) = bhCorrection(\%correctionArray, $fdr);
my %adjustedPValue = %{$adjustedPValueRef};
my %correctedFDR = %{$correctedFDRRef};
my @sortedOutput;
my $sortedOrder = 0;
foreach my $term (sort { $correctedFDR{$a} <=> $correctedFDR{$b} } keys %correctedFDR) {
	$sortedOutput[$sortedOrder++] = $term;
}

# Output.
open my $outputFH, '>', 'output/' . $outputPrefix . '.txt';
print $outputFH "term\tmatched\tbackground_size\tfold enrichment\tpvalue\tadj. pvalue\tbhfdr\tgenes\n";
for(my $i = 0; $i < scalar @sortedOutput; $i++) {
	if ($queryOutput{$sortedOutput[$i]}{'pValue'} == 0 ||
		$queryOutput{$sortedOutput[$i]}{'pValue'} < $correctedFDR{$sortedOutput[$i]}
	) {
	  print $outputFH "$sortedOutput[$i]\t";
	  print $outputFH "$queryOutput{$sortedOutput[$i]}{'matched'}\t";
		print $outputFH "$queryOutput{$sortedOutput[$i]}{'backgroundWithDomain'}\t";
		print $outputFH "$queryOutput{$sortedOutput[$i]}{'fEnrichment'}\t";
		print $outputFH "$queryOutput{$sortedOutput[$i]}{'pValue'}\t";
		print $outputFH "$adjustedPValue{$sortedOutput[$i]}\t";
		print $outputFH "$correctedFDR{$sortedOutput[$i]}\t";
		print $outputFH "$queryOutput{$sortedOutput[$i]}{'genes'}\n";
	}
}
close $outputFH;
