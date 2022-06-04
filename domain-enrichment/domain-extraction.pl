#!/usr/bin/perl

# 13/6/2017

use strict;
use warnings;

# modules
use List::MoreUtils qw(uniq);
use String::Util qw(trim);
use Text::CSV_XS;

# parameters
my $species = 'Homo sapiens';

# command line parameters
my $pfile = ''; # Pfam file
my $ufile = '';	# Uniprot database

if ($#ARGV == 0) {
	print "\nTakes a Uniprot file and from it grabs gene-to-uniprot map. Then using Pfam\n";
	print "and the UniProt ID, grabs domains.\n\n";
	print "\nusage:\n $0\n";
	print "-p [Pfam file]\n";
  print "-s [species: default Homo sapiens]\n";
	print "-u [UniProt/SwissProt file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV) {
		if($ARGV[$i] eq '-p'){
			$i++;
			$pfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-s'){
			$i++;
			$species = $ARGV[$i];
		} elsif($ARGV[$i] eq '-u'){
			$i++;
			$ufile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Read Pfam file
print STDERR "Reading Pfam information from $pfile \n";
open my $pfamFH, '<', $pfile or die "Could not open $pfile: $!";
my $pfamTSV = Text::CSV_XS->new({ sep_char => "\t" });
## Discard first 3 lines
for (my $i = 0; $i < 3; $i++) {
  $pfamTSV->getline($pfamFH);
}
my %pfam;
while(my $row = $pfamTSV->getline($pfamFH)) {
	push @{$pfam{@{$row}[0]}}, @{$row}[6];
}
close $pfamFH;

# Get details for each reviewed uniprot entry matching species.
print STDERR "Parsing UniProt data from $ufile\n";
my $currGene = '';
my $currUniprot = '';
my %geneDetails;
my $reviewedBoolean = 0;
my $speciesBoolean = 0;
open(UFILE, "<$ufile") || die "$ufile can't be opened: $!";
while(<UFILE>){
	if($_ =~ /^ID\s+\S+\s+Reviewed;/) {
    $reviewedBoolean = 1;
  } elsif($reviewedBoolean && !$speciesBoolean && $_ =~ /^OS   $species/) {
    $speciesBoolean = 1;
  } elsif($reviewedBoolean && !$currGene && $_ =~ /^GN   Name=([^;{]+); Synonyms=([^;]+);/) {
    $currGene = trim($1);
  } elsif($reviewedBoolean && !$currGene && $_ =~ /^GN   Name=([^;{]+)/) {
    $currGene = trim($1);
  } elsif($reviewedBoolean && !$currGene && $_ =~ /^GN   OrderedLocusNames=([^;{]+)/) {
    $currGene = trim($1);
  } elsif($reviewedBoolean && !$currGene && $_ =~ /^GN   ORFNames=([^;{]+)/) {
    $currGene = trim($1);
  } elsif($reviewedBoolean && !$currUniprot && $_ =~ /^AC\s+([A-Z0-9]+);/) {
		$currUniprot = $1;
	} elsif($reviewedBoolean && $_ =~ /^\/\//) {
		if ($speciesBoolean) {
			if (exists $pfam{$currUniprot}) {
				@{$geneDetails{$currUniprot}{'domains'}} = uniq @{$pfam{$currUniprot}};
			} else {
				@{$geneDetails{$currUniprot}{'domains'}} = ();
			}
			$geneDetails{$currUniprot}{'gene'} = $currGene;
		}
		# Reset some parameters.
		$currGene = '';
		$currUniprot = '';
		$reviewedBoolean = 0;
		$speciesBoolean = 0;
	}
}
close(UFILE);

# Print list of domains per gene.
open my $detailsFH, '>', 'output/pfam-domains.txt';
print $detailsFH "Gene\tUniProt\tdomains\n";
for my $id (keys %geneDetails) {
	print $detailsFH "$geneDetails{$id}{'gene'}\t";
	print $detailsFH "$id\t";
	my $domainString = join ";", @{$geneDetails{$id}{'domains'}};
	print $detailsFH "$domainString\n";
}
close $detailsFH;
