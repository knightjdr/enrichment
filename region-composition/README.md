# Region composition

Calculate percentage of a protein that has a specific region/motif (e.g. disorder). First used in Youn et al (2018), *High-Density Proximity Mapping Reveals the Subcellular Organization of mRNA-Associated Granules and Bodies*, see figure 6D.

## Download PFAM motifs

Fetch motifs (e.g. disorder) from PFAM. Motifs from PFAM can only be retrieved programmatically one-by-one. For example, information for the UniProt ID Q9NRA8 can be retrieved in JSON format at this URL: http://pfam.xfam.org/protein/Q9NRA8/graphic.

### Input

* -u: Reviewed UniProt entries from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz (make sure to unzip)

Run as `./pfam-fetch.pl -u uniprot_sprot.dat`

### Output

* `pfam-motifs.txt`: list motifs for each gene name/UniProt ID with no duplicate motifs per gene
* `pfam-motifs-complete.txt`: all motifs (including duplicates), along with start and end positions of the motif and length

## Fractional summary

This will produce stats on a motif in question for both the total proteome and a list of query genes. It will determine the mean and median fractional protein length for the region of interest (ROI) (a protein may have multiple instances of the region, so this is the fraction of the protein that has it, summing up the lengths of the separate regions), the mean and median longest region (for each protein the longest region is retained) and the mean and median total amino acid count with that motif (counts all amino acids comprising the motif of interest per protein).

### Input
* -b (optional): Background list with gene names and UniProt IDs in set to test for enrichment. If missing, uses all entries in Pfam file as background.
```
geneA P000001
geneB Q000002
```
* -l: List with gene names and UniProt IDs in set to test for enrichment
```
geneA P000001
geneB Q000002
```
* -m: motif of interest (default: disorder)
* -p: Pfam file with motif information for each protein
```
gene	uniprot	motif	start	end	motif length	protein length
ADAM23	O75077	disorder	1	21	21	832
ADAM23	O75077	disorder	27	37	11	832
ADAM23	O75077	disorder	64	65	2	832
```

Run as `./fractional-summary.pl -b background.txt -l genes.txt -m disorder -p pfam-motifs-complete.txt`

### Output

* `histograms.txt`:
  * bins for fractional protein length that is the region of interest
  * bins go from 0-5% to 95-100%
  * first column is for background and second for genes of interest
* `<region>-summary.txt`:
  * &lt;region&gt; will be replaced by the region/motif name being summarized
  * summary includes mean and median for background and genes of interest
* `wilcox_*_*.txt`:
  * list of values for all background proteins or genes of interest for performing Wilcoxian test

## Violin and box plots

The values in the `wilcox_*_*.txt` files can be used to make box or violin plots with `boxplot.R`.

For example, to create plots for the longest region, create a text file as below with the columns containing the values from the `wilcox` files

```
Background  List1  List2
0.02        0.5    0.17
0.57        0.34   0.18
...         ...     ...
```

Specify the name of the file in the script and the column headings to use.

### Output

* `longest_boxplot.pdf`: boxplot
* `longest_boxplot_noOutliers.pdf`: boxplot not showing outlier points
* `longest_violin.pdf`: violin plot

## Test for significance

To test for significance between the fractional region length of a query list versus the background (or the longest region instead of fraction), the `wilcox_*.*.txt` files can be used for performing a Wilcox test using the script `wilcox_test.R`

To run, open R and specify file names to import and compare.
