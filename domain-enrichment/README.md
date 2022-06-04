# Domain/motif enrichment

Determine what domains/motifs/etc are enriched in a input list relative to a background. The workflow below is for domains, but any terms can be used as long as the input file has the correct format.

## Download PFAM domains

First parse UniProt file for reviewed IDs and grab all domains for those entries.

### Input

* -p: Pfam domain file for species of interest (for humans: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/9606.tsv.gz.) Make sure to unzip.
* -s: species, default Homo sapiens
* -u: Reviewed UniProt entries from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz (make sure to unzip)

Run as `./domain-extraction.pl -p 9606.tsv -u uniprot_sprot.dat`

### Output

* `pfam-domains.txt`: gene name, UniProt ID and semi-colon separated list of domains

## Calculate enrichment

### Input

* -a: allow missing genes (absent from database) to contribute to enrichment score [1: true, 0: false (default)]
* -b (optional): background list of genes. If missing, uses all entries in term file (-t) as background.
```
geneA P000001
geneB Q000002
```
* -g: list of genes to test for enrichment
```
geneA P000001
geneB Q000002
```
* -t: list of genes with associated terms to test for enrichment
```
Gene	UniProt	domains
KLF15	Q9UIH9	zf-C2H2
TSACC	Q96A04	SSTK-IP
SNHG28	P0DPA3	
SLC38A9	Q8NBW4	Aa_trans
TRIM44	Q96DX7	zf-B_box
```

Run as `./enrichment.pl -g gene-list.txt -t pfam-domains.txt`

### Output

* `enrichment.txt`: enriched domains at a 1% FDR with statistics
