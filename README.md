# Bam Read Index (bti)

`bti` is a small tool to quickly retrieve alignments in a bam file by tag and is a minor modification and fork of the `bri` repo. With 10x genomics and other single cell data, it is sometimes desirable to extract all of the alignments for a particular cell or molecular barcode. `bti` is designed for this usecase and provides both a command line interface and htslib-inspired API. The index for ~15Gb of 10x Genomics single cell RNA-Seq data  with ~ 400K barcodes is under 20 MB, and querying takes ~ 1 sec per cellbarcode.

## Installation

```
git clone https://github.com/jts/bri.git
cd bri
make HTSDIR=/path/to/htslib
```

## Usage

First sort the bam file by the tag of interest (only string type "Z" tags are implemented fo now) using `samtools sort -t CB`

Then index the bam file by tag (hardcoded as CB for now):

```
> bri index test/cb_tag_sorted_small.bam
```

Extract the alignments for a particular tag value (output is in BAM):

```
> bri get test/cb_tag_sorted_small.bam AAACCTGGTACGCTGC-1 | samtools view | head -n 2

K00282:115:HJTMNBBXX:5:1223:0:2858381	0	1	88217701	255	98M	*	0	0AATGACATCCATGCTCGATCACAACCTGCTCTCAAAGTGACACGTTCCAGCTCTGGGGACAAAGGCCGCCATCTCTCTCAATGGCTTCTCACCTCCCA	AAFAFJJJJJJJJJJFAJJJFFAAJFJFJJJJJ<7AAFJJJFFJAFJA<-FJJFJJFF-A<-<AFJFFJFAAAFJJFJJFAJJFJJJJJJ7FFJJJJF	NH:i:1	HI:i:1	AS:i:94	nM:i:1	NM:i:1	CR:Z:AAACCTGGTACGCTGC	CY:Z:AAFFFJJJJJJJJJJJ	CB:Z:AAACCTGGTACGCTGC-1	UR:Z:AAAGCACAGA	UY:Z:JJJJJJJJJJ	UB:Z:AAAGCACAGA	BC:Z:GAAGGCTG	QT:Z:AAFFF-AJ	
K00282:115:HJTMNBBXX:5:2108:0:3201626	0	1	165292120	255	27S71M	*	0	0AAGCAGTGGTATCAACGCAGAGTACATGGGCGCCTTTAATCCCAGCACTCGGGAGGCAGAGGCAGGCGGATTTCTGAGTTAGAGGCCAGCCTGGTCTA	AA<AFFJFFJ7FJFJJJJJJAAJAJ7FJJJJJF77FF7<<<-7<7-77<7FJJ7JJ<7F7AJF-AFAJJ77AF7AJ-7FF7F7AJ7A7A7FFJJAAAA	NH:i:1	HI:i:1	AS:i:69	nM:i:0	NM:i:0	CR:Z:AAACCTGGTACGCTGC	CY:Z:AAFF-FFJJJFJFJJA	CB:Z:AAACCTGGTACGCTGC-1	UR:Z:TCTTGAGGGG	UY:Z:JFJJJFJJJJ	UB:Z:TCTTGAGGGG	BC:Z:TGGATTGC	QT:Z:AAFAFJJA	
```


