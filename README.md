# Bam Read Index (bri)

`bri` is a small tool to quickly retrieve alignments in a bam file by readname or tag. Most bam files are sorted by genomic position to allow fast sequential reading across regions of the genome. Bam does not allow random access to alignments by readname or tag as short reads would require very large readname-to-file position indices. With long read sequencing however, it is sometimes desirable to extract all of the alignments for a particular read. `bri` is designed for this usecase and provides both a command line interface and htslib-inspired API. The index for ~50X of nanopore data for a human genome is under 1GB.

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
> bri index  reads.sorted.bam
```

Extract the alignments for a particular tag value (output is in BAM):

```
> bri get reads.sorted.bam AAATGCCGTCCGTTAA-1  | samtools view | head -n 2

K00282:115:HJTMNBBXX:6:1118:0:1615302	0	1	5126662	255	98M	*	0	0	TAAGGCTGCTATGTACATAGTGGAGCATGTGTCCTTGTTATATGTTGGAGAATGACACATCAAAACGACCATTCACCATGATCAATTAGGCTTCATCC	AAA-FFJ<JFJJFJJFJFFJJJJA<JJJJJJJJJJJJJJ<FFFJJJJJFF7FFFFJFJFFJFJJJJJAJJFAAJFJFJFJ<FJFAFJJFJJFJJFJJJ	NH:i:1	HI:i:1	AS:i:96	nM:i:0	NM:i:0	CR:Z:AAATGCCGTCCGTTAA	CY:Z:AAFFFJJJJJJJJJJJ	CB:Z:AAATGCCGTCCGTTAA-1	UR:Z:AGCAAGGCGC	UY:Z:JJJJJJJJJJ	UB:Z:AGCAAGGCGC	BC:Z:GAAGGCTG	QT:Z:AAFFF-AF	
K00282:115:HJTMNBBXX:6:2101:0:3559433	0	1	5126705	255	98M	*	0	0	GTTGGAGAATGACACATCAAAACGACCATTCACCATGATCAATTAGGCTTCATCCCAGGGATACAGGGATGGTTCAATATATGGAAATCCATCAATGT	<77AFFFJJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJFJAFJJJJJFFAJJJJJJJFJJJ	NH:i:1	HI:i:1	AS:i:96	nM:i:0	NM:i:0	CR:Z:AAATGCCGTCCGTTAA	CY:Z:AAFFFJJJJJJJJJJJ	CB:Z:AAATGCCGTCCGTTAA-1	UR:Z:AGCAAGGCGC	UY:Z:JJJJJJJJJJ	UB:Z:AGCAAGGCGC	BC:Z:ATTCCGAT	QT:Z:AAFFFJFJ
```


