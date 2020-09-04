# HLAXPress
RepAn is an R package that implements a differential abundance analysis method for deep sequenced TCR immune repertoire (Repseq) datasets to identify enriched/expanded clonotypes associated with a condition/disease.

## Requirements
To install RepAn package in R, first install R devtools package and run its install_github function as follows: 

```
install.packages("devtools")
devtools::install_github("dyohanne/RepAn")
```

```
## Usage
RepAn is mainly designed to work with genomic TCR repertoire datasets and accepts [immunoseq](https://www.adaptivebiotech.com/immunoseq) formatted datasets. It also accepts MiXCR format data, in which case, if the data comes from cDNA libraries, the analysis can be interpreted as differential expression of TCRs.

To perform differential abundance analysis, first read in the repertoire samples using the readSample function (for immunoseq format data) or the readMiXCR function (for MiXCR format data). For example to read condition 1 samples named sample1.tsv, sample2.tsv,... , and condition 2 samples named sample1T.tsv, ..., do the following ; 

```
library(RepAn)

s1=readSample("sample1.tsv")
s2=readSample("sample2.tsv")
s3=readSample("sample3.tsv")

s1T=readSample("sample1T.tsv")
s2T=readSample("sample2T.tsv")
s3T=readSample("sample3T.tsv")
```


