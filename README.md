# HLAXPress
HLAXPress performs allele-specific expression estimation for HLA genes using unique molecular identifiers (UMIs) to quantify HLA transcripts in an unbiased manner. It can do this from both whole transcriptome data (RNAseq) or targetted HLA amplicon data. It typically expects illumina paired end data. It was also used on nanopore ONT targetted HLA RNAseq sequence data that utilized UMIs, thus it can also handles nanopore data with minor modifications depending on your library design. 


## Requirements
HLAXPress runs in a linux environment with the following installed:
- R (>= 2.15.1)
- LAST (for alignment), LAST can be downloaded and installed from (http://last.cbrc.jp/). After installing LAST, add its path to environment PATH variable:  

```
export PATH=installation_Directory/last/bin:$PATH
```

## Usage
Currently, HLAXPress can be run either from the command-line or from within R after loading all the necessary functions. To run from the command line, use R's Rscript by passing the required arguments: readDir, the directory where to find the sequence reads (in fastq or fastq.gz, file names MUST start with the same sampleName as in the HLA type information for that sample in hlaTypesFile followed by underscore), and hlaTypesFile, a file containing the HLA type information of each sample (one sample per row, first column is sampleName, next columns contain HLA alleles for the sample). Optional arguments like address to the last indexed HLA reference database and read type (illumina or nanopore) can also be passed. To do a test run, assuming the read data and HLA type information are present in the testData directory:  

```
usage: HLAXPress.R [--help] [--hlaRef HLAREF] [--samName SAMNAME] [--readType READTYPE] [--umiLen UMILEN] readDir hlaTypesFile

# To see the help:
Rscript HLAXPress.R -h

# For illumina read data. It performs the analysis for two samples present in the HLA types file testData/testDataHLAs (assuming their fastq files are also present)
# This downloads the reference HLA database from IMGT/HLA and indexes it for use with LAST
Rscript HLAXPress.R testData testData/testDataHLAs.txt --hlaRef totalHLAdb/hladb --readType illumina

# for nanopore data (assuming three separate fastq.gz files in testData for template, complement and 2D reads). Analysis is done for single sample named BC01
# We can use the already downloaded and indexed HLA reference data with --hlaRef
Rscript HLAXPress.R testData testData/testDataNanoporeDataHLAs.txt --samName BC01 --hlaRef totalHLAdb/hladb --readType nanopore

```
To run an analysis from within R:

```
# first load all the functions
source("HLAXPress_fxns.R")

# Run the analysis for sample45 (which is illumina read data)
getAlleleFreqCounts(sampleHlaData="testData/testDataHLAs.txt",sampleId="sample45",IlluminaDir="testData/testData")
```
For detailed description of the arguments required, read the description of the getAlleleFreqCounts function in HLAXPress_fxns.R and the command-line arguments in HLAXPress.R. 


## Brief description
Although the HLA region is highly variable, distinguishing between allele types from sequence reads aligning to the HLA region is a difficult task. The HLA alleles can differ by as little as one or few polymorphisms, as such, HLA alleles also share large parts of their sequence. For perfect identification of allele types, sequence reads with no sequencing errors are needed, which is not possible with the current technology. Thus, as reads originating from HLA transcripts typically multi-map to multiple reference alleles, careful discrimination between them is critical when performing allele-specific expression quantification. 

To achieve this, HLAXpress first builds a personalized HLA reference database per sample (the HLA types need to be provided, for now HLA typing should be done by other tools, future HLAXPress versions will come with HLA typing functionality), it then collects reads that align to the personalized HLA reference. It then proceeds to evaluate each UMI separately (all reads having that UMI -> originating from same transcript). UMIs that have reads aligning to only one specific allele are counted to that allele (just once), UMIs that have multi-mapping reads are assigned to the allele that gives the highest likelihood estimate P(UMI|Allele). The allele likelihoods are calculated based on how well UMIs resemble each allele at the key polymorphic sites, which are defined as the sites with shannon's entropy greater than 0.50 after multiple alignment of the alleles in the personalized HLA reference (The key polymorphic sites are high diversity positions that derive the difference between the alleles, the cutoff allows to pick highly variable sites with > ~10% alternative nucleotide types).


