# for illumina read data
Rscript HLAXPress.R testData testData/testDataHLAs --hlaRef totalHLAdb/hladb --readType illumina

# for nanopore data (three separate fastq.gz files in testData ), this is not yet tested, it had some error, check it
# hla allele info of sample47 from illumina used for BC01, for quick testing
Rscript HLAXPress.R testData testData/testDataHLAs --samName BC01 --hlaRef totalHLAdb/hladb --readType nanopore
