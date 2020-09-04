#!/usr/bin/env Rscript

# AUTHORS:
# Dawit A. Yohannes
# Immunomics Group Biomedicum, UH, Finland (https://www.helsinki.fi/en/researchgroups/immunomics)

# GENEREAL DESCRIPTION:
# This script estimates allele specific expression of HLA gene alleles for a sample from high-throughput sequencing data prepared with a UMI protocol.
# It requires the HLA allele types and can handle Illumina or nanopore read files.

# FEATURES:
# Makes all allele UMI counting at 4 digit HLA allele level instead of 6 digit HLA allele level
# Uses personalized HLA reference for allele-specific expression estimation

#________________________________________ Include all HLAXPress functions ___________________________________________________####
source("HLAXPress_fxns.R")

cat(".............................................................................................","\n")

#________________________________________ Parse-in arguments___________________________________________________####

#Arguments :

#readDir: 
  # The directory that contains Illumina or nanopore reads for samples. For both types of reads, HLAXPress expects the read files to start with the sample name followed by an underscore (sample12_).
  # For illumina reads, tool expects paired end R1 and R2 reads. For nanopore reads, all three types of reads per sample (2D, complement and template) should be combined into one .fastq or .fastq.gz reads file (e.g BC18_combined.fastq), or three separate files
  #	for 2D, template and complement reads can be provided which will be combined in HLAXPress into one file, or only one type of reads can be present in the reads directory for the sample, for example for sample BC18 (e.g BC18_2D.fastq).
  # Each sample must have a unique sampleName to identify it properly. The sample name must be present at the start of its reads file (followed by underscore), and must also be present in the first column hlaTypesFile (where sample's HLA allele types are provided).

#hlaTypesFile: 
  #A tab/space-delimited file containing the HLA allele types for each sample (in each row) to be analyzed. It contains colnames: sampleName, HLA-A, HLA-A, HLA-B,HLA-B,... (each HLA gene is entered twice with two allele entries)
  # First column (sampleName) should contain the sample name for each sample exactly as it appears in the beginning of its sequence read file in readDir for that sample.

#samName: 
  #provide sample name if analysis is to be done for a single sample, otherwise leave empty.

#hlaRef: 
  # Address of the HLA reference database file (the name after indexing with last aligner") (e.g hladb/totalHLADb). All other files produced when indexing with last aligner must be present in the same directory as this file. If not provided, the latest version of the HLA reference database from IMGT/HLA(ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta)
  # will be downloaded to the current directory, indexed with last, and used.

#readType: 
  #provide the read type, one of either nanopore or illumina read files

#umiLen: 
  #length of the umi (default 10bp), for now only umi length is required, TSO+primer sequences are hardcoded and may need to be modified in the script by user


# Create a parser
p <- arg_parser("This script estimates allele specific expression using UMIs after HLA-alleles are determined",name="HLAXPress")

# Add command line arguments
p <- add_argument(p, "readDir", help="Illumina reads or nanopore run directory")
p <- add_argument(p, "hlaTypesFile", help="HLA types file - space-delimited file having samplename HLA-A HLA-A HLA-B ...allele type records for each sample (each sample on one row)")

p <- add_argument(p, "--hlaRef", help="Address of the HLA reference database file (the name after indexing with last aligner)",default=NULL)
p <- add_argument(p, "--samName", help="Provide sample name if processing only for single sample is needed, otherwise all samples in hlaTypesFile get processed",default=NULL)
p <- add_argument(p, "--readType", help="Read type: either nanopore or illumina",default="illumina")
p <- add_argument(p, "--umiLen", help="length of the UMI",default=10)


# Parse the command line arguments
argv <- parse_args(p)

#________________________________________ Prepare to run analysis ____________________________________________ ####

#get sample reads 
entereddir <- argv$readDir

#get UMI length
umiLen <- argv$umiLen

#get read type
readType <- tolower(argv$readType)

#get HLA reference address
if(is.na(argv$hlaRef)){
  referenceDataFile <- NULL
}else{
  referenceDataFile <- argv$hlaRef
}

#get sample Name
if(is.na(argv$samName)){
  sid <- NULL
}else{
  sid <- argv$samName
  
}


#get sample HLA type information
fileInput = argv$hlaTypesFile

# HLA genes of interest for which the anlaysis is going to be done, this could by modified if needed
genesOfInterest = c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")
geneNames = as.character(sapply(genesOfInterest,function(x) return(substr(x,5,nchar(x)))))

#________________________________________ Run Allele-specific expression Analysis ____________________________________________ ####


if(readType=="illumina"){
  getAlleleFreqCounts(sampleHlaData=fileInput,sampleId=sid,IlluminaDir=entereddir,umiLength=umiLen,
                      hladbFile=referenceDataFile)
}else if(readType=="nanopore"){
  getAlleleFreqCounts(sampleHlaData=fileInput,sampleId=sid,nanDir=entereddir,umiLength=umiLen,
                      hladbFile=referenceDataFile)
}else{
  stop("No proper read type specified: either nanopore or illumina read types can be processed.")
}



