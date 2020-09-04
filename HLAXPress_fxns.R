# AUTHORS:
# Dawit A. Yohannes
# Immunomics Group Biomedicum, UH, Finland (https://www.helsinki.fi/en/researchgroups/immunomics)

# GENEREAL DESCRIPTION:
# This script estimates allele specific expression of HLA gene alleles for a sample from high-throughput sequencing data prepared with a UMI protocol.
# It requires the HLA allele types and can handle Illumina or nanopore read files.


#____________________________________ Load R packages needed, install them if they are not already installed___####

loadPacks <- function(package.list = c("ggplot2", "argparser","cluster","iterators","parallel","doParallel")){
  new.packages <-package.list[!(package.list %in% installed.packages()[,"Package"])]
  if(length(new.packages) > 0) install.packages(new.packages,repos='http://cran.us.r-project.org')
  loadSuccess <- lapply(eval(package.list), library, character.only=TRUE,quietly=TRUE)
}


loadBioconductorPacks <- function(package.list = c("Biostrings","ShortRead","Rsamtools","msa")){
  new.packages <-package.list[!(package.list %in% installed.packages()[,"Package"])]
  if(length(new.packages) > 0){
    if(R.Version()[["major"]] >= 3 & unlist(strsplit(R.Version()[["minor"]],"\\."))[1] >= 5){
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(new.packages)
      
    }else{
      source("http://bioconductor.org/biocLite.R")
      biocLite(new.packages)
    }
  }
  loadSuccess <- lapply(eval(package.list), library, character.only=TRUE,quietly=TRUE)

}


suppressMessages(loadPacks())
suppressMessages(loadBioconductorPacks())








#________________________________________ Analysis Functions ____________________________________________ ####

lastExists <- function(){
  
  if(Sys.which("lastal") != "")
    return(T)
  else
    return(F)
  
}


performAlignmentToTotalHLAData <- function(R1,R2,nR,referenceDataFile){
  
  # FIRST: perform alignment of all reads to the complete HLA reference database
  # The address to the total HLA reference file (the last alignment indexing should already 
  #be done and be present in the directory containing the reference file)
  
  # align reads to total reference database
  
  if(!lastExists()){
    stop("last aligner is not found in your system !!")
  }
  
  wd <- getwd()
  dbdir <- dirname(referenceDataFile)
  #setwd(dbdir)
  
  if(is.null(nR)){
    if(file.exists(R1) & file.exists(R2)){
      
      # -T 0 for local alignment--normally used
      # -T 1 for overal alignment--end to end---checking it later
      # -l 50, we use 50 bases to start the initial alignment seed lenght because the 300bp R2 quality drops very fast after that point. R1 reads are basically useless for this because they don't reach exon 2 and 3 
      # we tested both -l 50 and -l 100, in both cases small number of reads aligned to both alleles, because of the paired-end alignment in last
      
      # the paired end alignment in last works best (https://www.ncbi.nlm.nih.gov/pubmed/23413433), So we use it.
      
      system(paste("lastal -s 2 -T 0 -l 50 -m 1 -a 100 -Q 1 -i1",referenceDataFile,R2,paste("> ",dbdir,"/tempAlmt_R2",sep=""),sep=" "))
      system(paste("lastal -s 2 -T 0 -l 50 -m 1 -a 100 -Q 1 -i1",referenceDataFile,R1,paste("> ",dbdir,"/tempAlmt_R1",sep=""),sep=" "))
      system(paste("last-pair-probs",paste(dbdir,"/tempAlmt_R1",sep=""),paste(dbdir,"/tempAlmt_R2",sep=""),paste("> ",dbdir,"/tempAlmt",sep=""),sep=" ")) # merge the alignment of R1 and R2 into one file
      
      # Converting the maf file to sam commented out for now, not needed
      #try(system(paste("maf-convert sam","tempAlmt","> tempAlmt.sam",sep=" ")))
      
      
    }else{
      stop("Fastq read files not found for the sample! \n")
      
    }
    
  }else if(is.null(R1)){
    
    if(file.exists(nR)){
      
      # -T 0 for local alignment--normally used
      # -T 1 for overal alignment--end to end---checking it now
      
      system(paste("lastal -s 2 -T 0 -l 50 -m 1 -a 100 -Q 1",referenceDataFile,nR,paste("> ",dbdir,"/tempAlmt",sep=""),sep=" "))
      
      # Converting the maf file to sam commented out for now, not needed
      #try(system(paste("maf-convert sam","tempAlmt","> tempAlmt.sam",sep=" ")))
      
      
    }else{
      warning("Fastq file of nanopore reads was not found for the sample !\n")
      
    }
    
  }
  
  #setwd(wd)
  
  
  return(T)
  
}



prepareLastIndex <- function(databaseFile){
  
  if(!lastExists()){
    stop("last aligner is not found in your system !!")
  }
  
  wd <- getwd()
  dbdir <- dirname(databaseFile)
  dbname <- basename(databaseFile)
  setwd(dbdir)
  
  idxRes <- try(system(paste("lastdb -Q 0",dbname,dbname,sep=" ")),T)
  
  setwd(wd)
  
  if(class(idxRes) != "try-error"){
    return(T)
  }else{
    return(F)
  }
  
}


downloadHLAdb <- function(){
  
  # first create total hla db directory
  dir.create("totalHLAdb")
  
  
  imgthladbFile <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta"
  
  
  downloadedhladbFile <- download.file(imgthladbFile, paste("totalHLAdb/",basename(imgthladbFile),sep=""),quiet=T)
  
  if(downloadedhladbFile != 0){
    stop(paste("Error. HLA database file could not be downloaded from IMGT/HLA !",sep=""))
  }else{
    
    refDataFile <- paste("totalHLAdb/",basename(imgthladbFile),sep="")
    
    # read the downloaded fasta database file
    totalHLAref <- Biostrings::readDNAStringSet(refDataFile) 
    
    # edit the hla allele names. Remove the hla id numbers
    editedNames <- sapply(strsplit(names(totalHLAref)," "),function(x) x[2])
    names(totalHLAref) <- editedNames
    
    # write the hla database back
    refDataFileWritten <- paste("totalHLAdb/","hladb",sep="")
    Biostrings::writeXStringSet(totalHLAref, file=refDataFileWritten,format="fasta")
    
    #fileremoved <- file.remove(refDataFile,showWarnings = F)
    
    # after downloading the total hla db and indexing it with lastdb, we return the address of the total hla db
    lastindexing <- prepareLastIndex(refDataFileWritten)
    
    if(lastindexing){
      return(refDataFileWritten)
    }else{
      return(F)
    }
    
  }
  
  
}

# msf files not used
downloadHLAmsfFiles <- function(genes){
  
  # first create total hla db directory
  dir.create("msf")
  
  
  
  for(geneName in genes){
    hlaGeneMultipleAlignmentFile <- paste("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/",geneName,"_nuc.msf",sep="")
    
    downloadedMsfFile <- download.file(hlaGeneMultipleAlignmentFile, paste("msf/",geneName,"_nuc.msf",sep=""),quiet=T)
    
    if(downloadedMsfFile != 0){
      hlaGeneMultipleAlignmentFile <- paste("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/",geneName,"_gen.msf",sep="")
      downloadedMsfFile <- download.file(hlaGeneMultipleAlignmentFile, paste("msf/",geneName,"_nuc.msf",sep=""),quiet=T)
      
      if(downloadedMsfFile != 0)
        stop(paste("Error. Multiple alignment file for gene ",geneName," could not be downloaded !",sep=""))
    }
    
  }
  
  return(T)
  
}

downloadAlleleStatus <- function(){
  ## read allele_status information directly from imgt/hla and make selection using that information
  hdb <- read.table("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allele_status.txt",header=T,sep = ",",stringsAsFactors=F)
  
  return(hdb)
}


downloadAlleleFrequencies <- function(){
  ## read allele frequency information for European population from http://www.allelefrequencies.net/extaccess.asp
  # use rvest for scrapping the table from allelefrequencies webpage
  library("rvest")
  
  hlafreqDF <- NULL
  
  for(p in 1:30000){
    
    dataAddress <- paste("http://www.allelefrequencies.net/hla6006a.asp?page=",
                         p,"&hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=&hla_country=&hla_dataset=&hla_region=Europe&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=equal&hla_sample_year=&hla_sample_year_pattern=equal&hla_level=2&hla_level_pattern=equal&hla_show=&hla_order=order_2&standard=a",sep="")
    
    
    hlafreq <- dataAddress %>%
      read_html() %>%
      html_nodes(xpath='//*[@id="divGenDetail"]/table[1]') %>%
      html_table()
    
    if(nrow(hlafreq[[1]]) == 0)
      break
    
    hlafreqTemp <- hlafreq[[1]][,c(2,6)]
    hlafreqTemp[,1] <- sapply(strsplit(hlafreqTemp[,1]," "),function(x) x[length(x)])
    colnames(hlafreqTemp) <- c("Allele","AlleleFrequency")
    
    hlafreqDF <- rbind(hlafreqDF,hlafreqTemp)
    #print(paste("downloading page ",p,sep=""))
    
  }
  
  hlafreqDF$Allele <- sapply(hlafreqDF$Allele,function(x) paste(strsplit(x,"")[[1]][-length(strsplit(x,"")[[1]])],collapse=""))
  
  write.table(hlafreqDF,file="HLAlleleFrequencyEurope.txt",row.names=F,col.names=T)
  
  return(hlafreqDF)
}

prepareNanoporeData <- function(sampleName,nR,sampleResultDir){
  
  if(!grepl("\\.f.*q.*",nR)){
    
    nDir <- trimws(nR)
    
    #ontPath <- paste(nDir,"/downloads/pass/",sampleName,paste("/",sampleName,"_combinedPassFastq",sep=""),sep="")
    ontPath <- nDir
    #nRfiles <- list.files(path = ontPath, pattern=paste(".*",sampleName,".*","\\.fastq.gz",sep=""),ignore.case=TRUE,full.names=T,recursive=T)
    nRfiles <- list.files(path = ontPath, pattern=paste("^",sampleName,"_.*","\\.","(fastq.gz|fq.gz|fastq|fq)$",sep=""),ignore.case=TRUE,full.names=T,recursive=T)
    if(length(nRfiles)==0){
      ontPath <- paste(nDir,"/",sampleName,sep="")
      nRfiles <- list.files(path = ontPath, pattern=paste(".*",sampleName,"_.*","\\.","(fastq.gz|fq.gz|fastq|fq)$",sep=""),ignore.case=TRUE,full.names=T,recursive=T)
    }
    
    if(length(nRfiles)==0){
      stop(paste("Read files for sample ",sampleName," could not be found !\n",sep=""))  
    }
    
    # prepare read file address
    
    if(length(nRfiles)==1){
      nR_fileAddress <- nRfiles
    }else{
      #wd <- getwd()
      #resdir <- sampleResultDir
      #setwd(resdir)
      
      # combine all read files
      nR_fileAddress <- paste(sampleResultDir,"/",paste(sampleName,"tempComp2DReads.fastq",sep="_"),sep="")
      
      fcon <- file(nRfiles[1])
      if(summary(fcon)$class == "gzfile"){
        combineONTReads <- try(system(paste("zcat",nRfiles[1],nRfiles[2],nRfiles[3],">",nR_fileAddress,sep=" ")),T)
      }else{
        combineONTReads <- try(system(paste("cat",nRfiles[1],nRfiles[2],nRfiles[3],">",nR_fileAddress,sep=" ")),T)
      }
      close.connection(fcon)
      #setwd(wd)
      if(class(combineONTReads) == "try-error"){
        stop(paste("Unable to combine nanopore ONT reads for ",sampleName,"!",sep="")) 
      }
    }
    
    
    
    return(nR_fileAddress)
  }else{
    return(FALSE)
  }
  
}

getAllele4d <- function(allele){
  
  allele <- trimws(allele)
  
  if(grepl("^HLA",allele)){
    allele <- substr(allele,5,nchar(allele))
  }
  
  allele1_4digit <- ifelse(length(unlist(gregexpr(":",allele)))==1,
                           substr(allele,1,nchar(allele)),
                           substr(allele,1,unlist(gregexpr(":",allele))[2]-1))
  
  if(nchar(strsplit(allele1_4digit,"\\*")[[1]][2]) ==5)  
    allele <- allele1_4digit
  
  allele
}

findIlluminaReadAddress <- function(sampleName, IlluminaDir){
  
  IlluminaDir <- trimws(IlluminaDir)
  R1_read <- list.files(path = IlluminaDir, pattern=paste("^",sampleName,"_.*","R1",".*","\\.","(fastq.gz|fq.gz|fastq|fq)$",sep=""),ignore.case=TRUE,full.names=T,recursive=T)
  R2_read <- list.files(path = IlluminaDir, pattern=paste("^",sampleName,"_.*","R2",".*","\\.","(fastq.gz|fq.gz|fastq|fq)$",sep=""),ignore.case=TRUE,full.names=T,recursive=T)
  
  if(length(R1_read)==0 || length(R2_read)==0){
    stop(paste("Reads for sample ",sampleName," could not be found in the Reads directory!",sep=""))      
    
  }
  
  return(c(R1_read,R2_read))
}

getAlignedReadNamesAndReads <- function(tempAlignmentFile){
  
  allAlignments = readLines(tempAlignmentFile)
  if(length(allAlignments)==0){
    return(NULL)
  }
  
  headerEntries = grep("^#", allAlignments)
  
  allAlignments = allAlignments[-headerEntries]
  
  
  
  allAlignedRefAlleles = allAlignments[grep("^a", allAlignments) + 1]  #to get read names that aligned to either of the alleles
  
  allAlignedReads = allAlignments[grep("^a", allAlignments) + 2]  #to get read names that aligned to either of the alleles
  allAlignedReadQual = allAlignments[grep("^a", allAlignments) + 3]  #to get read names that aligned to either of the alleles
  
  allAlignedReadsSplit = strsplit(allAlignedReads," ")
  allAlignedReadQualSplit =  strsplit(allAlignedReadQual," ")
  allAlignedRefAllelesSplit = strsplit(allAlignedRefAlleles," ")
  
  allReadNames = sapply(allAlignedReadsSplit,function(x) return(x[2]))
  allReads =  sapply(allAlignedReadsSplit,function(x) return(x[length(x)]))
  allReadQ =  sapply(allAlignedReadQualSplit,function(x) return(x[length(x)]))
  
  allAlignedRefAllelesSplitStripped <- lapply(allAlignedRefAllelesSplit,function(k) k[k!=""])
  
  RefAlleles <- sapply(allAlignedRefAllelesSplitStripped,function(x) return(x[2]))
  startAlignRefPos = sapply(allAlignedRefAllelesSplitStripped,function(x) return(as.numeric(x[length(x)-4])+1)) # we add 1 because last uses 0-based coordinates
  endAlignRefPos = sapply(allAlignedRefAllelesSplitStripped,function(x) return(as.numeric(x[length(x)-4])+ 1 + as.numeric(x[length(x)-3]) )) # we add the spanned amount bases to the start position
  
  alignedR1Set <- DNAStringSet(allReads)
  names(alignedR1Set) <- allReadNames
  
  return(list(alignedR1Set,q=allReadQ,hitAlleles=RefAlleles,alStart=startAlignRefPos,alEnd=endAlignRefPos))
  
}


collectAllAlignedReads <- function(readType,writeToForSample){
  wd <- getwd()
  resdir <- writeToForSample
  setwd(resdir)
  
  if(readType=="illumina"){
    
    #get the reads aligning uniquely to the gene (instead of collecting all multimapping reads)
    allAligned_R1 <- getAlignedReadNamesAndReads("tempAlmt_R1")
    allAligned_R2 <- getAlignedReadNamesAndReads("tempAlmt_R2")
    allAlignedPE <- getAlignedReadNamesAndReads("tempAlmt")
    names(allAlignedPE[[1]]) <- sapply(strsplit(names(allAlignedPE[[1]]),"/"),function(x) return(x[1]))
    
    # Collect all reads and their qualities that aligned to the gene
    combinedR1R2AlignedReads <- c(allAligned_R1[[1]],allAligned_R2[[1]],allAlignedPE[[1]])
    combinedR1R2AlignedReads_Qual = c(allAligned_R1[[2]],allAligned_R2[[2]],allAlignedPE[[2]])
    combinedR1R2AlignedReads_Qual <- FastqQuality(BStringSet(combinedR1R2AlignedReads_Qual))
    
  }else{
    
    allAlignednR <- getAlignedReadNamesAndReads("tempAlmt")
    
    combinedR1R2AlignedReads <- c(allAlignednR[[1]])
    combinedR1R2AlignedReads_Qual = c(allAlignednR[[2]])
    combinedR1R2AlignedReads_Qual <- FastqQuality(BStringSet(combinedR1R2AlignedReads_Qual))
  }
  
  setwd(wd)
  
  return(list(combinedR1R2AlignedReads,combinedR1R2AlignedReads_Qual))
  
}

getPolymorphicSitesBetweenAlleles <- function(referenceSeqs=NULL){
  #get the variant positions in the alleles
  #invisible(capture.output(alleleAlmnt <- msa(referenceSeqs,order="input"))) #multiple align the reference sequences
  
  # for the multiple alignment below, error was occuring for some samples similar to this here (https://support.bioconductor.org/p/119291/)
  # some how the default clustalW alignment causes an error in some instances so we use clustalOmega
  
  #alleleAlmnt <- msa(referenceSeqs,order="input")
  
  invisible(capture.output(alleleAlmnt <- msa(referenceSeqs,method="ClustalOmega",order="input")))
  conMat <- consensusMatrix(alleleAlmnt)
  
  # Since the reference sequences have different length, the polymorhic positions that are found in both are selected
  # conMatPositionFreq = apply(conMat[!rownames(conMat) %in% c("-","+","."),],2,function(x) x/sum(x)) 
  
  # We just exclude conservation related symbols
  conMatPositionFreq = apply(conMat[!rownames(conMat) %in% c("+","."),],2,function(x) x/sum(x)) 
  
  posEntrop = apply(conMatPositionFreq,2,function(x) -sum(x[x > 0] * log2(x[x > 0]))) # entropy is calculated using elements in x that have frequency values more than 0 (the elements are frequencies of nucleotides)
  
  polymorphicPositions = which(posEntrop > 0.50)
  
  tempHLAreferencesAfterAlignment <- as.character(alleleAlmnt)
  
  return(list(tempHLAreferencesAfterAlignment,polymorphicPositions))
}

getReadQualityFromALignment <- function(alignedReads,alignedReadsQuality,readQualityEncoding=NULL,keyPositions=NULL){ 
  
  
  if(!is.null(keyPositions)){
    # get the corresponding quality at the polymorphic position for each read
    
    readOverlapStartPos <- apply(as.matrix(alignedReads),1,function(x) which(x!="-")[1]-1)
    
    allReadsAlmnt.atPositions_Quality = NULL
    for(rec in 1:nrow(alignedReads)){
      
      posinRead <- keyPositions - (readOverlapStartPos[rec])
      posinRead[posinRead <= 0 ] <- NA # remove zero and negative positions in read corresponding to early polymorphic positions which the read did not align to 
      
      
      BaseQualityAtPositions <- readQualityEncoding[strsplit(as.character(alignedReadsQuality[[rec]]),"")[[1]][posinRead]]
      
      errorProbAtBase = 10^(-BaseQualityAtPositions/10)
      
      allReadsAlmnt.atPositions_Quality <- rbind(allReadsAlmnt.atPositions_Quality,errorProbAtBase)
      
    }
    
    rownames(allReadsAlmnt.atPositions_Quality) <- rownames(alignedReads)
    return(allReadsAlmnt.atPositions_Quality) 
    
  }else{
    # get the corresponding quality for all sites in the aligning read
    
    allReadsAlmnt.atPositions_Quality = list()
    for(rec in 1:nrow(alignedReads)){
      
      alPos <- which(as.matrix(alignedReads)[rec,] != "-")
      
      posinRead <- alPos - (alPos[1] -1)
      posinRead[posinRead <= 0 ] <- NA # remove zero and negative positions in read corresponding to early polymorphic positions which the read did not align to 
      
      
      BaseQualityAtPositions <- readQualityEncoding[strsplit(as.character(alignedReadsQuality[[rec]]),"")[[1]][posinRead]]
      
      errorProbAtBase = 10^(-BaseQualityAtPositions/10)
      
      allReadsAlmnt.atPositions_Quality[[paste(rownames(alignedReads)[rec],rec,sep="_")]] <- errorProbAtBase
      
    }
    
    
    return(allReadsAlmnt.atPositions_Quality)
    
    
    
  }
  
}

getCoveredReferenceRegions <- function(alignedReads,refseq){
  readAligningReferenceRegions = list()
  for(rec in 1:nrow(alignedReads)){
    alPos <- which(as.matrix(alignedReads)[rec,] != "-")
    
    refSeqCoveredBases <- strsplit(as.character(refseq),"")[[1]][alPos]
    
    readAligningReferenceRegions[[paste(rownames(alignedReads)[rec],rec,sep="_")]] <- refSeqCoveredBases
    
  }
  
  return(readAligningReferenceRegions)
  
}

getAligningReadRegions <- function(alignedReads){
  
  readSegmentsAligningToRef = list()
  for(rec in 1:nrow(alignedReads)){
    alPos <- which(as.matrix(alignedReads)[rec,] != "-")
    
    refCoveringReadBases <- as.matrix(alignedReads)[rec,][alPos]
    
    readSegmentsAligningToRef[[paste(rownames(alignedReads)[rec],rec,sep="_")]] <- refCoveringReadBases
    
  }
  
  return(readSegmentsAligningToRef)
  
}

# Accuracy of UMI read in the aligned region
getUMIreadsAccuracy <- function(umiRead,referenceSequence){
  allRefBases <-  unlist(strsplit(as.character(referenceSequence),""))
  allReadBases <- as.character(umiRead)
  
  ReadBaseAvailIdx <- which(allReadBases != "-")
  
  if(length(ReadBaseAvailIdx)!=0){
    # match - mismatch
    alignedSegmentAccuracyToReference <- (sum(allReadBases[ReadBaseAvailIdx]==allRefBases[ReadBaseAvailIdx]) - sum(allReadBases[ReadBaseAvailIdx]!=allRefBases[ReadBaseAvailIdx])) / length(ReadBaseAvailIdx)
    
    alignedSegmentAccuracyToReference <- ifelse(alignedSegmentAccuracyToReference < 0,0,alignedSegmentAccuracyToReference)
    
  }else{
    alignedSegmentAccuracyToReference = 0
  }
  alignedSegmentAccuracyToReference
}

getUMIreadsKeySitesAccuracy <- function(umiRead,umiReadQuality,referenceSequence,polymorphicPositions){
  allRefBases <-  unlist(strsplit(as.character(referenceSequence),""))[polymorphicPositions]
  allReadBases <- as.character(umiRead)[polymorphicPositions]
  umiReadQuality <- as.numeric(umiReadQuality)
  
  ReadBaseAvailIdx <- which(allReadBases != "-")
  if(length(ReadBaseAvailIdx)!=0){
    
    matchAgg <- ifelse(length(which(allReadBases[ReadBaseAvailIdx]==allRefBases[ReadBaseAvailIdx])) > 0,
                       sum(1-umiReadQuality[which(allReadBases[ReadBaseAvailIdx]==allRefBases[ReadBaseAvailIdx])],na.rm=T),
                       0)
    mismatchAgg <- ifelse(length(which(allReadBases[ReadBaseAvailIdx]!=allRefBases[ReadBaseAvailIdx])) > 0,
                          sum(1-umiReadQuality[which(allReadBases[ReadBaseAvailIdx]!=allRefBases[ReadBaseAvailIdx])],na.rm=T),
                          0)
    
    
    #alignedSegmentAccuracyToReference <- sum(allReadBases[ReadBaseAvailIdx]==allRefBases[ReadBaseAvailIdx]) / length(ReadBaseAvailIdx)
    alignedSegmentAccuracyToReference <- (matchAgg - mismatchAgg) / length(ReadBaseAvailIdx)
    
    alignedSegmentAccuracyToReference <- ifelse(alignedSegmentAccuracyToReference < 0,0,alignedSegmentAccuracyToReference)
    
  }else{
    alignedSegmentAccuracyToReference = 0
  }
  
  alignedSegmentAccuracyToReference
  
}


# Allele likelihood. p(umi|Allele) = sum over n positions (1/n (product of p(base|allele) at each position)); 

probOfBaseGivenAlleleAtPosition <- function(readBases,readQualities,refBase){
  
  selectedIdx <- which(!is.na(readQualities))
  readQualities <- readQualities[selectedIdx]
  readBases <- readBases[selectedIdx]
  
  
  # remove unaligned reads
  #selIdx <- which(readBases != "-")
  #readQualities <- readQualities[selIdx]
  #readBases <- readBases[selIdx]
  
  if(length(readBases) == 0 || length(readBases)==sum(readBases == "-")){
    #no covering base at the position, the key position is not used in the likelihood calculation, no read covers it.
    pUMIbasesGivenAlleleBase = NA
  }else{
    # probability of the match being right is 1-e; while a mismatching base has e/3 prob 
    #probBasesOfEachRead <- sapply(1:length(readBases), function(x) ifelse(readBases[x]==refBase,(1-readQualities[x])* posWt,(readQualities[x]/3) * posWt )) 
    #pUMIbasesGivenAlleleBase <- prod(probBasesOfEachRead,na.rm=T) 
    
    # .... Using the most frequent base frequency at position. 
    # ..... We select the most frequent base among the reads aligning to that position, and use its probability 
    baseFreq <- table(readBases)/length(readBases)
    baseFreq <- baseFreq[names(baseFreq) != "-"]/sum(baseFreq[names(baseFreq) != "-"])
    mostFreqBase <- baseFreq[which.max(baseFreq)]
    
    # we select the most frequenty probability for the most frequenty base
    selectedBase <- names(mostFreqBase)
    selectedBaseProb <- 1- readQualities[which(readBases==selectedBase)]
    selectedBaseProb <- as.numeric(names(sort(table(selectedBaseProb),decreasing=T)[1]))
    
    pUMIbasesGivenAlleleBase <- ifelse(selectedBase==refBase,selectedBaseProb,(1-selectedBaseProb)/3)
    
    
    #.... Or we can calculate a product of the probabilities at the position
    
    # Here we instead of the most frequenty base aligning at the position, a product of the quality of all bases can be taken.
    # and a probability for all bases based on all read bases at the position can be calculated.
    
    refMatchingBasesProb <- function(rb){
      prod(unlist(sapply(1:length(readBases),function(i) ifelse(readBases[i] == rb,1-readQualities[i],readQualities[i]/3))))
    }
    
    
    # unmatchingBasesProb <- function(rb){
    #   prod(unlist(sapply(1:length(readBases),function(i) ifelse(readBases[i] != rb,1-readQualities[i],readQualities[i]/3))))
    # }
    # 
    
    probOfAllBases <- c(refMatchingBasesProb("A"),refMatchingBasesProb("C"),refMatchingBasesProb("G"),refMatchingBasesProb("T"),refMatchingBasesProb("-"))
    names(probOfAllBases) <- c("A","C","G","T","-")
    
    # for now we just use the previous pUMIbasesGivenAlleleBase calculated from most frequent base, but this second approach can be tested.
    #pUMIbasesGivenAlleleBase <- probOfAllBases[which(names(probOfAllBases)==refBase)]
    
  }
  
  pUMIbasesGivenAlleleBase
}


probOfBaseGivenAlleleAtPositionOld <- function(readBases,readQualities,refBase){
  
  # where a read does not cover a key position in the pileup at the position, it is given the least possible base probability
  
  #selectedIdx <- which(!is.na(readQualities))
  #readQualities <- readQualities[selectedIdx]
  #readBases <- readBases[selectedIdx]
  
  
  # remove unaligned reads
  #selIdx <- which(readBases != "-")
  #readQualities <- readQualities[selIdx]
  #readBases <- readBases[selIdx]
  
  
  # probability of the match being right is 1-e; while a mismatching base has e/3 prob 
  #probBasesOfEachRead <- sapply(1:length(readBases), function(x) ifelse(readBases[x]==refBase,(1-readQualities[x])* posWt,(readQualities[x]/3) * posWt )) 
  #pUMIbasesGivenAlleleBase <- prod(probBasesOfEachRead,na.rm=T) 
  
  # for cases where the reads don't cover the site  (no quality record and base at the position), we give all bases a probability of 
  # 7.943282e-05/4 ~ 1.98582e-05, that is the probability of other bases in best called bases with phred-score of 41 ~ 10^(-41/10)
  
  # base frequency at position
  baseFreq <- table(readBases)/length(readBases)
  baseFreq <- baseFreq[names(baseFreq) != "-"]/sum(baseFreq[names(baseFreq) != "-"])
  mostFreqBase <- baseFreq[which.max(baseFreq)]
  
  if(length(mostFreqBase) > 0){
    selectedBase <- names(mostFreqBase)
    selectedBaseProb <- 1- readQualities[which(readBases==selectedBase)]
    selectedBaseProb <- as.numeric(names(sort(table(selectedBaseProb),decreasing=T)[1]))
    
    #readbaseProb <- as.numeric(mostFreqBase)
    pUMIbasesGivenAlleleBase <- ifelse(selectedBase==refBase,selectedBaseProb,(1-selectedBaseProb)/3)
  }else{
    #No aligning base found at the position 
    pUMIbasesGivenAlleleBase <- 1.98582e-05
    
  }
  
  
  pUMIbasesGivenAlleleBase
}


readToAlleleLikelihood <- function(UMIReadsAtPositions,UMIReadQualitiesAtPositions,refAlleleBases){
  
  #readBases = UMIReadsAtPositions[,x]
  #readQualities = UMIReadQualitiesAtPositions[,x]
  #refBase = refAlleleBases[x]
  #UMIReadsAtPositions=allUMIReadsAlmnt.atPositions
  #UMIReadQualitiesAtPositions=allReadsAlmnt.atPositions_Quality
  #refAlleleBases=ReferenceBasesAllele
  
  pUMIbasesGivenAlleleBases <- sapply(1:ncol(UMIReadsAtPositions),function(x) probOfBaseGivenAlleleAtPosition(UMIReadsAtPositions[,x],UMIReadQualitiesAtPositions[,x],refAlleleBases[x]))
  
  if(sum(is.na(pUMIbasesGivenAlleleBases)) == length(pUMIbasesGivenAlleleBases)){
    # Allele likelihood could not be computed because no reads are covering the key sites in the allele.
    pUMIbasesGivenAllele <- 0 
  }else{
    pUMIbasesGivenAllele <- sum(pUMIbasesGivenAlleleBases,na.rm=T)/length(pUMIbasesGivenAlleleBases)
    
    # p(umi|AlleleKeyBases) <- Joint probability of p of base at each key site; product of them and then in phred scale; smallest is allele with the highest likelihood
    #print(-10 * log10(prod(pUMIbasesGivenAlleleBases,na.rm=T)))
    
    #pUMIbasesGivenAllele <- prod(pUMIbasesGivenAlleleBases,na.rm=T) # we don't have same number of key sites between genes. To make the results comparible we use the summation
  }
  
  
  pUMIbasesGivenAllele 
}

getAlleleLikelihoodsForGene <- function(alleleReferencesGene=NULL,tempHLAreferencesAfterAlignment=NULL,polymorphicPositions=NULL,combinedR1R2AlignedReadsForUMI=NULL,combinedR1R2AlignedReads_QualForUMI=NULL){
  
  if(is.null(alleleReferencesGene)){
    stop("Allele References for the gene not provided.")
  }
  
  AlleleLikGene <- c()
  
  for(a in 1:length(alleleReferencesGene)){
    
    #print(a)
    refseq <- as.character(tempHLAreferencesAfterAlignment)[a]
    
    ReferenceBasesAllele <- strsplit(as.character(tempHLAreferencesAfterAlignment)[a],"")[[1]][polymorphicPositions]
    
    allUMIReadsAlmnt <- pairwiseAlignment(pattern = as.character(combinedR1R2AlignedReadsForUMI), subject = as.character(tempHLAreferencesAfterAlignment)[a],type = "overlap")
    #allUMIReadsAlmnt <- pairwiseAlignment(pattern = as.character(contig), subject = as.character(tempHLAreferencesAfterAlignment)[a],type = "local")
    
    # alignment match accuracy for contig
    #checkidx <- which(as.matrix(allUMIReadsAlmnt)[1,] != "-")
    #querypos <- as.matrix(allUMIReadsAlmnt)[1,][checkidx]
    #refpos <- unlist(strsplit(as.character(tempHLAreferencesAfterAlignment)[a],""))[checkidx]
    #sum(querypos == refpos)/length(querypos)
    
    allUMIReadsAlmnt.atPositions <- as.matrix(allUMIReadsAlmnt)[,polymorphicPositions,drop=FALSE]
    
    allUMIReadsAlmnt <- as.matrix(allUMIReadsAlmnt)
    rownames(allUMIReadsAlmnt) <- names(combinedR1R2AlignedReadsForUMI)
    rownames(allUMIReadsAlmnt.atPositions) <- names(combinedR1R2AlignedReadsForUMI) 
    
    readQualityEncoding <- encoding(combinedR1R2AlignedReads_QualForUMI)
    
    
    
    allReadsAlmnt_Quality <- getReadQualityFromALignment(allUMIReadsAlmnt,combinedR1R2AlignedReads_QualForUMI,readQualityEncoding=readQualityEncoding)
    allReadsAlmnt.atPositions_Quality <- getReadQualityFromALignment(allUMIReadsAlmnt,combinedR1R2AlignedReads_QualForUMI,readQualityEncoding=readQualityEncoding,keyPositions=polymorphicPositions)
    
    
    # Accuracy for whole aligning segment needs to be improved, if the size of the aligned portion is too small, even just one correct base, it give high accuracy. 
    UMIreadsAccuracyToRef <- sapply(1:nrow(allUMIReadsAlmnt),function(x) getUMIreadsAccuracy(allUMIReadsAlmnt[x,],tempHLAreferencesAfterAlignment[a]))
    
    UMIreadsPolySitesAccuracyToRef <- sapply(1:nrow(allUMIReadsAlmnt),function(x) getUMIreadsKeySitesAccuracy(allUMIReadsAlmnt[x,],allReadsAlmnt.atPositions_Quality[x,],tempHLAreferencesAfterAlignment[a],polymorphicPositions))
    
    #print(gene)
    #print(dim(allUMIReadsAlmnt.atPositions))
    
    #print(dim(allReadsAlmnt.atPositions_Quality))
    
    #print(UMIreadsAccuracyToRef) # not a very good measure because UMI reads may align nicely to both alleles at different positions.
    #print(UMIreadsPolySitesAccuracyToRef) # a better measure
    
    #print(mean(UMIreadsAccuracyToRef)) # not a very good measure because UMI reads may align nicely to both alleles at different positions.
    #print(mean(UMIreadsPolySitesAccuracyToRef)) # a better measure
    
    AlleleLikelihood <- readToAlleleLikelihood(allUMIReadsAlmnt.atPositions,allReadsAlmnt.atPositions_Quality,ReferenceBasesAllele)
    
    
    # key sites from referene and aligned reads together
    # refSeqAtSites <- paste(ReferenceBasesAllele,collapse="")
    # alignedSeqs <- as.character(apply(allUMIReadsAlmnt.atPositions,1,paste,collapse=""))
    # DPA1_refReadsKeyPos <- c(refSeqAtSites,alignedSeqs)
    # 
    # refSeq <- paste(ReferenceBasesAllele,collapse="")
    # alignedSeqs <- as.character(apply(allUMIReadsAlmnt.atPositions,1,paste,collapse=""))
    # B1_refReadsKeyPos <- c(refSeq,alignedSeqs)
    # 
    # to check alignment in detail
    # refAndReadsAtKeySites <- rbind(ReferenceBasesAllele,allUMIReadsAlmnt.atPositions)
    # write.table(refAndReadsAtKeySites,file="DPA1_refAndReadsAtKeySites.txt")
    # 
    
    
    
    
    # Fix it starting from here:
    # For some reason the two alleles for a gene are getting the same likelihood, at least for this UMI, figure out what is going on and fix it, try another UMI also.
    
    #UMIReadsAtPositions <- allUMIReadsAlmnt.atPositions
    #UMIReadQualitiesAtPositions <- allReadsAlmnt.atPositions_Quality
    # refAlleleBases <- ReferenceBasesAllele
    
    #probOfBaseGivenAlleleAtPosition(UMIReadsAtPositions[,x],UMIReadQualitiesAtPositions[[x]],ReferenceBasesAllele[x])
    #a2 <- sapply(1:ncol(UMIReadsAtPositions),function(x) probOfBaseGivenAlleleAtPosition(UMIReadsAtPositions[,x],UMIReadQualitiesAtPositions[,x],ReferenceBasesAllele[x]))
    #a1 <- sapply(1:ncol(UMIReadsAtPositions),function(x) probOfBaseGivenAlleleAtPosition(UMIReadsAtPositions[,x],UMIReadQualitiesAtPositions[,x],ReferenceBasesAllele[x]))
    
    
    
    
    #AlleleLikelihood_AllSites <- readToAlleleLikelihood(aligningReadSegments,allReadsAlmnt_Quality,refCoveredRegions,keySitesOnly=F)
    
    
    AlleleLikGene <- c(AlleleLikGene,AlleleLikelihood)
    
    
  }
  
  AlleleLikGene
}


calculateUMIuniquenessPerAllele <- function(UMIsList=NULL,FromDir=NULL,fileKeyWord=NULL){
  
  getUniqueRates <- function(UMIsList){
    # get unique UMIs per Allele
    UniqueUMIrates <- c()
    
    for(fn in 1:length(UMIsList)){
      alleleUmis <- as.character(unlist(UMIsList[fn]))
      otherAlleleUmis <- as.character(unlist(UMIsList[-fn]))
      
      AlleleUniqueFraction <- sum(!(alleleUmis %in% otherAlleleUmis))/length(alleleUmis)
      UniqueUMIrates <- c(UniqueUMIrates,AlleleUniqueFraction)
    }
    
    names(UniqueUMIrates) <- names(UMIsList)
    
    UniqueUMIrates
    
  }
  
  if(!is.null(FromDir)){
    umiListfiles <- list.files(path = FromDir, pattern=fileKeyWord,ignore.case=TRUE,full.names=T,recursive=T)
    if(length(umiListfiles) == 0){
      stop("Files could not be found")
    }else{
      UMIsList <- list()
      for(fn in 1:length(umiListfiles)){
        fnName <- unlist(strsplit(umiListfiles[fn],"_"))[2]
        UMIsList[[fnName]] <- read.table(umiListfiles[fn],header=T,stringsAsFactors=F)[,1]
      }
      
      getUniqueRates(UMIsList)
    }
    
  }else{
    
    if(is.null(UMIsList)){
      stop("No list of UMIs given")
    }
    # get uniqueness from UMIsList directly
    getUniqueRates(UMIsList)
  }
  
}

#umiUniquenessInAlleles <- calculateUMIuniquenessPerAllele(FromDir=FromDir,fileKeyWord=fileKeyWord)

#write.table(umiUniquenessInAlleles,file="UMI_uniqueness_in_Alleles_withBeforeReviewScript.txt",sep="\t")


extractUMIsFromReadName <- function(readnames){
  umis <- sapply(strsplit(readnames,"_"),function(x) return(x[2]))
  umis
}

extractUMIsFromReadNameONT <- function(AlignedReadnames,allRawReads,umiLen=10){
  
  # Primers on 5'
  universalPrimerAndTso <- "TTTCTGTTGGTGCTGATATTGCCAGTGGTATCAACGCAGAGT"
  universalPrimerAndTso_RevComp <- "ACTCTGCGTTGATACCACTGGCAATATCAGCACCAACAGAAA"
  
  
  
  
  allReadNames = as.character(id(allRawReads))
  allReadNames = sapply(allReadNames,function(x) as.character(strsplit(x," ")[[1]][1]))
  
  rawReadsForAligningReads = allRawReads[match(unique(AlignedReadnames),allReadNames)]
  
  alleleReadSet = sread(rawReadsForAligningReads)
  
  selectedAligningReadnames <- as.character(id(rawReadsForAligningReads)) # these are aligned read names for which the UMI was detected in the reads.
  selectedAligningReadnames = as.character(sapply(selectedAligningReadnames,function(x) as.character(strsplit(x," ")[[1]][1])))
  
  
  # for now we do this allowing 4 mismatches and no indels (4 mismatch from 10% error rate of the 42bp long primers)
  # first we look for antisense reads and reverse complement them. Then look for all sense reads (which will include reverse complemented antisense reads)
  
  antiSensehits_idx=vcountPattern(universalPrimerAndTso_RevComp,alleleReadSet,max.mismatch=4,with.indels=F)
  alleleReadSet[which(antiSensehits_idx > 0)] <- reverseComplement(alleleReadSet[which(antiSensehits_idx > 0)])
  
  # Now we find all sense reads (which includes the antisense as well)
  sensehits_idx=vcountPattern(universalPrimerAndTso,alleleReadSet,max.mismatch=4,with.indels=F)
  
  primerReads = alleleReadSet[sensehits_idx == 1]
  primerReadsNames = selectedAligningReadnames[sensehits_idx == 1]
  
  noPrimerReads = alleleReadSet[sensehits_idx == 0]
  #noPrimerReadsNames = selectedAligningReadnames[sensehits_idx == 0]
  
  twoPrimerReads = alleleReadSet[sensehits_idx > 1]
  #twoPrimerReadsNames = selectedAligningReadnames[sensehits_idx > 1]
  
  # get the end Of TSO positions in the reads in which the TSO was found
  sensehit_startEnd = vmatchPattern(universalPrimerAndTso,primerReads,max.mismatch=4,with.indels=F)
  
  #endOfTSO = unlist(sapply(sensehit_startEnd,function(x) if(length(x) > 1)end(x[1,]) else end(x)))  This works fine when we want to include reads with more than one TSO content
  
  sensehit_startEnd = unlist(sensehit_startEnd)
  endOfTSO = end(sensehit_startEnd)
  
  # Do UMI extraction based on reads that have the primerTSO sequences. The pattern search above can be improved.
  # select those with reasonable end positions of TSO (they have to be below position 150); This can be left out or changed.
  
  #selectedReads_idx = which(sensehits_idx == 1)
  #senseReadsForAlignment = rawReadsForAligningReads[selectedReads_idx]
  
  selectedEndTSOs = endOfTSO[endOfTSO < 150]
  #senseReadsSelected = senseReadsForAlignment[endOfTSO < 150]
  #selectedReadSet = sread(senseReadsSelected)
  
  senseSelectedReadSet = primerReads[endOfTSO < 150]
  primerReadsNamesSelected <- primerReadsNames[endOfTSO < 150]
  
  selectedReadSet_df = as.data.frame(senseSelectedReadSet)
  
  #collec all UMIs
  
  TotalUMIs = c()
  
  for(j in 1:nrow(selectedReadSet_df)){
    
    # collect umis
    umiStart = selectedEndTSOs[j]+1
    umiWidth = umiLen-1 # passed parameter - 1
    umi=substr(selectedReadSet_df[j,1],umiStart,umiStart+umiWidth)
    TotalUMIs = c(TotalUMIs,umi) 
  }
  
  names(TotalUMIs) <- primerReadsNamesSelected 
  
  TotalUMIs
  
}

collapseUMIs <- function(allUMIs,distCutOff=1){
  # returns a frequency vector of collapsed UMIs. The names of the vector holds the umis. 
  
  #allUMIs <- unique(allAlignedUMIs)
  
  tempumis <- allUMIs
  #tempumiOrigFreq <-sort(table(tempumis),decreasing=T)
  
  umisDedup <- c()
  umiDedupFreq <- c()
  
  tempumiFreq <-sort(table(tempumis),decreasing=T)
  tempumis <- names(tempumiFreq)
  
  while(length(tempumis) > 0){
    umisDedup <- c(umisDedup,tempumis[1])
    umid <- as.vector(adist(tempumis[1],tempumis))
    
    frequencyAfterCollapsing <- sum(tempumiFreq[names(tempumiFreq) %in% tempumis[which(umid <= distCutOff)]])
    umiDedupFreq <- c(umiDedupFreq,frequencyAfterCollapsing)
    
    tempumiFreq <- tempumiFreq[umid > distCutOff]
    tempumis <- names(tempumiFreq)
    
  }
  
  names(umiDedupFreq) <- umisDedup
  
  umiDedupFreq
}


#________________________________________ Run Allele-specific expression Analysis function ____________________________________________ ####

#' Allele-specific HLA gene expression estimation function
#' 
#' This is the main function that performs allele-specific expression estimating for HLA gene. UMIs in the reads are utilized to perform the estimation. 
#' 
#' @param sampleHlaData a space/tab delimited file that contains HLA type information for each every sample. It should contain the sample name on first column, and HLA types 
#' of the sample for the HLA genes for which allele-specific expression estimation is required. 
#' @param sampleId sample name of a single sample, use only when evaluation is needed only for a single sample. If NULL, function is run for all samples in sampleHlaData.  
#' @param IlluminaDir directory where illumina fastq reads data resides. These files should start with the sample name in sampleHlaData followed by underscore. If your data is nanopore reads, this should be NULL (the default)
#' @param nanDir directory where nanopore fastq reads data resides (for each sample, it must be either three files, each for 2D, template,and complement data or one single file combining all three). 
#' These files should start with the sample name in sampleHlaData followed by underscore. If your data is illumina reads, this should be NULL (the default)
#' @param umiLength the length of the UMI sequence (default 10). For illumina reads it should be present at the end of all read names in the fastq files following an underscore. For nanopore ONT reads,the extraction of 
#' the UMI is done in this script assuming the reads contain a universal primer and template switch oligo (TSO) of "TTTCTGTTGGTGCTGATATTGCCAGTGGTATCAACGCAGAGT" at the beginning of the reads (as was in our libraries) and the UMI is found following the TSO sequence.
#' This extraction is done in function extractUMIsFromReadNameONT. Modifications for the TSO and other ways of UMI extraction can be done in that function when needed.  
#' @param genes a vector of HLA genes for which alleles-specific expression estimation is to be done. By default, the following genes are evaluated: c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"). 
#' @param hladbFile Address of the HLA reference database file after indexing with last aligner (e.g hladb/totalHLADb). All other files produced when indexing with last aligner must be present in the same directory as this file. If this address is not provided, the latest version of the HLA reference database from IMGT/HLA(ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta)
#' will be downloaded to the current directory, indexed with last, and used.  
#' @param writeTo directory where results should be written. If not given, current working directory is used.


#' @return Intermediate and final Results are written in a directory inside writeTo for each sample. Function does not return other values. The main result, allele-specific expression counts, is in file named sampleName_UMI_countingReport.txt for each sample (on column named nAlleleUMIs_keySitesFromAllGenes). Detailed information about each UMI evaluated is provided in sampleName_UMIdetails.txt file.
#'
#' 

getAlleleFreqCounts <- function(sampleHlaData=NULL,sampleId=NULL,IlluminaDir=NULL,nanDir = NULL,umiLength=10,genes=NULL,hladbFile=NULL,writeTo=NULL){
  
  if(is.null(sampleHlaData)){
    stop("\t\t No samples data provided. Please provide the address of the file containing sample information to sampleHlaData.")
  }else{
    sampleHlaData = read.table(file=sampleHlaData,header=T,stringsAsFactors=F,sep="\t",na.strings=c(""," ","NA"),check.names=FALSE)
  }
  
  if(is.null(sampleId)){
    sids <- sampleHlaData[,1] # first column in sampleHlaData should contain the sample names. 
  }else{
    sids <- sampleId
  }
  
  if(is.null(genes)){
    # HLA genes of interest for which the anlaysis is going to be done, this could by provided by user
    genes = c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")
  }else if(sum(grepl("^HLA",genes)==F) > 0){
    genes[grepl("^HLA",genes)==F] <- paste("HLA",genes[grepl("^HLA",genes)==F],sep="-")
  }
  
  
  if(is.null(IlluminaDir)){
    if(is.null(nanDir)){
      stop("\t\t No address to read directories provided. Please provide either the nanopore or Illumina reads directory!!")
    }else{
      readType <- "nanopore"
    }
  }else{
    readType <- "illumina"
  }
  
  
  if(is.null(hladbFile)){
    cat("HLA database file not given. HLA database file will downloaded from IMGT/hla, and indexed with lastdb...\n")
    hladbFile <- downloadHLAdb() 
  }
  
  # read in total HLA db
  totalHLAreference <- Biostrings::readDNAStringSet(hladbFile) 
  
  
  
  if(is.null(writeTo)){
    writeTo = getwd()
  }else{
    writeTo = writeTo
  }
  
  
  #register more cores for parallel run
  registerDoParallel(cores=detectCores())
  
  
  # Start anaysis for sample
  
  sampleASE = foreach(sid=sids,.combine=c) %dopar% {
  
  print(paste("Processing sample",sid,sep=" "))
  #sid=sids
  sampleUMICountReport = NULL
  
  sampleName <- toString(trimws(sid))
  
  ### create intermediate and results files directory for sample
  sampleDir <- toString(sampleName)
  dir.create(sampleDir)
  
  writeToForSample <- paste(writeTo,"/",sampleDir,"/",sep="")
  
  
  ### Prepare HLA alleles reference file for sample
  
  AlleleIndicesInRef <- c()
  AllelesForSample <- c()
  
  for(gene in genes){
    
    geneAlleleInTable <- grep(gene,colnames(sampleHlaData))
    
    if(length(geneAlleleInTable) == 0){
      next
    }
    
    
    if(length(geneAlleleInTable) == 2){
      
      
      allele1 <- sampleHlaData[sampleHlaData[,1]==sampleName,geneAlleleInTable[1]]
      allele2 <- sampleHlaData[sampleHlaData[,1]==sampleName,geneAlleleInTable[2]]
      
      if(is.na(allele1) || is.na(allele2)){
        next 
      }
      
    }else if(length(geneAlleleInTable) == 1){
      cat("Less than two alleles found in the Input file for gene ",gene,".This gene is skipped from analysis.\n")
      next
    }
    
    
    
    # Change alleles to 4d
    allele1 <- getAllele4d(allele1)
    allele2 <- getAllele4d(allele2)
    
    if(allele1 == allele2){
      sameAlleles4DigitLevel = T
    }else{
      sameAlleles4DigitLevel = F
    }
    
    
    #find indices of allele sequences in the total HLA reference data
    
    allele1_idx=grep(allele1,names(totalHLAreference),fixed=T)
    allele2_idx=grep(allele2,names(totalHLAreference),fixed=T)
    
    
    if(length(allele1_idx) == 0 | length(allele2_idx) == 0){
      
      cat("\t\t Warning: Allele of gene",gene,"does not exist in the HLA reference database.","\n")
      next 
      
    }
    
    # record allele indices in total reference for all HLA genes
    # we just record the first reference hit for each allele; which is the standard HLA reference sequence at the first hit of the allele at 4d
    
    if(sameAlleles4DigitLevel){
      AlleleIndicesInRef <- c(AlleleIndicesInRef,allele1_idx[1])
      AllelesForSample <- c(AllelesForSample,allele1)
    }else{
      AlleleIndicesInRef <- c(AlleleIndicesInRef,allele1_idx[1],allele2_idx[1])
      AllelesForSample <- c(AllelesForSample,allele1,allele2)
      
    }
    
  }
  
  
  tempHLAReference = totalHLAreference[AlleleIndicesInRef]
  
  tempReferenceFile <- paste(writeToForSample,"tempHLADb",sep="")
  
  writeXStringSet(tempHLAReference, file=tempReferenceFile,format="fasta")
  
  lastDbIndexed <- prepareLastIndex(tempReferenceFile)
  
  
  ### Prepare and align sample reads to the reference Allele sequences; 
  if(readType=="illumina"){
    
    
    ReadAdds <- findIlluminaReadAddress(sampleName, IlluminaDir)
    
    R1_reads <- ReadAdds[1]
    R2_reads <- ReadAdds[2]

    # Uncomment if read file are needed to also be present in the result directory in fastq format.
    #R1con <- file(R1_reads)
    #R2con <- file(R2_reads)
    
    #if(summary(R1con)$class == "gzfile" & summary(R2con)$class== "gzfile"){
    
    #  system(paste("zcat",R1_reads,">",paste(writeToForSample,"R1.fastq",sep=""),sep=" "))
    #  system(paste("zcat",R2_reads,">",paste(writeToForSample,"R2.fastq",sep=""),sep=" "))
    
    #  R1_reads <- paste(writeToForSample,"R1.fastq",sep="")
    #  R2_reads <- paste(writeToForSample,"R2.fastq",sep="")
    #}
    
    #close.connection(R1con)
    #close.connection(R2con)
    
    
    alignmentToTotalHLA <- performAlignmentToTotalHLAData(R1_reads,R2_reads,nR=NULL,tempReferenceFile)
    allRawR1Reads=readFastq(R1_reads)
    allRawR2Reads=readFastq(R2_reads)
    
    
  }else{
    
    nanRFile <- prepareNanoporeData(sampleName,nanDir,writeToForSample)
    
    alignmentToTotalHLA <- performAlignmentToTotalHLAData(R1=NULL,R2File,nanRFile,tempReferenceFile)
    
    #nan_reads <- paste(writeToForSample,nanRFile,sep="")
    nan_reads <- nanRFile
    allRawR1Reads=readFastq(nan_reads)
    
  }
  
  
  ### extract aligning read segments and their quality
  
  if(readType=="illumina"){
    allAligned_R1 <- getAlignedReadNamesAndReads(paste(writeToForSample,"tempAlmt_R1",sep=""))
    allAligned_R2 <- getAlignedReadNamesAndReads(paste(writeToForSample,"tempAlmt_R2",sep=""))
    allAlignedPE <- getAlignedReadNamesAndReads(paste(writeToForSample,"tempAlmt",sep=""))
    #print(is.null(allAlignedPE))
    if(!is.null(allAlignedPE))
      names(allAlignedPE[[1]]) <- sapply(strsplit(names(allAlignedPE[[1]]),"/"),function(x) return(x[1]))
    
    # Collect all reads and their qualities that aligned to the gene
    combinedR1R2AlignedReads <- DNAStringSet()
      if(!is.null(allAligned_R1)) 
        combinedR1R2AlignedReads = c(combinedR1R2AlignedReads,allAligned_R1[[1]])
      if(!is.null(allAligned_R2)) 
        combinedR1R2AlignedReads = c(combinedR1R2AlignedReads,allAligned_R2[[1]])
      if(!is.null(allAlignedPE)) 
        combinedR1R2AlignedReads = c(combinedR1R2AlignedReads,allAlignedPE[[1]])
    
                            
    combinedR1R2AlignedReads_Qual <- c() 
      if(!is.null(allAligned_R1)) 
        combinedR1R2AlignedReads_Qual = c(combinedR1R2AlignedReads_Qual,allAligned_R1[[2]])
      if(!is.null(allAligned_R2)) 
        combinedR1R2AlignedReads_Qual = c(combinedR1R2AlignedReads_Qual,allAligned_R2[[2]])
      if(!is.null(allAlignedPE)) 
        combinedR1R2AlignedReads_Qual = c(combinedR1R2AlignedReads_Qual,allAlignedPE[[2]])
      
    combinedR1R2AlignedReads_Qual <- FastqQuality(BStringSet(combinedR1R2AlignedReads_Qual))
    
    
    combinedR1R2AlignedRead_RefAlleles <- c()
      if(!is.null(allAligned_R1)) 
        combinedR1R2AlignedRead_RefAlleles = c(combinedR1R2AlignedRead_RefAlleles,allAligned_R1[[3]])
      if(!is.null(allAligned_R2)) 
        combinedR1R2AlignedRead_RefAlleles = c(combinedR1R2AlignedRead_RefAlleles,allAligned_R2[[3]])
      if(!is.null(allAlignedPE)) 
        combinedR1R2AlignedRead_RefAlleles = c(combinedR1R2AlignedRead_RefAlleles,allAlignedPE[[3]])
      
  
    combinedR1R2AlignedRead_alRefStart <- c()
      if(!is.null(allAligned_R1)) 
        combinedR1R2AlignedRead_alRefStart = c(combinedR1R2AlignedRead_alRefStart,allAligned_R1[[4]])
      if(!is.null(allAligned_R2)) 
        combinedR1R2AlignedRead_alRefStart = c(combinedR1R2AlignedRead_alRefStart,allAligned_R2[[4]])
      if(!is.null(allAlignedPE)) 
        combinedR1R2AlignedRead_alRefStart = c(combinedR1R2AlignedRead_alRefStart,allAlignedPE[[4]])
      
    
    combinedR1R2AlignedRead_alRefEnd <- c()
      if(!is.null(allAligned_R1)) 
        combinedR1R2AlignedRead_alRefEnd = c(combinedR1R2AlignedRead_alRefEnd,allAligned_R1[[5]])
      if(!is.null(allAligned_R2)) 
        combinedR1R2AlignedRead_alRefEnd = c(combinedR1R2AlignedRead_alRefEnd,allAligned_R2[[5]])
      if(!is.null(allAlignedPE)) 
        combinedR1R2AlignedRead_alRefEnd = c(combinedR1R2AlignedRead_alRefEnd,allAlignedPE[[5]])
     

  }else{
    
    allAlignednR <- getAlignedReadNamesAndReads(paste(writeToForSample,"tempAlmt",sep=""))
    
    combinedR1R2AlignedReads <- c(allAlignednR[[1]])
    combinedR1R2AlignedReads_Qual = c(allAlignednR[[2]])
    combinedR1R2AlignedReads_Qual <- FastqQuality(BStringSet(combinedR1R2AlignedReads_Qual))
    
    combinedR1R2AlignedRead_RefAlleles <- c(allAlignednR[[3]])
    combinedR1R2AlignedRead_alRefStart <- c(allAlignednR[[4]])
    combinedR1R2AlignedRead_alRefEnd <- c(allAlignednR[[5]])
    
    
  }
  
  
  # The above can be done with the following three lines as well
  
  #allAlignedReadsCollected <- collectAllAlignedReads(readType,writeToForSample)
  #combinedR1R2AlignedReads <- allAlignedReadsCollected[[1]]
  #combinedR1R2AlignedReads_Qual = allAlignedReadsCollected[[2]]
  
  totalNumberOfAlignedReads <- length(unique(names(combinedR1R2AlignedReads)))
  

  # UMI assignment is done for each UMI
  
  ### First extract UMIs from the aligning reads, collapse UMIs
  
  if(totalNumberOfAlignedReads != 0){
    
    # there are aligned reads to HLA genes. UMI extraction can be done.
    allAlignedReadNames <- names(combinedR1R2AlignedReads)
    
    if(readType=="illumina"){
      allAlignedUMIs <- extractUMIsFromReadName(allAlignedReadNames)
    }else{
      allAlignedUMIs <- extractUMIsFromReadNameONT(AlignedReadnames=allAlignedReadNames,allRawReads=allRawR1Reads,umiLen=umiLength)
      
      aligningReadsWithUMIDetected <- names(allAlignedUMIs)
      toBekeptReadIdx <- allAlignedReadNames %in% aligningReadsWithUMIDetected
      
      # For ONT samples, keep as aligning reads only those for which the UMI was detected in the reads.
      # and update all relevant variables based on this data
      
      combinedR1R2AlignedReads <- combinedR1R2AlignedReads[toBekeptReadIdx]
      combinedR1R2AlignedReads_Qual = combinedR1R2AlignedReads_Qual[toBekeptReadIdx]
      
      combinedR1R2AlignedRead_RefAlleles <- combinedR1R2AlignedRead_RefAlleles[toBekeptReadIdx]
      combinedR1R2AlignedRead_alRefStart <- combinedR1R2AlignedRead_alRefStart[toBekeptReadIdx]
      combinedR1R2AlignedRead_alRefEnd <- combinedR1R2AlignedRead_alRefEnd[toBekeptReadIdx]
      
      # update variables that are used later
      allAlignedReadNames <- names(combinedR1R2AlignedReads)
      allAlignedUMIs1 <- as.character(sapply(allAlignedReadNames,function(x) allAlignedUMIs[match(x,aligningReadsWithUMIDetected)]))
      allAlignedUMIs <- allAlignedUMIs1 
      
    }
    
    
    
    nOfAlleles <- length(AllelesForSample)
    totalNUmberOfReadsInSample <- length(allRawR1Reads)
    totalNumberOfAlignedReads <- length(unique(names(combinedR1R2AlignedReads)))
    
    
    nAllUMIs <- length(unique(allAlignedUMIs))
    
    
    dedupUMIs <- collapseUMIs(allAlignedUMIs,distCutOff=1)
    
    totalNumberOfUniqueUMIs <- length(dedupUMIs)
    
    #umiDuplicationRateAtGeneLevel <- nAllUMIs/totalNumberOfUniqueUMIs    # this is wrong because nAllUMIs counts multiple reads with same umi as one.
    umiDuplicationRateAtGeneLevel <- totalNumberOfAlignedReads/totalNumberOfUniqueUMIs
    
    # Report total reads in sample,  total number of aligned reads to all HLA alleles, total uncollapsed and collapsed unique UMIs from all HLA alleles.
    sampleUMICountReport <- cbind(sampleUMICountReport,AllelesForSample,rep(totalNUmberOfReadsInSample,nOfAlleles),rep(totalNumberOfAlignedReads,nOfAlleles),
                                  rep(nAllUMIs,nOfAlleles),rep(totalNumberOfUniqueUMIs,nOfAlleles),rep(umiDuplicationRateAtGeneLevel,nOfAlleles)) 
    
    
    # get read indices of same/similar UMIs for each deduplicated UMI
    dedupUMIseqs <- names(dedupUMIs)
    

    if(readType=="illumina"){
      allRawR1ReadsNames <- as.character(id(allRawR1Reads))
      allRawR1ReadsNames = sapply(allRawR1ReadsNames,function(x) as.character(strsplit(x," ")[[1]][1]))
      
      alignedRawReadsIdx <- match(unique(allAlignedReadNames),allRawR1ReadsNames)
      
      allRawAlignedR1ReadsTemp = allRawR1Reads[alignedRawReadsIdx]
      allRawAlignedR2ReadsTemp = allRawR2Reads[alignedRawReadsIdx]
      
      allRawAlignedR1Reads <- sread(allRawAlignedR1ReadsTemp)
      names(allRawAlignedR1Reads) <- id(allRawAlignedR1ReadsTemp)
      
      allRawAlignedR2Reads <- sread(allRawAlignedR2ReadsTemp)
      names(allRawAlignedR2Reads) <- id(allRawAlignedR2ReadsTemp)
      
    }
    
    

    ### Get polymorphic sites between alleles of each gene; do this here and keep it in a list. Saves time. 
    
    ReferenceAlignmentsAndPolymorphicSites <- list()
    for(gn in genes){
      
      gn <- ifelse(grepl("-",gn),unlist(strsplit(gn,"-"))[2],gn)
      
      referencesForGene <- tempHLAReference[grep(paste("^",gn,sep=""),names(tempHLAReference),fixed=F)]
      
      #print(gn)
      
      #print(referencesForGene)
      if(length(referencesForGene)==0){
        # no given alleles for the gene
        next
      }else if(length(referencesForGene)==1){
        
        # homozygous
        # We find other alleles for the gene that are different from this sample's allele at 4d. We define the key sites that make this alleles different by multiple alignment.
        
        homozygousAllele <- AllelesForSample[grep(paste("^",gn,sep=""),AllelesForSample)]
        homoAlleleLength <- as.numeric(unlist(strsplit(names(referencesForGene)," "))[2])
        
        geneInRef_idx=grep(paste("^",gn,sep=""),names(totalHLAreference),fixed=F)
        allAllelesForGene <- totalHLAreference[geneInRef_idx]
        differentAllelesForgene <- allAllelesForGene[!grepl(homozygousAllele,names(allAllelesForGene),fixed=T)]
        
        diffAllelesLength <- sapply(1:length(differentAllelesForgene),function(x) as.numeric(unlist(strsplit(names(differentAllelesForgene[x])," "))[2]))
        sameLengthdifferentAllelesForgene <- differentAllelesForgene[which(diffAllelesLength==homoAlleleLength)]
        
        selectedDifferentAlleleForgene <- differentAllelesForgene[1]
        if(length(sameLengthdifferentAllelesForgene) > 0)
          selectedDifferentAlleleForgene <- sameLengthdifferentAllelesForgene[1]
        
        
        referencesForGeneAndOtherAlleles <- c(referencesForGene,selectedDifferentAlleleForgene)
        
        
        
        #print(referencesForGeneAndOtherAlleles)
        
        
        alignmentAndPolySitesForReferences <- getPolymorphicSitesBetweenAlleles(referenceSeqs=referencesForGeneAndOtherAlleles)
        
        ReferenceAlignmentsAndPolymorphicSites[[gn]] <- alignmentAndPolySitesForReferences
        
        # referenceAlignment <- alignmentAndPolySitesForReferences[[1]]
        # differingPositions <- alignmentAndPolySitesForReferences[[2]]
        # 
        # geneAllelesLik <- getAlleleLikelihoodsForGene(alleleReferencesGene=referencesForGene,tempHLAreferencesAfterAlignment=referenceAlignment,
        #                                               polymorphicPositions=differingPositions,
        #                                               combinedR1R2AlignedReadsForUMI=combinedR1R2AlignedReadsForUMI,
        #                                               combinedR1R2AlignedReads_QualForUMI=combinedR1R2AlignedReads_QualForUMI)
        # 
        
      }else{
        
        # heterozygous
        #referencesForGeneAndOtherAlleles <- c(referencesForGene,selectedDifferentAlleleForgene)
        
        alignmentAndPolySitesForReferences <- getPolymorphicSitesBetweenAlleles(referenceSeqs=referencesForGene)
        
        ReferenceAlignmentsAndPolymorphicSites[[gn]] <- alignmentAndPolySitesForReferences
        
        
        # referenceAlignment <- alignmentAndPolySitesForReferences[[1]]
        # differingPositions <- alignmentAndPolySitesForReferences[[2]]
        # 
        # reference bases at the key polymorphic positions between the two reference alleles
        #ReferenceBasesAllele1 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[1],"")[[1]][polymorphicPositions]
        #ReferenceBasesAllele2 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[2],"")[[1]][polymorphicPositions]
        
        
        # geneAllelesLik <- getAlleleLikelihoodsForGene(alleleReferencesGene=referencesForGene,tempHLAreferencesAfterAlignment=referenceAlignment,
        #                                               polymorphicPositions=differingPositions,
        #                                               combinedR1R2AlignedReadsForUMI=combinedR1R2AlignedReadsForUMI,
        #                                               combinedR1R2AlignedReads_QualForUMI=combinedR1R2AlignedReads_QualForUMI)
        # 
        # 
        
        
        
      }
      
      
      
      
    }
    
    
    
    ###
    
    
    
    AssignedUMIsList <- list()
    AssignedUMIsList_usingUniqueAligningUMIs <- list()
    AssignedUMIsList_usingAllReferencekeySites <- list()
    
    UnassignedUMIs <- c() # UMI does not align to key sites in all alleles and gets unassigned.
    UniqueAligningUMIs <- 0
    
    UMIdetails <- NULL
    
    #length(dedupUMIseqs)
    #print("Assigning UMIs to alleles:")
    for(umidx in 1:length(dedupUMIseqs)){
      
      referencesForGene <- NULL
      
      #print(umidx)
      #currentUMI <- "CATGGCGCCC"
      currentUMI <- dedupUMIseqs[umidx]
      #print(currentUMI)
      
      
      
      sameUMI_idx <- which(allAlignedUMIs == currentUMI) 
      allAlignedReadsForUMI <- allAlignedReadNames[sameUMI_idx]
      #similarUMIFreq <- table(extractUMIsFromReadName(allAlignedReadsForUMI))
      
      # aligned read fragments for UMI reads
      combinedR1R2AlignedReadsForUMI <- combinedR1R2AlignedReads[sameUMI_idx]
      combinedR1R2AlignedReads_QualForUMI <- combinedR1R2AlignedReads_Qual[sameUMI_idx]
      
      combinedR1R2AlignedRead_RefAllelesForUMI <- combinedR1R2AlignedRead_RefAlleles[sameUMI_idx]
      combinedR1R2AlignedRead_alRefStartForUMI <- combinedR1R2AlignedRead_alRefStart[sameUMI_idx]
      combinedR1R2AlignedRead_alRefEndForUMI <- combinedR1R2AlignedRead_alRefEnd[sameUMI_idx]
      
      # evaluate UMI. Check if it looks real or not.
      if(readType=="illumina"){
        # raw reads for UMI reads.
        allAlignedRawReadsForUMI_idx <- unique(sapply(1:length(allAlignedReadsForUMI),function(x) grep(allAlignedReadsForUMI[x],names(allRawAlignedR1Reads),fixed=T)))
        
        allAlignedR1RawReadsForUMI <- allRawAlignedR1Reads[allAlignedRawReadsForUMI_idx]
        allAlignedR2RawReadsForUMI <- allRawAlignedR2Reads[allAlignedRawReadsForUMI_idx]
        
        
        
        # is UMI found in the reference sequence? This is highly unlikely, thus the UMI might be wrong or part of the real sequence.
        umiMatchingRefs <- sapply(1:length(tempHLAReference),function(x) grepl(currentUMI,as.character(tempHLAReference[[x]]),fixed=T))
        
        nUMIMatchingRefs <- sum(umiMatchingRefs)
        
        # is UMI found in in the total aligned read sequences? could also be that it is a real UMI but unextracted for some reason
        
        umiMatchingReads <- sapply(1:length(combinedR1R2AlignedReads),function(x) grepl(currentUMI,as.character(combinedR1R2AlignedReads[[x]]),fixed=T))
        
        umiMatchingR1RawReads <- sapply(1:length(allRawAlignedR1Reads),function(x) grepl(currentUMI,as.character(allRawAlignedR1Reads[[x]]),fixed=T))
        umiMatchingR2RawReads <- sapply(1:length(allRawAlignedR2Reads),function(x) grepl(currentUMI,as.character(allRawAlignedR2Reads[[x]]),fixed=T))
        
        # readnames are not unique in the combinedR1R2AlignedReads
        nUMIMatchingReads <- max(length(unique(names(combinedR1R2AlignedReads[umiMatchingReads]))),sum(umiMatchingR1RawReads),sum(umiMatchingR2RawReads))
        
        
        
        # is UMI found in in the aligned reads with the UMI? It doesnt have to be because it is extracted as UMI from these reads.
        umiMatchingReadsForUMI <- sapply(1:length(combinedR1R2AlignedReadsForUMI),function(x) grepl(currentUMI,as.character(combinedR1R2AlignedReadsForUMI[[x]]),fixed=T))
        
        umiMatchingR1RawReadsForUMI <- sapply(1:length(allAlignedR1RawReadsForUMI),function(x) grepl(currentUMI,as.character(allAlignedR1RawReadsForUMI[[x]]),fixed=T))
        umiMatchingR2RawReadsForUMI <- sapply(1:length(allAlignedR2RawReadsForUMI),function(x) grepl(currentUMI,as.character(allAlignedR2RawReadsForUMI[[x]]),fixed=T))
        
        nUMIMatchingReadsForUMI <- max(length(unique(names(combinedR1R2AlignedReadsForUMI[umiMatchingReadsForUMI]))),sum(umiMatchingR1RawReadsForUMI),sum(umiMatchingR2RawReadsForUMI))
        
        # do all reads having the UMI start with the 3 Gs ?
        
        gggStartMatchingReadsForUMI <- sapply(1:length(allAlignedR1RawReadsForUMI),function(x) grepl("^GGG",as.character(allAlignedR1RawReadsForUMI[[x]]),fixed=F))
        
        nNonGGGStartMatchingUMIReads <- sum(gggStartMatchingReadsForUMI == F) # all reads should start with GGG or zero number of non GGG starting reads is expected.
        
        # combine hit counts from all
        UMIAgg <- sum(nUMIMatchingRefs,nUMIMatchingReads,nUMIMatchingReadsForUMI,nNonGGGStartMatchingUMIReads)
        UMIstatus <- ifelse(UMIAgg==0,"Real","Not Real")
        
        toUMIdetails <- c(currentUMI,nUMIMatchingRefs,nUMIMatchingReads,nUMIMatchingReadsForUMI,nNonGGGStartMatchingUMIReads,UMIstatus)
        
      }else{
        
        # for ONT reads extracted UMIs are always considered real
        nUMIMatchingRefs <- 0
        nUMIMatchingReads <- 0
        nUMIMatchingReadsForUMI <- 0
        nNonGGGStartMatchingUMIReads <- 0
        
        UMIAgg <- sum(nUMIMatchingRefs,nUMIMatchingReads,nUMIMatchingReadsForUMI,nNonGGGStartMatchingUMIReads)
        UMIstatus <- ifelse(UMIAgg==0,"Real","Not Real")
        
        toUMIdetails <- c(currentUMI,nUMIMatchingRefs,nUMIMatchingReads,nUMIMatchingReadsForUMI,nNonGGGStartMatchingUMIReads,UMIstatus)
        
      }
      
      
      if(UMIAgg > 0){
        toUMIdetails <- c(toUMIdetails,NA,NA)
        UMIdetails <- rbind(UMIdetails,toUMIdetails)
        next
      }
      
      
      hitRefAllelesForUMI <- unique(combinedR1R2AlignedRead_RefAllelesForUMI)
      
      HLAgenesInReference <- sort(unique(sapply(strsplit(hitRefAllelesForUMI,"\\*"),function(x) x[1])))
      
      selectedAllelesForSample <- AllelesForSample[unlist(lapply(HLAgenesInReference,function(x) grep(paste("^",x,sep=""),AllelesForSample,fixed=F)))]
      
      
      #print("Allele hit by aligner:")
      #print(hitRefAllelesForUMI)
      
      if(length(hitRefAllelesForUMI) == 1){
        
        
        # UMI aligns only to one allele. No need to evaluate other alleles.
        mostProbAleleForUMI <- selectedAllelesForSample[which(sapply(selectedAllelesForSample,function(x) grepl(x,hitRefAllelesForUMI,fixed=T)))]
        mostProbAleleForUMI_usingAllReferencekeySites <- mostProbAleleForUMI
        
        UniqueAligningUMIs <- UniqueAligningUMIs + 1
        
        if(mostProbAleleForUMI %in% names(AssignedUMIsList)){
          AssignedUMIsList[[mostProbAleleForUMI]] <- c(AssignedUMIsList[[mostProbAleleForUMI]],currentUMI)
        }else{
          AssignedUMIsList[[mostProbAleleForUMI]] <- currentUMI
        }
        
        # count using only unique aligning UMIs
        if(mostProbAleleForUMI %in% names(AssignedUMIsList_usingUniqueAligningUMIs)){
          AssignedUMIsList_usingUniqueAligningUMIs[[mostProbAleleForUMI]] <- c(AssignedUMIsList_usingUniqueAligningUMIs[[mostProbAleleForUMI]],currentUMI)
        }else{
          AssignedUMIsList_usingUniqueAligningUMIs[[mostProbAleleForUMI]] <- currentUMI
        }
        
        
        # using key sites from all references
        if(mostProbAleleForUMI_usingAllReferencekeySites %in% names(AssignedUMIsList_usingAllReferencekeySites)){
          AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]] <- c(AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]],currentUMI)
        }else{
          AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]] <- currentUMI
        }
        
        
        
        
        toUMIdetails <- c(toUMIdetails,paste(hitRefAllelesForUMI,sep=";"),mostProbAleleForUMI)
        UMIdetails <- rbind(UMIdetails,toUMIdetails)
        
        next
        
      }
      
      #HLAgenesInReference <- unique(sapply(strsplit(AllelesForSample,"\\*"),function(x) x[1]))
      
      #AllelesForSample <- sort(AllelesForSample)
      
      
      
      # Prepare UMI reads for assembly
      # allRawR1ReadsNames <- as.character(id(allRawR1Reads))
      # allRawR1ReadsNames = sapply(allRawR1ReadsNames,function(x) as.character(strsplit(x," ")[[1]][1]))
      # rawRRawR1ForUMI = allRawR1Reads[match(unique(allAlignedReadsForUMI),allRawR1ReadsNames)]
      # 
      # allRawR2ReadsNames <- as.character(id(allRawR2Reads))
      # allRawR2ReadsNames = sapply(allRawR2ReadsNames,function(x) as.character(strsplit(x," ")[[1]][1]))
      # rawRRawR2ForUMI = allRawR2Reads[match(unique(allAlignedReadsForUMI),allRawR2ReadsNames)]
      # 
      # writeUMIR1 <- paste(writeToForSample,currentUMI,"_R1.fastq.gz",sep="")
      # writeFastq(rawRRawR1ForUMI, writeUMIR1)
      # 
      # writeUMIR2 <- paste(writeToForSample,currentUMI,"_R2.fastq.gz",sep="")
      # writeFastq(rawRRawR2ForUMI, writeUMIR2)
      
      # Velvet contigs
      #contig <- "TGGTATCAACGCAGAGTCATGGCGCCCCGAACCCTCCTCCTGCTGCTCTCGGGAGCCCTGGCCCTGACCGAGACCTGGGCCTGCTCCCACTCCATGAGGTATTTCTACACCGCTGTGTCCCGGCCCAGCCGCGGAGAGCCCCACTTCATCGCAGTGGGCTACGTGGACGACACGC"
      #cont2 <- "TGGTATCAACGCAGAGTCATGGCGCCCCGAACCCTCCTCCTGCTTCTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCC"
      #contig <- "AGTGAACCTGCGGAAACTGCGCGGCTACTACAACCAGAGCGAGGACGGGTCTCACACCCTCCAGTGGATGTTTGGCTGCGACCTGGGGCCGGACGGGCGCCTCCTCCGCGGGTATAACCAGTTCGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGATCTGCGCTCCTGGACCGCCGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGAGGC"
      
      # use key sites gene by gene.
      AlleleLikelihoods_keySites <- c()
      
      #print("Key site likelihood calculation")
      
      for(gene in HLAgenesInReference){
        
        #print(gene)
        
        referencesForGene <- tempHLAReference[grep(paste("^",gene,sep=""),names(tempHLAReference),fixed=F)]
        
        #print("number of references")
        
        #print(length(referencesForGene))
        # 
        # if(length(referencesForGene)==1){
        #   
        # homozygous
        # We find other alleles for the gene that are different from this sample's allele at 4d. We define the key sites that make this alleles different by multiple alignment.
        
        # homozygousAllele <- AllelesForSample[grep(paste("^",gene,sep=""),AllelesForSample)]
        # homoAlleleLength <- as.numeric(unlist(strsplit(names(referencesForGene)," "))[2])
        # 
        # geneInRef_idx=grep(paste("^",gene,sep=""),names(totalHLAreference),fixed=F)
        # allAllelesForGene <- totalHLAreference[geneInRef_idx]
        # differentAllelesForgene <- allAllelesForGene[!grepl(homozygousAllele,names(allAllelesForGene),fixed=T)]
        # 
        # diffAllelesLength <- sapply(1:length(differentAllelesForgene),function(x) as.numeric(unlist(strsplit(names(differentAllelesForgene[x])," "))[2]))
        # sameLengthdifferentAllelesForgene <- differentAllelesForgene[which(diffAllelesLength==homoAlleleLength)]
        # 
        # selectedDifferentAlleleForgene <- differentAllelesForgene[1]
        # if(length(sameLengthdifferentAllelesForgene) > 0)
        #   selectedDifferentAlleleForgene <- sameLengthdifferentAllelesForgene[1]
        # 
        # 
        # referencesForGeneAndOtherAlleles <- c(referencesForGene,selectedDifferentAlleleForgene)
        # 
        # 
        
        #print(referencesForGeneAndOtherAlleles)
        
        
        #alignmentAndPolySitesForReferences <- getPolymorphicSitesBetweenAlleles(referenceSeqs=referencesForGeneAndOtherAlleles)
        
        alignmentAndPolySitesForReferences <- ReferenceAlignmentsAndPolymorphicSites[[gene]]
        referenceAlignment <- alignmentAndPolySitesForReferences[[1]]
        differingPositions <- alignmentAndPolySitesForReferences[[2]]
        
        geneAllelesLik <- getAlleleLikelihoodsForGene(alleleReferencesGene=referencesForGene,tempHLAreferencesAfterAlignment=referenceAlignment,
                                                      polymorphicPositions=differingPositions,
                                                      combinedR1R2AlignedReadsForUMI=combinedR1R2AlignedReadsForUMI,
                                                      combinedR1R2AlignedReads_QualForUMI=combinedR1R2AlignedReads_QualForUMI)
        
        
        # }else{
        
        # heterozygous
        #referencesForGeneAndOtherAlleles <- c(referencesForGene,selectedDifferentAlleleForgene)
        
        #alignmentAndPolySitesForReferences <- getPolymorphicSitesBetweenAlleles(referenceSeqs=referencesForGene)
        
        # alignmentAndPolySitesForReferences <- ReferenceAlignmentsAndPolymorphicSites[[gene]]
        # referenceAlignment <- alignmentAndPolySitesForReferences[[1]]
        # differingPositions <- alignmentAndPolySitesForReferences[[2]]
        # 
        # # reference bases at the key polymorphic positions between the two reference alleles
        # #ReferenceBasesAllele1 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[1],"")[[1]][polymorphicPositions]
        # #ReferenceBasesAllele2 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[2],"")[[1]][polymorphicPositions]
        # 
        # 
        # geneAllelesLik <- getAlleleLikelihoodsForGene(alleleReferencesGene=referencesForGene,tempHLAreferencesAfterAlignment=referenceAlignment,
        #                             polymorphicPositions=differingPositions,
        #                             combinedR1R2AlignedReadsForUMI=combinedR1R2AlignedReadsForUMI,
        #                             combinedR1R2AlignedReads_QualForUMI=combinedR1R2AlignedReads_QualForUMI)
        # 
        # 
        
        
        
        # }
        
        AlleleLikelihoods_keySites <- c(AlleleLikelihoods_keySites,geneAllelesLik)
        
        
        
      }
      
      #print("Now do all key result")
      
      # use key polymorphic sites from all alleles across genes to which the UMI aligned
      getAlleleLikelihoodsUsingAllAlleles <- function(){
        
        allAlleleReferences <- tempHLAReference[unlist(lapply(1:length(HLAgenesInReference),function(x) grep(paste("^",HLAgenesInReference[x],sep=""),names(tempHLAReference),fixed=F)))]
        
        #print(allAlleleReferences)
        alignmentAndPolySitesForReferences <- getPolymorphicSitesBetweenAlleles(referenceSeqs=allAlleleReferences)
        
        referenceAlignment <- alignmentAndPolySitesForReferences[[1]]
        differingPositions <- alignmentAndPolySitesForReferences[[2]]
        
        # reference bases at the key polymorphic positions between the two reference alleles
        #ReferenceBasesAllele1 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[1],"")[[1]][polymorphicPositions]
        #ReferenceBasesAllele2 <- strsplit(as.character(tempHLAreferencesAfterAlignment)[2],"")[[1]][polymorphicPositions]
        
        
        allAllelesLik <- getAlleleLikelihoodsForGene(alleleReferencesGene=allAlleleReferences,tempHLAreferencesAfterAlignment=referenceAlignment,
                                                     polymorphicPositions=differingPositions,
                                                     combinedR1R2AlignedReadsForUMI=combinedR1R2AlignedReadsForUMI,
                                                     combinedR1R2AlignedReads_QualForUMI=combinedR1R2AlignedReads_QualForUMI)
        
        allAllelesLik
        
      }
      
      AlleleLikelihoods_usingAllReferencekeySites <- getAlleleLikelihoodsUsingAllAlleles()
      
      #print("all key result")
      #print(AlleleLikelihoods_usingAllReferencekeySites)
      
      # get most likely alleles
      ProbAllelesGivenData <- AlleleLikelihoods_keySites/sum(AlleleLikelihoods_keySites)
      
      ProbAllelesGivenData_usingAllReferencekeySites <- AlleleLikelihoods_usingAllReferencekeySites/sum(AlleleLikelihoods_usingAllReferencekeySites)
      
      
      # here check if more than one alleles have maximum possible probability and decide what to do, because right now the first one is being picked.
      mostProbAleleForUMI <- selectedAllelesForSample[which.max(ProbAllelesGivenData)]
      
      mostProbAleleForUMI_usingAllReferencekeySites <- selectedAllelesForSample[which.max(ProbAllelesGivenData_usingAllReferencekeySites)]
      
      
      if(length(mostProbAleleForUMI) == 0 && length(mostProbAleleForUMI_usingAllReferencekeySites)==0){
        # UMI is unassigned as its best likeliy allele could not be found. It did not align to key sites well.
        UnassignedUMIs <- currentUMI
        next
      }else{
        if(length(mostProbAleleForUMI) > 0){
          if(mostProbAleleForUMI %in% names(AssignedUMIsList)){
            #print(mostProbAleleForUMI %in% names(AssignedUMIsList))
            AssignedUMIsList[[mostProbAleleForUMI]] <- c(AssignedUMIsList[[mostProbAleleForUMI]],currentUMI)
          }else{
            AssignedUMIsList[[mostProbAleleForUMI]] <- currentUMI
          }
        }
        
        if(length(mostProbAleleForUMI_usingAllReferencekeySites) > 0){
          # using key sites from all references
          if(mostProbAleleForUMI_usingAllReferencekeySites %in% names(AssignedUMIsList_usingAllReferencekeySites)){
            AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]] <- c(AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]],currentUMI)
          }else{
            AssignedUMIsList_usingAllReferencekeySites[[mostProbAleleForUMI_usingAllReferencekeySites]] <- currentUMI
          }
        }
        
      }
      
      
      # allele hits from aligner are separated by ;. And allele hits from likelihood calculation with key sites determined at gene and all genes levels separated by |.
      toUMIdetails <- c(toUMIdetails,paste(hitRefAllelesForUMI,collapse=";"),paste(mostProbAleleForUMI,mostProbAleleForUMI_usingAllReferencekeySites,sep="|"))
      #print(toUMIdetails)
      UMIdetails <- rbind(UMIdetails,toUMIdetails)
      
      #print("************")
      
      
    }
    

    
    # add number of single allele aligning UMIs
    sampleUMICountReport <- cbind(sampleUMICountReport,rep(UniqueAligningUMIs,nOfAlleles))
    
    # Add allele UMI counts (using gene level likelihood calculation) to each allele
    alleleUmiCounts <- sapply(AssignedUMIsList,length)
    uc <- alleleUmiCounts[match(sampleUMICountReport[,1],names(alleleUmiCounts))]
    sampleUMICountReport <- cbind(sampleUMICountReport,uc)
    
    # Add allele UMI counts (using likelihood calculation using all alleles across genes) to each allele
    alleleUmiCounts_allRefs <- sapply(AssignedUMIsList_usingAllReferencekeySites,length)
    uc_allRefs <- alleleUmiCounts_allRefs[match(sampleUMICountReport[,1],names(alleleUmiCounts_allRefs))]
    sampleUMICountReport <- cbind(sampleUMICountReport,uc_allRefs)
    
    # Add allele UMI counts (using unique aligning UMIs only) to each allele
    alleleUmiCounts_uniqueAligning <- sapply(AssignedUMIsList_usingUniqueAligningUMIs,length)
    uc <- alleleUmiCounts_uniqueAligning[match(sampleUMICountReport[,1],names(alleleUmiCounts_uniqueAligning))]
    sampleUMICountReport <- cbind(sampleUMICountReport,uc)
    
    
    
    rownames(sampleUMICountReport) <- AllelesForSample
    colnames(sampleUMICountReport) <- c("Alleles","nReadsInSample","nAlignedReadsToHLA","nAllUniqueUMIs","nCollapsedUMIs","UmiDuplicationRate","nUMIAligningToSingleAllele","nAlleleUMIs_keySitesInGene","nAlleleUMIs_keySitesFromAllGenes","nAlleleUMIs_uniqueAligningUMIs")
    
    write.table(sampleUMICountReport,file=paste(writeToForSample,sid,"_UMI_countingReport.txt",sep=""),row.names=F,sep="\t")
    
    
    colnames(UMIdetails) <- c("UMI","#FoundInRefs","#FoundInAlignedReadSeqs","#FoundInUMIReadSeqs","#nonGGGStartingUMIReads","UMIStatus","AlleleHitsFromAligner","MostLikelyAlleles")
    write.table(UMIdetails,file=paste(writeToForSample,sid,"_UMIdetails.txt",sep=""),row.names=F,sep="\t")
    
    
    #calculateUMIuniquenessPerAllele(UMIsList=AssignedUMIsList)
    
    cat("\t","Total number of unique HLA transcripts (UMIs) evaluated in",sid,":",as.character(length(dedupUMIseqs)),"\n")

    
  }else{ 
    # no alignment found for sample
    cat("\t","No alignment found for sample",sid,".\n")
    
  }
  
  
  ### For each UMI reads, find all alleles to which it aligned. For alleles where it doesnt align,give those alleles likelihood of 0. For all other alleles to which it aligns, calculate allele likelihoods
  
  
  
  
  
  
  
  
  }
  
  
  
}















