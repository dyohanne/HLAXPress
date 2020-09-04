
source("../HLAXPress_fxns.R")
ptm <- proc.time()


getAlleleFreqCounts(sampleHlaData="testDataHLAs",
                    sampleId="sample45",IlluminaDir="testData",nanDir=NULL)

print(proc.time() - ptm)


