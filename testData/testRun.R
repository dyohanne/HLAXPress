
source("../HLAXPress_fxns.R")
ptm <- proc.time()


getAlleleFreqCounts(sampleHlaData="/scratch/project_2000416/dawit/TiiraOutputFiles/HLAXpressAlleleSpecExp_changesAfterReview/ReRuns_2020/HLAxpress_forGithub/testData/testDataHLAs",
                    sampleId="sample45",IlluminaDir="/scratch/project_2000416/dawit/TiiraOutputFiles/HLAXpressAlleleSpecExp_changesAfterReview/ReRuns_2020/HLAxpress_forGithub/testData",nanDir=NULL)

print(proc.time() - ptm)


