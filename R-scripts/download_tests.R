setwd(system("find $HOME -name download_tests.R -exec dirname {} \\;",intern=T))
cds.path <- system("find $HOME -name cds.R -exec dirname {} \\;",intern=T)
source(cds.path)
require(esd)

#Test for downloading ERA-interim monthly data from the public database
# system("conda activate py27")
# Sys.setenv(PATH = paste("/home/ubuntu/miniconda3/envs/py27/bin/python", Sys.getenv("PATH"), sep=":")) 
# system("python --version")
# test <- getERA("pr",start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",destfile=NULL)
# saveRDS(test,"era.pr.Rdata")
# test <- getERA("tas",start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",destfile=NULL)

#Test for downloading monthly CFSR data from Climate explorer
cfsr.tas <- getCFSR(variable="tas")

#Test for downloading daily E-OBS data and then aggregating in to monthly values.
eobs.pr <- getEOBS(variable="pr")