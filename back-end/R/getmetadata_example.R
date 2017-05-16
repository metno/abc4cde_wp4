## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.doblerd@met.no  Oslo, Norway, 2017-05-16
##
## An example script, colleting metadata from files located on OPeNDAP servers
## and save/add the information to an Rdata file

#Load function
source('cds.R')

#Define lists with data URLs
CORDEXList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc")

CORDEXAdjList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/tas/hist/tasAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/pr/hist/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")

###########
#Check if metadata file exists
#If yes: load
#If no: generate an empty list and get first entry
if (file.exists("metaextracted.rda"))
{
  load("metaextracted.rda")
} else {
  mdList <- list()
  dataurl <- CORDEXList[1]
  mdList[[length(mdList)+1]] <- metaextract_opendap(dataurl)
}

#Get the metadata from the first list
#Skip the entries with matching URLs
#and add the new ones
for (dataurl in CORDEXList)
{
  if (sum(sapply(mdList,function(X) X$dataurl == dataurl)) ==0 ){
    print(paste("get metadata from",dataurl))
    mdList[[length(mdList)+1]] <- metaextract_opendap(dataurl)
  }
}

#Add the the metadata from the second list
#Skip the entries with matching URLs
#and add the new ones
for (dataurl in CORDEXAdjList)
{
  if (sum(sapply(mdList,function(X) X$dataurl == dataurl)) ==0 ){
    print(paste("get metadata from",dataurl))
    mdList[[length(mdList)+1]] <- metaextract_opendap(dataurl)
  }
}
#Feedback message
print("All metadata here!")

#Update the metadata file
save(mdList,file="metaextracted.rda")

####Example query & output
#Do an example query and print results
query <- which(sapply(mdList,function(X) ("tas" %in% names(X$varlist) || "pr" %in% names(X$varlist)) && (X$globat$model_id == "SMHI-RCA4")))

print("Your query matches:")
for (q in query)
{
  print(mdList[[q]]$dataurl)
  # print(mdList[[q]]$globat$driving_model_id)
}
#Print an overview
# print("###ALL###")
# for (q in 1:length(mdList))
# {
#   print(mdList[[q]]$globat$driving_model_id)
#   print(mdList[[q]]$dimlist$time$calendar)
# }

