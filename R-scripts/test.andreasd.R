## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
##
## An example script, calculating correlation and standard deviations for a list of (local) files
## and save the information to an Rdata file

library(esd)

## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
##
## An example script, colleting metadata from files located on OPeNDAP servers
## and save/add the information to an Rdata file

#Load function
#source('cds.R')

#Define lists with data URLs
# CORDEXList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", 
#                 "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc")
# 
# CORDEXAdjList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/tas/hist/tasAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc", 
#                    "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/pr/hist/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")

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
#Define lists with data URLs
CORDEXList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc")


CORDEXAdjList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/tas/hist/tasAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/pr/hist/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")

#Make a model selection list
ML <- c(CORDEXAdjList[1], CORDEXList[1], CORDEXList[2], CORDEXAdjList[2], CORDEXList[17])

#Function to calculate standard deviation (as defined in taylor.diagram() in the plottrix package)
SD <- function(x, subn) {
  meanx <- mean(x, na.rm = TRUE)
  devx <- x - meanx
  ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) -
                                                subn))
  return(ssd)
}

#Get spatial fields of reference and model data
nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/timmean_tg_EUR11-grid_v14.0.nc")
tref <- ncvar_get(nc,"tg")
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/maskedtimmean_tasAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")
tmod1 <- ncvar_get(nc,"tasAdjust")-273.15
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/maskedtimmean_tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc")
tmod2 <- ncvar_get(nc,"tas")-273.15
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/maskedtimmean_tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc")
tmod3 <- ncvar_get(nc,"tas")-273.15
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/timmean_rr_EUR11-grid_v14.0.nc")
pref <- ncvar_get(nc,"rr")
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/maskedtimmean_prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")
prmod1 <- ncvar_get(nc,"prAdjust")*86400
nc_close(nc)

nc <- nc_open("/run/user/44131/gvfs/sftp:host=xvis-m3b/home/andreasd/scripts_and_eval/ABC4CDE/abc4cde_wp4.git/branches/andreasd/data/maskedtimmean_pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc")
prmod2 <- ncvar_get(nc,"pr")*86400
nc_close(nc)

#Calculate spatial standard deviatons
sdtr <- SD(tref,FALSE)
sdt1 <- SD(tmod1,FALSE)
sdt2 <- SD(tmod2,FALSE)
sdt3 <- SD(tmod3,FALSE)

sdpr <- SD(pref,FALSE)
sdp1 <- SD(prmod1,FALSE)
sdp2 <- SD(prmod2,FALSE)

#Calculate spatial correlations with reference
cort1 <- cor(as.vector(tref),as.vector(tmod1),use = "pairwise")
cort2 <- cor(as.vector(tref),as.vector(tmod2),use = "pairwise")
cort3 <- cor(as.vector(tref),as.vector(tmod3),use = "pairwise")

corp1 <- cor(as.vector(pref),as.vector(prmod1),use = "pairwise")
corp2 <- cor(as.vector(pref),as.vector(prmod2),use = "pairwise")

#Create data fram with URLs, spatial standard deviations (SSD), reference SSD and correlations.
RSD_data <- data.frame(DataURL= ML,
                       SSD= round(c(sdt1,sdt2,sdt3,sdp1,sdp2),3),
                       SSDRef=round(c(sdtr,sdtr,sdtr,sdpr,sdpr),3),
                       CorrWithRef= round(c(cort1,cort2,cort3,corp1,corp2),3))

#Save to RData file
save(RSD_data,file="sd_cor.rda")
