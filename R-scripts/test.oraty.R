## test.ratyo.R
## Test Olle's back-end functions

file.meta <- find.file("metaextract.rda")
file.meta <- file.meta[grep("abc4cde",file.meta)]

file.shape <- find.file("*.shp")
file.shape <- file.shape[grep("abc4cde",file.shape)]
X <- get.shapefile(file.shape,with.path=TRUE)

# missing script getMonthlyERA.py
#getERA("tas",start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",destfile="test.nc")
#python.getEra(1979,2016,"tas",0,type,stream,outfile)

getCFSR()