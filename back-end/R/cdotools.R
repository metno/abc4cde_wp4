# Calculate the mean value over time with CDO. Faster than in R.
# Is this function affected by the error that Abdelkader discovered with -timmean?

cdo.mean <- function(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE,
                     is.temp=TRUE,verbose=FALSE) {
  
  commands <- c("-fldmean","-timmean","-selyear")
  input <- c("","",paste(period,collapse="/"))
  
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timmean","-yseasmean")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model.file,out.file)
  
  command <- ("output")
  input <- c("")
  
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=TRUE))
  if(seasonal) {
    names(out) <- c("djf","mam","jja","son")
  } else {
    names(out) <- "ann"
  } 
  # If applying to e.g. slp data, set is.temp to FALSE to skip this correction:
  if(out>200 & is.temp) out <- out-273.15 
  system(paste("rm",out.file,sep=" "))
  invisible(out)
}

## Calculate the temporal standard deviation with cdo
cdo.timeSd <- function(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE) {
  
  commands <- c("-timstd","-fldmean","-ymean","-selyear")
  input <- c("","",paste(period,collapse="/"))
  
  if(!is.null(mask)) {
    commands <- append(commands,"-maskregion",after=3)
    input <- append(input,mask,after=3) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-ymean","-seasmean")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model.file,out.file)
  
  command <- ("output")
  input <- c("")
  
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=TRUE))
  if(seasonal) {
    names(out) <- c("djf","mam","jja","son")
  } else {
    names(out) <- "ann"
  }
  system(paste("rm",out.file,sep=" "))
  invisible(out)
}

# Calculate the spatial standard deviation with cdo
cdo.spatSd <- function(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE) {
  
  commands <- c("-fldstd","-timmean","-selyear")
  input <- c("","",paste(period,collapse="/"))
  
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timmean","-yseasmean")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model.file,out.file)
  
  command <- ("output")
  input <- c("")
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=TRUE))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  }
  system("rm tmp.nc")
  invisible(out)
}

# Calculate the spatial correlation of two gridded data sets with cdo
cdo.gridcor <- function(model.file,reference.file,period=c(1981,2010),mask=NULL,seasonal=FALSE) {
  
  commands <- c("-timavg","-selyear")
  input <- c("",paste(period,collapse="/"))
  
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timavg","-yseasavg")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model.file,out.file)
  
  out.file <- "tmp2.nc"
  cdo.command(commands,input,reference.file,out.file)
  
  commands <- c("fldcor")
  input <- c("")
  in.file <- c("tmp.nc tmp2.nc")
  out.file <- "tmp_cor.nc"
  cdo.command(commands,input,in.file,out.file)
  
  command <- ("output")
  input <- c("")
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=TRUE))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  }
  system("rm tmp.nc tmp2.nc tmp_cor.nc")
  invisible(out)
}


#Apply a set of cdo commands on a grib/netcdf file. Several commands can be piped.
cdo.command <- function(commands,input,infile,outfile,intern=FALSE) {
  cdo.coms <- array()
  separators <- array(" ",dim=length(commands))
  separators[which(is.na(match(input,"")))] <- ","
  for(i in 1:length(separators)){
    cdo.coms[i]  <- paste(commands[i],input[i],sep=separators[i])
  }
  system.command <- paste("cdo",paste(cdo.coms,collapse=" "),infile,outfile,sep=" ")
  
  if(intern) {
    output <- system(system.command,wait=TRUE,intern=TRUE)
    return(output)
  } else {
    system(system.command,wait=TRUE)
  }
}

#Unzip a gz package
gunzip <- function(filename) {
  system.command <- paste("gunzip",filename)
  system(system.command,wait=TRUE)
}
