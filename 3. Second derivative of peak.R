################################################################################
########################## Notes on provided code ##############################
################################################################################

# The code provide is for determining the wavenumber range for carbon-deuterium 
# peak in the silent region of Raman spectra by applying a spline smoothing 
# algorithm and calculating the second derivative for the smoothed line. The
# input data should be the processed data output from the R code 2 file found 
# on the same Github page. The code will calculate the second derivative for all
# spectra and provide the average wavenumbers to use as a range.

################################################################################
rm(list=ls()); graphics.off()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

palette('R3')

pkgs=c('dplyr','doBy','readxl','pspline')

  for (p in pkgs) {
    if(!require(p,character.only = TRUE)) {
      install.packages(p)
      require(p,character.only = TRUE)}
  }

###################### Functions for following code ############################

# Determine the operating system
get_os = function(){
  sysinf = Sys.info()
  if (!is.null(sysinf)){
    os = sysinf['sysname']
    if (os == 'Darwin')
      os = "osx"
  } else {
    os = .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os = "osx"
    if (grepl("linux-gnu", R.version$os))
      os = "linux"
  }
  tolower(os)
}

# Plot function
plot_x11 = function(display="",width=7,height=7,pointsize=12,
                gamma=1,bg="transparent",canvas="white",xpos=NA,ypos=NA,
                title=""){
  
  (opsys = as.character(get_os()))
  if(opsys!="osx"){
    X11(display=display,width=width,height=height,
        pointsize=pointsize,gamma=gamma,bg=bg,canvas=canvas,
        xpos=xpos,ypos=ypos,title=title)
  }
  else {
    quartz(title=title,width=width,height=height,
           pointsize=pointsize,bg=bg,canvas=canvas)
  }
}

################################################################################

# Load Raman spectra that have been preprocessed
Raman_spectra = as.data.frame(read.csv('Raman_SNIP_normalized.csv'))
head(Raman_spectra)
names(Raman_spectra)

################################################################################

# Quick visualization of data
x11()
plot(Raman_spectra$wavenumber,Raman_spectra$normalized,pch=20,cex=0.5)

# Truncate to region of interest in spectra
CD_peak=c(2000,2400)

abline(v=CD_peak,col=2,lwd=3)

################################################################################

# Create new column for unique ID
Raman_spectra$cell=Raman_spectra$FileID
(cells=(unique(Raman_spectra$cell)))

################################################################################

DATA=list(NULL)
j=1
for(j in 1:length(cells)){
graphics.off()
(celln=cells[j])

# Pull celln data
tmp_data=subset(Raman_spectra,cell==celln)
dim(tmp_data)

# Process first peak data
tmp_data=subset(tmp_data,wavenumber>=CD_peak[1]&wavenumber<=CD_peak[2])

head(tmp_data);tail(tmp_data)
dim(tmp_data)

# Change df value to set the degrees of freedome (DOF) for pspline function
smobj5=sm.spline(tmp_data$wavenumber,tmp_data$normalized, df=5) 

plot_x11()
plot(tmp_data$wavenumber,tmp_data$normalized,pch=20,cex=0.75,col=1, main='x1')
points(smobj5$x,smobj5$ysmth,type='l',col=1,lwd=3)

(pick0=which(smobj5$ysmth==max(smobj5$ysmth)))
(x0=smobj5$x[pick0])
abline(v=x0)

x=tmp_data$wavenumber

################################

# Change as needed for changing df in pspline
smobj=smobj5

################################

d2=predict(smobj,x,nderiv=2)
class(d2)
dim(d2)
x11()
plot(x,d2,type='l',lwd=6,main='x1 smobj5 deriv 2')
abline(v=x0)

pick1=which(d2==max(d2[x<x0]))
x2L=x[pick1]
abline(v=x2L,lty=2)

pick2=which(d2==max(d2[x>x0]))
x2R=x[pick2]
abline(v=x2R,lty=2)

c(x2L,x2R)

plot_x11()
plot(tmp_data$wavenumber,tmp_data$normalized,pch=20,cex=0.75,col=1, 
     main='smobj5 deriv 2')
points(smobj5$x,smobj5$ysmth,type='l',col=1,lwd=3)

abline(v=x0)
abline(v=c(x2L,x2R),lty=2,lwd=2,col=2)

P1D2L=x2L; P1D2R=x2R

(dfx=as.data.frame(cbind(P1D2L,P1D2R)))
dfx$j=j; dfx$celln=celln
DATA[[j]]=dfx

print(paste('loop',j,'of',length(cells)))
} 

###############################################################################

dfX=bind_rows(DATA)
summary(dfX)
(tmp=dfX[,c('P1D2L','P1D2R')])

(dfx=apply(tmp,2,mean))

x11()
(ylim=range(Raman_spectra$normalized))
plot(Raman_spectra$wavenumber,Raman_spectra$normalized,pch=20,cex=0.25,ylim=ylim)

abline(v=as.vector(dfx),lty=2,col=rep(c(2,3,4,6),each=2),lwd=2)

