# Land-Atmosphere Interactions Model (LAIM) 
# By John C. Lin (John.Lin@utah.edu)


#################################################
# Flags to Turn On/Off Processes 
vegcontrolTF <- TRUE    # vegetation control?
atmrespondTF <- TRUE    # does atmosphere respond to surface fluxes?
ABLTF<- TRUE            # does ABL growth or decay, according to surface heat fluxes?
cloudTF <- TRUE         # cloud physics response to relative humidity
soilWTF <- TRUE         # turn on soil moisture feedbacks?
#################################################

#################################################
# Model timestep and duration
dt <- 60           # model timestep [s]
t.day <- 2         # run time in days
tmax <- t.day*24*3600  #maximum time [s]
DTtol <- 0.01      # tolerance for change in T when solving numerically (if T is within this range, then stop iterating) [oK]
countTmax <- 1000  # max number of times to iterate T calculation
#################################################

#################################################
# Land surface characteristics
gvmax <- 1/50      # max vegetation conductance [m/s]; reciprocal of vegetation resistance
albedo.c <- 0.1    # surface albedo
albedo <- albedo.c # surface albedo
z0 <- 0.5          # roughness length [m]
epsilon.s <- 0.97  # surface emissivity for forest, according to Jin & Liang (2006)
LAI <- 3.0         # average leaf area index; for a forest like Harvard Forest, ~3.0 over the year [.]
Kb <- 0.5          # extinction coefficient within plant canopy [.]; average value ~0.5:  https://link.springer.com/article/10.1007/s11707-014-0446-7
# heat capacity of land surface
# a) heat capacity based on soil
D <- 0.1*(1/24)          # the depth of soil that temp fluctuations would penetrate [m]; 0.1m is roughly the depth that would penetrate on diurnal timescales
Cp.soil <- 1921          # specific heat of soil organic material [J/kg/K]
rho.soil <- 1300         # density of soil organic material [kg/m3]
Cs <- Cp.soil*rho.soil*D # heat capacity of organic soil [J/K/m2]
# b) heat capacity based on vegetation
# Hveg <- 10               # height of vegetation [m]
# rho.veg <- 0.7E6         # wood density [g/m3]
# Cp.veg <- 3000           # bulk heat capacity of above-ground vegetation [J/kg/K];  Sect. 7.2 of Bonan (2019)
# Cs <- Cp.veg*(rho.veg/1000)*Hveg/100  # heat capacity of vegetation [J/K/m2]

Lambda <- 5.9       # thermal diffusivity of skin layer [.]
# two-layer (force-restore) soil model, from de Arellano et al. (2015)
Tsoil2 <- 286       # T of deep soil layer [K] that is constant
Tsoil1 <- Tsoil2    # T of top soil layer [K] that varies w/ time
Wsoil2 <- 0.21      # volumetric water concent of deep soil layer [m3/m3]
Wsoil1 <- Wsoil2    # volumetric water concent of top soil layer [m3/m3]
d1 <- 0.1
# soil parameters from Clapp & Hornberger (1978)
# a) Sandy loam soil
Wsat <- 0.472  # saturated volumetric water content [m3/m3]
Wfc <- 0.323   # volumetric water content at field capacity [m3/m3]
Wwilt <- 0.171 # volumetric water content at wilting point [m3/m3]
aa <- 0.219    # Clapp & Hornberger (1978) retention parameter a
bb <- 4.9      # Clapp & Hornberger (1978) retention parameter b
pp <- 4        # Clapp & Hornberger (1978) retention parameter c
CGsat <- 3.56E-6 # saturated soil conductivity for heat [???]
C1sat <- 0.132
C2ref <- 1.8 
#################################################

#################################################
# Physical constants
Cp <- 1005.7;Cv <- 719 # heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel (1994)
g <- 9.80665 # standard surface gravity [m/s2]
Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel (1994))
Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson (1999)
sigma <- 5.670373E-8    # Stefan-Boltzmann constant [W/m2/K4]
Md <- 28.97  #molar mass of dry air [g/mole]
rho.W <- 1000 # density of water [kg/m3]
#################################################

#################################################
# ---------- External forcing ------------------#
# Downward shortwave radiation
t.hr<-0:24
# a) hourly varying
SWdn<--15*(t.hr-12)^2+800 # hourly downward shortwave radiation [W/m2]
names(SWdn)<-t.hr
SWdn[as.character(c(0:5,19:24))]<-0
# b) constant
# SWdn[1:length(SWdn)]<-1000

# Downward longwave radiation
LWdn <- SWdn; LWdn[1:length(LWdn)] <- 350 # constant downward longwave radiation [W/m2]

# Air temperature
Ta.c<- -0.5*(t.hr-12)^2+30  # PRESCRIBED air temperature [deg-C]
names(Ta.c) <- t.hr
Ta.c[1:length(Ta.c)] <- 5    # override with CONSTANT air temperature [deg-C]

# specific humidity of air:  determine from RH, air temperature
#qair<-2/1000      #specific humidity of air [g/g]
#JHD
#  RH<-1.0     	   #Edited this to force precipitation
RH<-0.8
e<-RH*satvap(mean(Ta.c))/100  #vapor pressure [hPa]
Psurf<-1000     #surface pressure [hPa] 
qair<-(Rd/Rv)*e/Psurf   #specific humidity [g/g]

# atmospheric conditions
hmin <- 100       # minimum height of atmospheric boundary layer [m]
thetavM0<-(Ta.c[1]+273.15)*(1+0.61*qair) # initial virtual potential temperature [K]; Eq. 1.5.1b of Stull [1988]
Beta <- 0.2       # closure hypothesis:  fraction of surface virtual potential temperature flux that determines entrainment heat flux
gamma <- 5/1000   # slope of thetav above growing ABL [K/m]
qabove <- qair/5  # specific humidity of air above ABL [g/g] changed from 5 to 1
W <- 0            # subsidence rate [m/s]
Ur <- 1           # reference windspeed [m/s] at reference height zr
zr <- 50          # reference height [m] where Ur applies
Cair <- 400       # atmospheric CO2 concentration [umole/mole, or ppm]
#################################################

#################################################
# --------------Functions-----------------#
if(vegcontrolTF)source("Ball_Berry_Farquhar.r")  #load Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model  

latentheat <- function(T.c){
  # Takes temperature [C] and returns value of latent heat of vaporization [J/g]
  lambda <- 2.501 - 0.0024 * T.c
  return(lambda * 1000)
} #latentheat<-function(T.c){

satvap <- function(T.c){
  # Takes temp in Celsius as argument and returns saturation vapor pressure [Pa]
  # NOTE:  calls upon function 'latentheat'
  # 3/19/1998
  kelvin <- T.c+273.15
  #return value in [J/g], so multiply by 1000 to convert into [J/kg]
  lambda <- 1000*latentheat(T.c)
  #461 is Gas Constant for water vapor [J deg-1 kg-1]
  #611 & 273.15 are reference vapor pressure and reference temp., respectively
  saturated <- 611*exp((lambda/461)*((1/273.15)-(1/kelvin)))
  return(saturated)
} #satvap<-function(T.c){

# function to calculate vegetation resistance using Jarvis-type empirical model
rv.f <- function(T,VPD,gvmax,Tmin.c=0,Tmax.c=60,Topt.c=30,VPDmin=1500,VPDmax=7500){
  # Jarvis [1976] type model for stomatal control
  # 1) temperature function
  T.c <- T-273.15  
  fT <- (T.c-Tmin.c)*(T.c-Tmax.c)/((T.c-Tmin.c)*(T.c-Tmax.c)-(T.c-Topt.c)^2)
  fT[fT<0] <- 0
  # 2) vapor pressure deficit (VPD) function
  fVPD <- rep(1.0,length(VPD))
  fVPD[VPD>VPDmax] <- 0      # hi VPD--stomata completely closed
  sel <- VPD>=VPDmin & VPD<=VPDmax
  slope<- -1/(VPDmax-VPDmin)
  fVPD[sel] <- 1+slope*(VPD[sel]-VPDmin)
  
  # final penalty function
  beta <- fT*fVPD
  gv <- beta*gvmax
  rv <- 1/gv   # vegetation resistance [s/m]
  rv[rv>1E6] <- 999999   # when infinite, just assign a really large number
  return(rv)
} #rv.f<-function()

rvs <- rv.f(T=273.15+20,VPD=0:4000,gvmax=gvmax)

# function to calculate aerodynamic resistance
ra.f <- function(Ur=1,zr=50,z0=z0,d=0,rho=1){
  #arguments:  Ur is reference windspeed [m/s] at reference height zr
  #            zr is reference height [m] where Ur applies
  #            z0 is roughness length [m]; 0.01m is typical value for crop
  #            d is displacement height [m]
  #            rho is air density [kg/m^3]
  k <- 0.4  # von Karman constant

  CD <- (k^2)/(log((zr-d)/z0))^2  # aerodynamic transfer coefficient
  ra <- 1/(CD*Ur)                 # aerodynamic resistance [s/m]
  return(ra)
} #ra.f<-function(){
#################################################

# initialize T with equilibrium value (determined through "uniroot")
f<-function(T,Ta,SWdn,LWdn,albedo,epsilon.s,Tsoil1,Ur,zr,z0,gvmax=gvmax){
  # --------------Physical constants--------#
  Cp <- 1005.7; Cv <- 719 # heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel [1994])
  g <- 9.80665 # standard surface gravity [m/s2]
  Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel [1994])
  Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson [1999])
  sigma <- 5.670373E-8 # Stefan-Boltzmann constant [W/m2/K4]
  # --------------Physical constants--------#
  LWup <- epsilon.s*sigma*T^4
  SWup <- albedo*SWdn
  Rnet <- SWdn-SWup+LWdn-LWup
  
  # determine sensible heat flux
  rho <- 1   # air density [kg/m3]
  ra <- ra.f(Ur=Ur,zr=zr,z0=z0,rho=rho)
  H <- (Cp*rho/(ra))*(T-Ta)   # [W/m2]

  # determine latent heat flux
  lambda <- 1000*latentheat(T-273.15)  # latent heat of vaporization [J/kg]
  esat <- satvap(T-273.15)/100 # saturation specific humidity [hPa]
  e <- qair*Psurf/(Rd/Rv)      # vapor pressure [hPa]
  VPD <- 100*(esat-e)          # vapor pressure deficit [Pa]
  qstar <- (Rd/Rv)*esat/Psurf  # saturation specific humidity [g/g]
  if (vegcontrolTF) {
    # a) Jarvis-type empirical model for vegetation resistance [s/m]
    #rv <- rv.f(T=T,VPD=VPD,gvmax=gvmax)   
    # b) Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model for vegetation resistance [s/m]
    hs <- e/esat  # fractional humidity (=1/RH) at leaf surface [.]
    cs <- Cair  # CO2 concentration at leaf surface [umole/mole]
    BBFout <- BBF(SW=SWdn.t,Tleaf.C=T-273.15,hs=hs,beta=1.0,cs=cs,Psurf=Psurf)  
    gsw <- BBFout["gsw"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]  
    rho.mole <- rho*1000/Md # air density [kg/m3] => molar density [moles/m3]
    gsw <- gsw/rho.mole   # [mole/m2/s] => [m/s]
    rv <- 1/gsw           # vegetation resistance [s/m]
  } else {
    rv <- 1/gvmax
  } # if(vegcontrolTF){
  # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
  scale.canopy<-(1-exp(-Kb*LAI))/Kb
  An <- An*scale.canopy
  gv <- (1/rv)*scale.canopy
  rv <- 1/gv
  LE <- (lambda*rho/(ra+rv))*(qstar-qair) #[W/m2]

  # determine ground heat flux from two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
  G <- Lambda * (T - Tsoil1)
  
  # this should =0 when T is at equilibrium value
  return(Rnet-H-LE-G)
} # f<-function(T,Ta,SWdn,LWdn,albedo,epsilon.s){

xinterv <- Ta.c[1]+273.15+c(-50,50)  # interval over which to search for equil temperature
# use initial radiation, temps to solve for initial equil. temperature
Tinit <- uniroot(f,interval=xinterv,Ta=Ta.c[1]+273.15,SWdn=SWdn[1],LWdn=LWdn[1],Tsoil1=Tsoil1,
                 albedo=albedo,epsilon.s=epsilon.s,Ur=Ur,zr=zr,z0=z0,gvmax=gvmax)$root
tmp <- f(T=Tinit,Ta=Ta.c[1]+273.15,SWdn=SWdn[1],LWdn=LWdn[1],Tsoil1=Tsoil1,
         albedo=albedo,epsilon.s=epsilon.s,Ur=Ur,zr=zr,z0=z0,gvmax=gvmax)
print(paste("Tinit [oC]:",signif(Tinit-273.15,5),"; ",signif(tmp,5)))
# Impose perturbation
# Tinit <- Tinit+10

# ---------------------------------------Time Loop START------------------------------------------#
tcurr <- 0
T <- Tinit
result_s <- NULL
h <- hmin    # initialize with minimum 
thetaM <- Ta.c[1]+273.15
qM <- qair   # initialize with prescribed specific humidity [g/g]
thetavM <- thetaM*(1+0.61*qM)   # virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]
countp <- 0
tcc <- 0   # total cloud fraction
while (tcurr<tmax) {
  countp <- countp+1
  if (countp>500) {
    print(paste("-------- tcurr=",round(tcurr/(3600*24)*10)/10,"[days] --------"))
    countp <- 0
  } #if (countp>500)
  
  LWup <- epsilon.s*sigma*T^4   # upward longwave radiation [W/m2]
  LWdn.t <- approx(x=as.numeric(names(LWdn))*3600,y=LWdn,xout=tcurr%%(24*3600))$y  # downward shortwave radiation [W/m2]
  if (cloudTF) {
	  RHe <- 0.5			# fitting parameter sets TCC to gain at approximately 60% humidity
	  tcc <- min(exp((qM/qstar-1)/(1-RHe)),1) 	# Walcek 1994 model equation (1) 	
	  LWdn.t<- LWup-100+80*(tcc)
	  albedo <- albedo.c+0.75*(tcc)
	} #if (cloudTF)
  SWdn.t <- approx(x=as.numeric(names(SWdn))*3600,y=SWdn,xout=tcurr%%(24*3600))$y  # downward shortwave radiation [W/m2]
  SWup <- albedo*SWdn.t
  # determine net radiation
  Rnet <- SWdn.t-SWup+LWdn.t-LWup

  countT <- 0; iterateT <- TRUE
while (iterateT) {   #iterate until convergence
  countT <- countT + 1
  if(countT > countTmax)stop("T does not converge")
    
  if (atmrespondTF) {
    Ta<-thetaM   #air temperature [K] is the one from zero-order jump model
    qa<-qM
  } else {
    Ta <- approx(x=as.numeric(names(Ta.c))*3600,y=Ta.c,xout=tcurr%%(24*3600))$y+273.15  #use prescribed value
    qa <- qair
  } # if(atmrespondTF){
  # determine sensible heat flux
  rho <- 1   # air density [kg/m3]
  ra <- ra.f(Ur=Ur,zr=zr,z0=z0,rho=rho)
  H <- (Cp*rho/(ra))*(T-Ta)   # [W/m2]

  # determine latent heat flux
  beta <- 1   # water stress parameter (dependent on soil moisture)
  lambda <- 1000*latentheat(T-273.15)  # latent heat of vaporization [J/kg]
  esat <- satvap(T-273.15)/100 # saturation specific humidity [hPa]
  e <- qa*Psurf/(Rd/Rv)        # vapor pressure [hPa]
  VPD <- 100*(esat-e)            # vapor pressure deficit [Pa]
  qstar <- (Rd/Rv)*esat/Psurf    # saturation specific humidity [g/g]
  if (vegcontrolTF) {
    if (soilWTF) {
      # Eq. (12.56) of Bonan (2019)
      beta <- (Wsoil1 - Wwilt)/(Wfc - Wwilt)
      if (Wsoil1 >= Wfc) beta <- 1.0
      if (Wsoil1 <= Wwilt) beta <- 0
      # print(paste("Wsoil1, beta:",round(Wsoil1,4),round(beta,4)))
    } # if (soilWTF)
    # a) Jarvis-type empirical model for vegetation resistance [s/m]
    # rv <- rv.f(T=T,VPD=VPD,gvmax=gvmax)
    # b) Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model for vegetation resistance [s/m]
    hs <- e/esat  # fractional humidity (=1/RH) at leaf surface [.]
    cs <- Cair  # CO2 concentration at leaf surface [umole/mole]
    BBFout <- BBF(SW=SWdn.t,Tleaf.C=T-273.15,hs=hs,beta=beta,cs=cs,Psurf=Psurf)  
    gsw <- BBFout["gsw"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]  
    rho.mole <- rho*1000/Md # air density [kg/m3] => molar density [moles/m3]
    gsw <- gsw/rho.mole   # [mole/m2/s] => [m/s]
    rv <- 1/gsw        # vegetation resistance [s/m]
    An <- BBFout["An"]    # Net photosynthesis [umole/m2/s]
    ci <- BBFout["ci"]    # intercellular CO2 [umole/mole]
  } else {
    rv <- 1/gvmax
    An <- NA; ci <- NA
  } # if(vegcontrolTF){
  
  # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
  scale.canopy<-(1-exp(-Kb*LAI))/Kb
  An <- An*scale.canopy
  gv <- (1/rv)*scale.canopy
  rv <- 1/gv
  LE <- (lambda*rho/(ra+rv))*(qstar-qa) #[W/m2]
  
  # determine ground heat flux 
  # use two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
  G <- Lambda * (T - Tsoil1)
  
  Storage <- Rnet - LE - H - G
  # update temperature 
  DT <- (Storage/Cs)*dt
  T <- T+DT
  #print(paste("iterating so that T converges:",countT,paste("T =",signif(T,5)),signif(DT,4)))
  iterateT <- abs(DT)>DTtol   # continue iterating until T converges
} # while (iterateT) {   #iterate until converge

  # heat transport between surface and deep soil layer to update Tsoil1 from CLASS model
  CG <- CGsat * (Wsat/Wsoil2)^(bb/(2*log(10)))
  Tsoil1 <- Tsoil1 + (CG*G - (2*pi/(86400))*(Tsoil1 - Tsoil2))*dt  #Eq. (9.32) of de Arellano et al. (2015)
  
  if (soilWTF) {
    # update soil water content, based on CLASS model
    C1 <- C1sat*(Wsat/Wsoil1)^(bb/2 + 1)  #Eq. (9.35) of de Arellano et al. (2015)
    Wsmall <- 1E-3
    C2 <- C2ref*(Wsoil2/(Wsat - Wsoil2 + Wsmall))  #Eq. (9.36) of de Arellano et al. (2015)
    Wsoil1eq <- Wsoil2 - aa*Wsat*((Wsoil2/Wsat)^pp)*(1-(Wsoil2/Wsat)^(8*pp))  #Eq. (9.37) of de Arellano et al. (2015)
    # Eq. (9.34) of de Arellano et al. (2015); NOTE:  use LE instead of LEsoil as in (9.34), and -1 multiplied by C1 that is missing in (9.34)
    dWsoil <- ((-C1/(rho.W*d1))*(LE/lambda) - (C2/86400)*(Wsoil1 - Wsoil2))*dt  
    dWsoil.a <- ((-C1/(rho.W*d1))*(LE/lambda))*dt
    dWsoil.b <- (- (C2/86400)*(Wsoil1 - Wsoil2))*dt
    Wsoil1 <- Wsoil1 + dWsoil
    if (Wsoil1 < 0) Wsoil1 <- 0
    # print(paste("Wsoil:",round(Wsoil1,4),signif(C1,4),signif(C2,4),round(Wsoil1eq,4),round(dWsoil,4),signif(dWsoil.a,4),signif(dWsoil.b,4)))
  } #if (soilWTF) {
  
  # if want atmosphere to respond
  # Based on "zero-order jump" or "slab" model of convective boundary layer, described in Pg. 151~155 of Garratt [1992]
  if (atmrespondTF) {
    #update ABL-averaged q
    #calculate surface virtual heat flux
    lambda <- latentheat(T-273.15)  # latent heat of vaporization [J/g]
    E <- LE/lambda   # surface moisture flux [g/m^2/s] 
    F0theta <- H/Cp  # potential heat flux [K-kg/m^2/s]
    F0thetav <- F0theta+0.073*lambda*E/Cp # virtual heat flux [K-kg/m^2/s]
    if (ABLTF) {
      Fhthetav <- -1*Beta*F0thetav   # closure hypothesis (Eq. 6.15 of Garratt [1992])
      # calculate ABL growth rate
      dh.dt<-(1+2*Beta)*F0thetav/(gamma*h)
      if (F0thetav<=0.00) dh.dt <- 0 # ABL collapses
    } else {
      dh.dt <- 0
      Fhthetav <- 0
    } # if(ABLTF){
    
    # calculate entrainment flux of humidity
    deltaq <- qabove - qM
    Fhq <- -1*rho*deltaq*(dh.dt-W)*1000  # entrainment flux of humidity [g/m2/s] NOTE:  assume CONSTANT air density!
    dq.dt <- (E-Fhq)/(rho*1000*h) # change of humidity in ABL [1/s]; NOTE:  assume CONSTANT air density!
    qM <- qM+dq.dt*dt             # updated qM [g/g]
    # update ABL-averaged thetav
    dthetavM.dt <- (F0thetav-Fhthetav)/h  # change of thetav in ABL [K-kg/m^3/s]
    dthetavM.dt <- dthetavM.dt/rho        # [K-kg/m^3/s]=>[K/s]
    thetavM <- thetavM+dthetavM.dt*dt
    thetaM <- thetavM/(1+0.61*qM)   # potential temperature, from virtual pot temp [K];  Eq. 1.5.1b of Stull [1988]
    # update ABL height
    h <- h+dh.dt*dt
    if (F0thetav<=0.00&ABLTF) h<-hmin # override value:  ABL collapses
  } # if(atmrespondTF){
  
  tmp <- c(tcurr,T,Ta,Tsoil1,Wsoil1,beta,LWup,Rnet,H,LE,G,DT,qstar,qa,ra,rv,An,ci,h,qM,thetavM,thetaM)
  result_s <- rbind(result_s,tmp)
  tcurr <- tcurr + min(c(dt,tmax-tcurr))
} # while(tcurr<ttmax){
#---------------------------------------Time Loop END ------------------------------------------#
result <- result_s[-1,]   # remove initial value
dimnames(result) <- list(NULL,c("time","T","Ta","Tsoil1","Wsoil1","beta","LWup","Rnet","H","LE","G","DT",
                              "qstar","qair","ra","rv","An","ci","h","qM","thetavM","thetaM"))
filenm <- "result.csv"
write.csv(result,file=filenm)
print(paste(filenm,"written out"))

# text on plot 
xmain <- paste("vegcontrolTF=",vegcontrolTF)
xmain <- paste(xmain,"  atmrespondTF=",atmrespondTF)
xmain <- paste(xmain,"\nABLTF=",ABLTF)
xmain <- paste(xmain,"  soilWTF=",soilWTF)
xmain <- paste(xmain,"  cloudTF=",cloudTF)
# regenerate VPD from qstar and qair 
e <- result[,"qair"]*Psurf/(Rd/Rv)      # vapor pressure [hPa]
esat <- result[,"qstar"]*Psurf*(Rv/Rd)  # saturation vapor pressure [hPa]
VPD <- 100*(esat-e)                     # vapor pressure deficit [Pa]

#################################################
# -------------------- Plotting ----------------#

# generate 4 different plots in one window, with 2*2 configuration
dev.new(); par(mfrow=c(2,2),cex.main=0.7)   
ylims <- range(result[,c("T","Ta","Tsoil1")]-273.15)
plot(result[,"time"]/3600,result[,"T"]-273.15,type="l",xlab="Time [hour]",ylab="Temperature [deg-C]",
     cex.axis=1.3,cex.lab=1.3,lwd=2,ylim=ylims,main=xmain)
lines(result[,"time"]/3600,result[,"Ta"]-273.15,lty=2,lwd=2)
lines(result[,"time"]/3600,result[,"Tsoil1"]-273.15,lty=3,lwd=2)
legend(x="topright",c("Tsurf","Tair","Tsoil"),lwd=2,lty=c(1,2,3))

plot(result[,"time"]/3600,VPD,type="l",xlab="Time [hour]",ylab="VPD [Pa]",lwd=2,
     cex.axis=1.3,cex.lab=1.3,main=xmain)

ylims <- range(result[,c("qstar","qair")]*1000)
plot(result[,"time"]/3600,result[,"qstar"]*1000,type="l",xlab="Time [hour]",ylab="Specific humidity [g/kg]",
     cex.axis=1.3,cex.lab=1.3,lwd=2,ylim=ylims,main=xmain)
lines(result[,"time"]/3600,result[,"qair"]*1000,lty=2,lwd=2)
legend(x="topright",c("qstar","qair"),lwd=2,lty=c(1,2))

ylims <- range(result[,c("ra","rv")])
plot(result[,"time"]/3600,result[,"ra"],type="l",xlab="Time [hour]",ylab="Resistances [s/m]",
     cex.axis=1.3,cex.lab=1.3,lwd=2,ylim=ylims,main=xmain,col="yellow")
lines(result[,"time"]/3600,result[,"rv"],lty=1,lwd=2,col="lightgreen")
legend(x="topright",c("raero","rveg"),lwd=2,lty=c(1),col=c("yellow","lightgreen"))
dev.copy(png,"T_q_r_VPD.png");dev.off()

# plot with energy fluxes 
dev.new()
matplot(result[,"time"]/3600,result[,c("Rnet","LWup","H","LE","G")],type="l",lty=c(1,2,1,1,1),
        cex.axis=1.5,cex.lab=1.5,col=c("black","black","orange","blue","darkgreen"),lwd=2,xlab="Time [hr]",ylab="Energy Fluxes [W/m2]")
legend(x="topright",c("Rnet","LWup","H","LE","G"),col=c("black","black","orange","blue","darkgreen"),lwd=2,lty=c(1,2,1,1,1))
title(main=xmain)
dev.copy(png,"Energyfluxes.png");dev.off()

# plot soil water content 
if (soilWTF) {
  dev.new()
  plot(result[,"time"]/3600,result[,"Wsoil1"],type="l",xlab="Time [hour]",ylab="Soil Volumetric Water Content [m3/m3]",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain,col="black",ylim=c(0,Wfc))
  abline(h=Wsoil2,lty=3,lwd=2)
  lines(result[,"time"]/3600,result[,"beta"],type="l",col="lightgreen",lwd=2)
  legend(x="topright",c("Wsoil1","Wsoil2","beta"),col=c("black","black","lightgreen"),lwd=2,lty=c(1,3,1))
  dev.copy(png,"Wsoil.png");dev.off()
} #if (soilWTF)

if (atmrespondTF) {
  # plot time series of ABL height 
  dev.new()
  plot(result[,"time"]/3600,result[,"h"],type="l",xlab="Time [hour]",ylab="ABL height [m]",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain)
  dev.copy(png,"ABLht.png");dev.off()
} #if(atmrespondTF){

if (vegcontrolTF) {
  # plot time series of photosynthetic uptake 
  dev.new()
  plot(result[,"time"]/3600,result[,"An"],type="l",xlab="Time [hour]",ylab="Net Photosynthesis [umole/m2/s]",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain)
  dev.copy(png,"PSN.png");dev.off()
} #if(vegcontrolTF){
#################################################
