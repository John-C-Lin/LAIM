# Land-Atmosphere Interactions Model (LAIM) 
# By John C. Lin (John.Lin@utah.edu)

require("deSolve")   #load deSolve package to access function "ode"

#################################################
# Flags to Turn On/Off Processes 
atmrespondTF <- TRUE    # does atmosphere respond to surface fluxes?
ABLTF <- TRUE           # does ABL grow or decay, according to surface heat fluxes?
cloudTF <- TRUE         # does cloud cover change as function of atmospheric humidity?
vegcontrolTF <- TRUE    # vegetation control?
soilWTF <- TRUE         # turn on soil moisture feedbacks?
co2budgetTF <- TRUE     # track atmospheric CO2, based on surface and entrainment fluxes? 
if (!atmrespondTF & ABLTF) stop ("atmrespondTF needs to be TRUE to allow ABL to grow and decay")
if (!vegcontrolTF & co2budgetTF) stop ("vegcontrolTF needs to be TRUE to track CO2")
if (!vegcontrolTF & soilWTF) stop ("vegcontrolTF needs to be TRUE for soil moisture feedback to work")
LWdnTF <- TRUE          # does LWdn respond dynamically?  
#################################################

#################################################
# Model timestep and duration
dt <- 20           # model timestep [s]
t.day <- 2         # run time in days
tmax <- t.day*24*3600  #maximum time [s]
times <- seq(0,tmax,dt) #vector of time steps [s]
DTtol <- 0.01      # tolerance for change in T when solving numerically (if T is within this range, then stop iterating) [oK]
countTmax <- 1000  # max number of times to iterate T calculation
#################################################

#################################################
# Load in functions
if(vegcontrolTF){
  if(!file.exists("Ball_Berry_Farquhar.r"))stop(paste("Can not find 'Ball_Berry_Farquhar.r' in working directory:",getwd()))
  source("Ball_Berry_Farquhar.r")  #load Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model  
} # if(vegcontrolTF){

latentheat <- function(T.c){
  # Takes temperature [C] and returns value of latent heat of vaporization [J/g]
  Lv <- 2.501 - 0.0024 * T.c
  return(Lv * 1000)
} #latentheat<-function(T.c){

satvap <- function(T.c){
  # Takes temp in Celsius as argument and returns saturation vapor pressure [Pa]
  # NOTE:  calls upon function 'latentheat'
  # 3/19/1998
  kelvin <- T.c+273.15
  #return value in [J/g], so multiply by 1000 to convert into [J/kg]
  Lv <- 1000*latentheat(T.c)
  #461 is Gas Constant for water vapor [J deg-1 kg-1]
  #611 & 273.15 are reference vapor pressure and reference temp., respectively
  saturated <- 611*exp((Lv/461)*((1/273.15)-(1/kelvin)))
  return(saturated)
} #satvap<-function(T.c){

# function to calculate aerodynamic resistance
raero.f <- function(Ur=1,zr=50,z0=z0,d=0,rho=1){
  #arguments:  Ur is reference windspeed [m/s] at reference height zr
  #            zr is reference height [m] where Ur applies
  #            z0 is roughness length [m]; 0.01m is typical value for crop
  #            d is displacement height [m]
  #            rho is air density [kg/m^3]
  k <- 0.4  # von Karman constant
  
  CD <- (k^2)/(log((zr-d)/z0))^2  # aerodynamic transfer coefficient
  raero <- 1/(CD*Ur)              # aerodynamic resistance [s/m]
  return(raero)
} #raero.f<-function(){
#################################################

#################################################
# Land surface characteristics
gvmax <- 1/50      # max vegetation conductance [m/s] (reciprocal of vegetation resistance) when vegcontrol is FALSE;  when TRUE, calculated by BBF function
albedo.surf <- 0.1    # surface albedo
albedo <- albedo.surf # surface albedo
z0 <- 0.5          # roughness length [m]
epsilon.s <- 0.97  # surface emissivity for forest, according to Jin & Liang (2006)
LAI <- 3.0         # average leaf area index; for a forest like Harvard Forest, ~3.0 over the year [.]
Kb <- 0.5          # extinction coefficient within plant canopy [.]; average value ~0.5:  https://link.springer.com/article/10.1007/s11707-014-0446-7
Q10 <- 2           # temperature-dependence of respiration:  what factor does rate increase for 10-deg C increase
Resp25 <- 2        # respiration rate at 25-deg C [umole CO2/m2/s]
# heat capacity of land surface
# a) heat capacity based on soil
# D <- 0.1*(1/24)          # the depth of soil that temp fluctuations would penetrate [m]; 0.1m is roughly the depth that would penetrate on diurnal timescales
# Cp.soil <- 1921          # specific heat of soil organic material [J/kg/K]
# rho.soil <- 1300         # density of soil organic material [kg/m3]
# Cs <- Cp.soil*rho.soil*D # heat capacity of organic soil [J/K/m2]
# b) heat capacity based on vegetation
Hveg <- 10               # height of vegetation [m]
rho.veg <- 0.7E6         # wood density [g/m3]
Cp.veg <- 3000           # bulk heat capacity of above-ground vegetation [J/kg/K];  Sect. 7.2 of Bonan (2019)
Cs <- Cp.veg*(rho.veg/1000)*Hveg/100  # heat capacity of vegetation [J/K/m2]

# soil parameters from Clapp & Hornberger (1978); taken from CLASS model (https://github.com/classmodel/modelgui/blob/master/landsoil.cpp)
# select soil type from one of below
soiltype <- "Sandy loam"
# soiltype <- "Sand"
# soiltype <- "Clay"
if (tolower(soiltype)=="sandy loam") {
  # a) Sandy loam soil
  Wsat <- 0.472  # saturated volumetric water content [m3/m3]
  Wfc <- 0.323   # volumetric water content at field capacity [m3/m3]
  Wwilt <- 0.171 # volumetric water content at wilting point [m3/m3]
  aa <- 0.219    # Clapp & Hornberger (1978) retention parameter a
  bb <- 4.9      # Clapp & Hornberger (1978) retention parameter b
  pp <- 4        # Clapp & Hornberger (1978) retention parameter c
  rTsoil.sat <- 3.56E-6 # saturated soil thermal insulance factor [K m2 J-1]
  C1sat <- 0.132
  C2ref <- 1.8 
} else if (tolower(soiltype)=="sand") {
  # b) Sand
  Wsat <- 0.403  # saturated volumetric water content [m3/m3]
  Wfc <- 0.244   # volumetric water content at field capacity [m3/m3]
  Wwilt <- 0.059 # volumetric water content at wilting point [m3/m3]
  aa <- 0.387    # Clapp & Hornberger (1978) retention parameter a
  bb <- 4.05     # Clapp & Hornberger (1978) retention parameter b
  pp <- 4        # Clapp & Hornberger (1978) retention parameter c
  rTsoil.sat <- 3.222E-6 # saturated soil thermal insulance factor [K m2 J-1]
  C1sat <- 0.082
  C2ref <- 3.9 
} else if (tolower(soiltype)=="clay") {
  # c) Clay
  Wsat <- 0.614  # saturated volumetric water content [m3/m3]
  Wfc <- 0.541   # volumetric water content at field capacity [m3/m3]
  Wwilt <- 0.335 # volumetric water content at wilting point [m3/m3]
  aa <- 0.083    # Clapp & Hornberger (1978) retention parameter a
  bb <- 11.4     # Clapp & Hornberger (1978) retention parameter b
  pp <- 12        # Clapp & Hornberger (1978) retention parameter c
  rTsoil.sat <- 3.6E-6 # saturated soil thermal insulance factor [K m2 J-1]
  C1sat <- 0.342
  C2ref <- 0.3 
} else { 
  stop (paste("Need to select valid soil type:",soiltype))  
}

Lambda <- 5.9       # thermal transfer coefficient of surface layer [W/m2/K]
# Initialize two-layer (force-restore) soil model, from de Arellano et al. (2015)
Tsoil2 <- 286       # T of deep soil layer [K] that is constant
Tsoil1 <- Tsoil2    # T of top soil layer [K] that varies w/ time
Wsoil2 <- Wfc       # volumetric water content of deep soil layer [m3/m3]
Wsoil1 <- Wsoil2    # volumetric water content of top soil layer [m3/m3]
d1 <- 0.1           # soil depth to which diurnal variations in moisture penetrates [m], from Deardorff [1977]
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
# a) hourly varying SWdn
SWdn<--15*(t.hr-12)^2+700 # hourly downward shortwave radiation [W/m2]
names(SWdn)<-t.hr
SWdn[as.character(c(0:5,19:24))]<-0 # night time:  set to 0
SWdn_DAY <- SWdn

# b) constant SWdn
# SWdn[1:length(SWdn)]<-1000

# Downward longwave radiation (over-written with dynamically varying LWdn when LWdnTF set to TRUE)
LWdn <- SWdn; LWdn[1:length(LWdn)] <- 300 # constant downward longwave radiation [W/m2]
LWdn_DAY <- LWdn

# -----------Atmospheric conditions----------#
# Air temperature
Ta.c<- -0.5*(t.hr-12)^2+30  # PRESCRIBED air temperature [deg-C]
names(Ta.c) <- t.hr
Ta.c[1:length(Ta.c)] <- 5    # override with CONSTANT air temperature [deg-C]
Ta.c_DAY <- Ta.c

# specific humidity of air:  determine from RH, air temperature
RH <- 0.9
e <- RH*satvap(mean(Ta.c))/100  #vapor pressure [hPa]
Psurf <- 1000     #surface pressure [hPa] 
qa.presc <- (Rd/Rv)*e/Psurf   #prescribed specific humidity [g/g]
Hscale <- 8000    # scale height of atmosphere--i.e., height at which Psurf decays to (1/e) [m]
hmin <- 200       # minimum height of atmospheric boundary layer [m]
thetavM0<-(Ta.c[1]+273.15)*(1+0.61*qa.presc) # initial virtual potential temperature [K]; Eq. 1.5.1b of Stull [1988]
Beta <- 0.2       # closure hypothesis:  fraction of surface virtual potential temperature flux that determines entrainment heat flux
gamma <- 5/1000   # slope of thetav above growing ABL [K/m]
qabove <- qa.presc/5  # specific humidity of air above ABL [g/g]
W <- 0            # subsidence rate [m/s]
Ur <- 1           # reference windspeed [m/s] at reference height zr
zr <- 50          # reference height [m] where Ur applies
Cair <- 400       # atmospheric CO2 concentration [umole/mole, or ppm]; this is also the initial CO2 value within ABL if co2budgetTF = TRUE
Cfree <- 400      # CO2 concentration [ppm] in free troposphere (not modified by values in ABL)
Cabove <- Cfree   # CO2 concentration [ppm] above ABL (later modified by value in residual layer)
albedo.cloud <- 0.5 # albedo of cloud

# parameters determining CO2 greenhouse effect
CO2.SENSITIVITY <- 3.7  # CO2 doubling sensitivity [W/m2 per doubling of CO2]  
CO2.baseline <- 280     # baseline to determine doubling (pre-industrial CO2 concentration [ppm])
#################################################


# ave CO2 in atmospheric column, using scale height as weighting (i.e., density follows exponential decay)
CO2.colave <- Cair + (Cfree - Cair)*exp(-hmin/Hscale)           
GHG.FORCE <- CO2.SENSITIVITY*log(CO2.colave/CO2.baseline)/log(2) # GHG forcing--from CO2 elevated above CO2base.ppm [W/m2]

# initialize T with equilibrium value (determined through "uniroot")
f <- function(T, Ta, SWdn, LWdn, albedo.cloud, albedo.surf, epsilon.s, Tsoil1, Ur,
              zr, z0, gvmax=gvmax, RH=RH, CO2=Cair, Psurf=1000, Hscale=8000, GHG.FORCE=GHG.FORCE){  
  # --------------Physical constants--------#
  Cp <- 1005.7; Cv <- 719 # heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel [1994])
  g <- 9.80665 # standard surface gravity [m/s2]
  Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel [1994])
  Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson [1999])
  sigma <- 5.670373E-8 # Stefan-Boltzmann constant [W/m2/K4]
  # --------------Physical constants--------#
  
  if(cloudTF){
    # diagnose cloud fraction based on Eq. 3 of Slingo [1987]:  "The development and verification of a cloud prediction scheme for the ECMWF model"
    RHcrit <- 0.8
    tmp <- (RH-RHcrit)/(1-RHcrit)
    tmp[tmp<0] <- 0
    cloud <- tmp^2
    if(cloud > 1.0)cloud <- 1.0
  } else { cloud <- 0 } # if(cloudTF){
  albedo <- (1-cloud)*albedo.surf + cloud*albedo.cloud
  SWup <- albedo*SWdn
  Ta.c <- Ta - 273.15
  e <- RH*satvap(mean(Ta.c))/100  #vapor pressure [hPa]
  LWup <- epsilon.s*sigma*T^4
  if (LWdnTF) {
    # empirical formula of downward longwave radiation based on Yang et al. (2023): https://doi.org/10.5194/acp-23-4419-2023
    epsilon.clr <- 0.532 + 0.808*((e/Ta)^(1/3))  # clear-sky emissivity
    epsilon.all <- epsilon.clr*(1-0.201*cloud^0.796) + 0.088*(cloud^1.038)*((RH*100)^0.221)  # all-sky emissivity
    LWdn <- epsilon.all * sigma * Ta^4 
    LWdn <- LWdn + GHG.FORCE  # add GHG forcing
  } # if (LWdnTF) {
  Rn <- SWdn - SWup + LWdn - LWup
  
  # determine sensible heat flux
  rho.surf <- Psurf*100/(Rd*T)   # surface air density [kg/m3]
  raero <- raero.f(Ur=Ur,zr=zr,z0=z0,rho=rho.surf)
  H <- (Cp*rho.surf/(raero))*(T-Ta)   # [W/m2]
  
  # determine latent heat flux
  Lv <- 1000*latentheat(T-273.15)  # latent heat of vaporization [J/kg]
  esat <- satvap(T-273.15)/100 # saturation specific humidity [hPa]
  VPD <- 100*(esat-e)          # vapor pressure deficit [Pa]
  qsat <- (Rd/Rv)*esat/Psurf  # saturation specific humidity [g/g]
  if (vegcontrolTF) {
    # Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model for vegetation resistance [s/m]
    hs <- e/esat  # fractional humidity (=1/RH) at leaf surface [.]
    cs <- CO2     # CO2 concentration at leaf surface [umole/mole]
    BBFout <- BBF(SW=SWdn,Tleaf.C=T-273.15,hs=hs,beta.W=1.0,cs=cs,Psurf=Psurf)  
    gsv <- BBFout["gsv"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]  
    rho.mole <- rho.surf*1000/Md # air density [kg/m3] => molar density [moles/m3]
    gsv <- gsv/rho.mole   # [mole/m2/s] => [m/s]
    rveg <- 1/gsv           # vegetation resistance [s/m]
    An <- BBFout["An"]    # Net photosynthesis [umole/m2/s]
    ci <- BBFout["ci"]    # intercellular CO2 [umole/mole]
  } else {
    rveg <- 1/gvmax
    An <- NA; ci <- NA
  } # if(vegcontrolTF){
  # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
  scale.canopy<-(1-exp(-Kb*LAI))/Kb
  An <- An*scale.canopy
  gveg <- (1/rveg)*scale.canopy
  rveg <- 1/gveg
  LE <- (Lv*rho.surf/(raero + rveg))*(qsat - qa.presc) #[W/m2]
  if(LE<0) LE <- 0
  
  # determine ground heat flux from two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
  G <- Lambda * (T - Tsoil1)
  
  # this should =0 when T is at equilibrium value
  return(Rn - H - LE - G)
} # f<-function(T,Ta,SWdn,LWdn,albedo,epsilon.s){

xinterv <- Ta.c[1]+273.15+c(-50,50)  # interval over which to search for equil temperature
# use initial radiation, temps to solve for initial equil. temperature
Tinit <- uniroot(f,interval=xinterv,Ta=Ta.c[1]+273.15,SWdn=SWdn[1],LWdn=LWdn[1],Tsoil1=Tsoil1,albedo.cloud=albedo.cloud,
                 albedo.surf=albedo.surf,epsilon.s=epsilon.s,Ur=Ur,zr=zr,z0=z0,gvmax=gvmax,RH=RH,Psurf=Psurf,Hscale=Hscale,GHG.FORCE=GHG.FORCE)$root
tmp <- f(T=Tinit,Ta=Ta.c[1]+273.15,SWdn=SWdn[1],LWdn=LWdn[1],Tsoil1=Tsoil1,albedo.cloud=albedo.cloud,
         albedo.surf=albedo.surf,epsilon.s=epsilon.s,Ur=Ur,zr=zr,z0=z0,gvmax=gvmax,RH=RH,Psurf=Psurf,Hscale=Hscale,GHG.FORCE=GHG.FORCE)
print(paste("Tinit [oC]:",signif(Tinit-273.15,5),";   (Rn-H-LE-G) =",signif(tmp,4),"[W/m2]"))

#############################################################################################################
#--------------------------------- ODE ODE ODE ODE ODE ODE ODE ODE ODE ODE ---------------------------------#
#############################################################################################################
# initialize state variables
thetaM <- Ta.c[1]+273.15
qa <- qa.presc   # initialize with prescribed specific humidity [g/g]
thetavM <- thetaM*(1+0.61*qa)   # virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]

yini <- c(T=Tinit, Ta=Ta.c[1]+273.15, qa=qa, thetavM=thetavM,
          Tsoil1=Tsoil1, Wsoil1=Wsoil1, h=hmin, CO2=Cair) 
names(yini) <- c("T","Ta","qa","thetavM","Tsoil1","Wsoil1","h","CO2")

########################################################
# initialize parameters
# 0. numerical parameters
parms <- c(dt=dt,DTtol=DTtol,countTmax=countTmax)
# 1.  flags
parms <- c(parms,vegcontrolTF=vegcontrolTF,atmrespondTF=atmrespondTF,ABLTF=ABLTF,
           soilWTF=soilWTF,co2budgetTF=co2budgetTF)
# 2.  atmospheric conditions 
parms <- c(parms,Psurf=Psurf,qa.presc=qa.presc,Hscale=Hscale,hmin=hmin,Beta=Beta,
           gamma=gamma,qabove=qabove,W=W,Ur=Ur,zr=zr,Cabove=Cabove,Cfree=Cfree,albedo.cloud=albedo.cloud)
# 3.  land surface characteristics
parms <- c(parms,gvmax=gvmax,albedo.surf=albedo.surf,z0=z0,epsilon.s=epsilon.s,
           LAI=LAI,Kb=Kb,Hveg=Hveg,rho.veg=rho.veg,Cp.veg=Cp.veg,Cs=Cs,Resp25=Resp25,Q10=Q10)
# 4.  soil characteristics
parms <- c(parms,Wsat=Wsat,Wfc=Wfc,Wwilt=Wwilt,aa=aa,bb=bb,pp=pp,rTsoil.sat=rTsoil.sat,
           C1sat=C1sat,C2ref=C2ref,Lambda=Lambda,Tsoil2=Tsoil2,Wsoil2=Wsoil2,d1=d1)

########################################################
# define LAIM model function (what happens each time step)
LAIM <-function(time,state,parms,SWdn_DAY,LWdn_DAY,Ta.c_DAY){
  #------------------#
  # Physical constants
  Cp <- 1005.7;Cv <- 719 # heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel (1994)
  g <- 9.80665 # standard surface gravity [m/s2]
  Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel (1994))
  Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson (1999)
  sigma <- 5.670373E-8    # Stefan-Boltzmann constant [W/m2/K4]
  Md <- 28.97  #molar mass of dry air [g/mole]
  rho.W <- 1000 # density of water [kg/m3]
  #------------------#
  
  if(((time/3600)%%1)==0) print(paste("Running model: time=",time/3600,"[hr]"))
  with(as.list(c(state,parms)),{
    
    if(cloudTF){
      # diagnose cloud fraction based on Eq. 3 of Slingo [1987]:  "The development and verification of a cloud prediction scheme for the ECMWF model"
      RH <- e/(satvap(Ta - 273.15)/100)
      RHcrit <- 0.8
      tmp <- (RH-RHcrit)/(1-RHcrit)
      tmp[tmp<0] <- 0
      cloud <- tmp^2
      if(cloud > 1.0)cloud <- 1.0
    } else { cloud <- 0 } # if(cloudTF){
    albedo <- (1-cloud)*albedo.surf + cloud*albedo.cloud
    SWdn.t <- approx(x=as.numeric(names(SWdn_DAY))*3600,y=SWdn_DAY,xout=time%%(24*3600))$y  # downward shortwave radiation [W/m2]
    SWup <- albedo*SWdn.t
  
    # ave CO2 in atmospheric column, using scale height as weighting (i.e., density follows exponential decay)
    #     NOTE: ignore the variation in CO2 within shallow residual layer (represented by updated Cabove) 
    CO2.colave <- CO2 + (Cfree - CO2)*exp(-h/Hscale)                 
    GHG.FORCE <- CO2.SENSITIVITY*log(CO2.colave/CO2.baseline)/log(2) # GHG forcing--from CO2 elevated above CO2.baseline

    LWup <- epsilon.s*sigma*T^4   # upward longwave radiation [W/m2]
    LWdn.t <- approx(x=as.numeric(names(LWdn_DAY))*3600,y=LWdn_DAY,xout=time%%(24*3600))$y  # downward shortwave radiation [W/m2]
    if (LWdnTF) {
      # empirical formula of downward longwave radiation based on Yang et al. (2023): https://doi.org/10.5194/acp-23-4419-2023
      epsilon.clr <- 0.532 + 0.808*((e/Ta)^(1/3))  # clear-sky emissivity
      epsilon.all <- epsilon.clr*(1-0.201*cloud^0.796) + 0.088*(cloud^1.038)*((RH*100)^0.221)  # all-sky emissivity
      LWdn.t <- epsilon.all * sigma * Ta^4 
      LWdn.t <- LWdn.t + GHG.FORCE  # add GHG forcing
    } # if (LWdnTF) {
    
    # determine net radiation
    Rn <- SWdn.t-SWup+LWdn.t-LWup
    
  countT <- 0; iterateT <- TRUE
  while (iterateT) {   #iterate until convergence
    countT <- countT + 1
    if(countT > countTmax)stop("T does not converge")
      
    if (!atmrespondTF) {
      Ta <- approx(x=as.numeric(names(Ta.c_DAY))*3600,y=Ta.c_DAY,xout=time%%(24*3600))$y+273.15  #use prescribed value
      qa <- qa.presc
    } # if(atmrespondTF){
    
    # determine sensible heat flux
    rho.surf <- Psurf*100/(Rd*T)   # surface air density [kg/m3]
    raero <- raero.f(Ur=Ur,zr=zr,z0=z0,rho=rho.surf)
    H <- (Cp*rho.surf/(raero))*(T-Ta)   # [W/m2]
      
    # determine latent heat flux
    beta.W <- 1   # water stress parameter (dependent on soil moisture)
    Lv <- 1000*latentheat(T-273.15)  # latent heat of vaporization [J/kg]
    esat <- satvap(T-273.15)/100     # saturation specific humidity [hPa]
    e <- qa*Psurf/(Rd/Rv)            # vapor pressure [hPa]
    VPD <- 100*(esat-e)              # vapor pressure deficit [Pa]
    qsat <- (Rd/Rv)*esat/Psurf       # saturation specific humidity [g/g]
    if (vegcontrolTF) {
      if (soilWTF) {
        # Eq. (12.56) of Bonan (2019)
        beta.W <- (Wsoil1 - Wwilt)/(Wfc - Wwilt)
        if (Wsoil1 >= Wfc) beta.W <- 1.0
        if (Wsoil1 <= Wwilt) beta.W <- 0
      } # if (soilWTF)
      # Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model for vegetation resistance [s/m]
      hs <- e/esat  # fractional humidity (=1/RH) at leaf surface [.]   
      if(hs<0.7) hs <- hs + 0.3   #!!! quick adjustment that ensures leaf surface is not too dry...accounts for higher humidity within canopy  !!!#
      cs <- CO2    # CO2 concentration at leaf surface [umole/mole]
      BBFout <- BBF(SW=SWdn.t,Tleaf.C=T-273.15,hs=hs,beta.W=beta.W,cs=cs,Psurf=Psurf)  
      gsv <- BBFout["gsv"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]  
      rho.mole <- rho.surf*1000/Md # air density [kg/m3] => molar density [moles/m3]
      gsv <- gsv/rho.mole   # [mole/m2/s] => [m/s]
      rveg <- 1/gsv         # vegetation resistance [s/m]
      An <- BBFout["An"]    # Net photosynthesis [umole/m2/s]
      ci <- BBFout["ci"]    # intercellular CO2 [umole/mole]
    } else {
      rveg <- 1/gvmax
      An <- NA; ci <- NA
    } # if(vegcontrolTF){
      
    # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
    scale.canopy<-(1-exp(-Kb*LAI))/Kb
    An <- An*scale.canopy
    gv <- (1/rveg)*scale.canopy
    rveg <- 1/gv
    LE <- (Lv*rho.surf/(raero+rveg))*(qsat-qa) #[W/m2]
    if(LE<0) LE <- 0
    
    # determine respirational flux of CO2 to atmosphere
    Resp <- Resp25*(Q10^((T-298.15)/10))  # respiration based on Q10 formulation
          
    # determine ground heat flux 
    # use two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
    G <- Lambda * (T - Tsoil1)
      
    Storage <- Rn - LE - H - G
    # update temperature 
    DT <- (Storage/Cs)*dt
    T <- T+DT
    #print(paste("iterating so that T converges:",countT,paste("T =",signif(T,5)),signif(DT,4)))
    iterateT <- abs(DT)>DTtol   # continue iterating until T converges
  } # while (iterateT) {   #iterate until converge

  
  # heat transport between surface and deep soil layer to update Tsoil1 from CLASS model
  rTsoil <- rTsoil.sat * (Wsat/Wsoil2)^(bb/(2*log(10)))
  dTsoil1.dt <- (rTsoil*G - (2*pi/(86400))*(Tsoil1 - Tsoil2)) #Eq. (9.32) of de Arellano et al. (2015)
  
  if (soilWTF) {
    # update soil water content, based on CLASS model
    C1 <- C1sat*(Wsat/Wsoil1)^(bb/2 + 1)  #Eq. (9.35) of de Arellano et al. (2015)
    Wsmall <- 1E-3
    C2 <- C2ref*(Wsoil2/(Wsat - Wsoil2 + Wsmall))  #Eq. (9.36) of de Arellano et al. (2015)
    Wsoil1eq <- Wsoil2 - aa*Wsat*((Wsoil2/Wsat)^pp)*(1-(Wsoil2/Wsat)^(8*pp))  #Eq. (9.37) of de Arellano et al. (2015)
    # Eq. (9.34) of de Arellano et al. (2015); NOTE:  use LE instead of LEsoil as in (9.34), and -1 multiplied by C1 that is missing in (9.34)
    dWsoil1.dt <- ((-C1/(rho.W*d1))*(LE/Lv) - (C2/86400)*(Wsoil1 - Wsoil1eq))
    # make sure that Wsoil1 does not dip below Wwilt;  NOTE:  this does NOT conserve water (since could stll have residual E from minimum gv)
    if(Wsoil1 < Wwilt){dWsoil1.dt <- (Wwilt-Wsoil1)/dt;Wsoil1 <- Wwilt}  
    if (Wsoil1 < 0) {dWsoil1.dt <- (0-Wsoil1)/dt;Wsoil1 <- 0}
  } else {
    dWsoil1.dt <- 0
  } #if (soilWTF) {
  
  # if want atmosphere to respond
  # Based on "zero-order jump" or "slab" model of convective boundary layer, described in Pg. 151~155 of Garratt [1992]
  CO2flux.veg <- NA; CO2flux.ent <- NA; CO2flux.tot <- NA
  if (atmrespondTF) {
    #calculate surface virtual heat flux
    Lv <- latentheat(T-273.15)  # latent heat of vaporization [J/g]
    E <- LE/Lv   # surface moisture flux [g/m^2/s] 
    F0theta <- H/Cp  # potential heat flux [K-kg/m^2/s]
    F0thetav <- F0theta+0.073*Lv*E/Cp # virtual heat flux [K-kg/m^2/s]
    if (ABLTF) {
      Fhthetav <- -1*Beta*F0thetav   # closure hypothesis (Eq. 6.15 of Garratt [1992])
      # calculate ABL growth rate
      dh.dt<-(1+2*Beta)*F0thetav/(gamma*h)
      if (F0thetav<=0.00){dh.dt <- (hmin - h)/dt} # override value:  ABL collapses
    } else {
      dh.dt <- 0
      Fhthetav <- 0
    } # if(ABLTF){
    
    rhobar <- rho.surf*(1-exp(-h/Hscale))*(Hscale/h)  # determine ABL-averaged air density [kg/m3]
    
    # calculate entrainment flux of humidity
    deltaq <- (qabove - qa)
    Fhq <- 0
    if(dh.dt>=0)Fhq <- -1*rhobar*deltaq*(dh.dt-W)*1000  # entrainment flux of humidity [g/m2/s] NOTE:  assume CONSTANT air density!
    dq.dt <- (E - Fhq)/(rhobar*1000*h) # change of humidity in ABL [1/s]
    
    # update ABL-averaged thetav
    dthetavM.dt <- (F0thetav - Fhthetav)/h   # change of thetav in ABL [K-kg/m^3/s]
    dthetavM.dt <- dthetavM.dt/rhobar        # [K-kg/m^3/s]=>[K/s]
    
    # update ABL-averaged CO2
    dC.dt <- 0
    if (co2budgetTF) {
      CO2flux.veg <- (-1*An + Resp)  # surface CO2 flux [umole/m2/s]; photosynthesis is a negative flux (removal from atmosphere)
      CO2flux.tot <- CO2flux.veg
      CO2flux.ent <- 0
      if(dh.dt>0){
        CO2flux.ent<-(rhobar/(Md/1000))*(dh.dt - W)*(Cabove - CO2)   # entrainment flux of CO2 [umole/m2/s]
        CO2flux.tot <- CO2flux.veg + CO2flux.ent
      } # if(dh.dt>0){
      dC.dt <- CO2flux.tot*(Md/1000)/(rhobar*h)  # dilute surface flux in box of height h to generate change in CO2 [ppm/s]
    } # if (co2budgetTF) {
    
  } else{
    Lv <- latentheat(T-273.15)  # latent heat of vaporization [J/g]
    E <- LE/Lv   # surface moisture flux [g/m^2/s] 
    Fhq <- 0
    deltaq <- 0
    dC.dt <- 0
    dthetavM.dt <- 0
    dq.dt <- 0
    dh.dt <- 0
  } # if(atmrespondTF){
  
  
  # derivatives of variables--need to be returned as part of call to 'ode'
  DT <- DT
  DTa <- dthetavM.dt/(1+0.61*qa)
  Dqa <- dq.dt
  DthetavM <- dthetavM.dt
  DTsoil1 <- dTsoil1.dt
  DWsoil1 <- dWsoil1.dt
  Dh <- dh.dt
  DCO2 <- dC.dt 
    
  #variables that aren't integrated with time and aren't returned as derivatives
  vars2<-c(SWdn=SWdn.t,LWdn=LWdn.t,GHG.FORCE=GHG.FORCE,Rn=Rn,LWup=as.numeric(LWup),H=as.numeric(H),LE=as.numeric(LE),G=G,RH=RH,cloud=cloud,albedo=albedo,
           qsat=as.numeric(qsat),An=as.numeric(An),rveg=as.numeric(rveg),raero=raero,beta.W=as.numeric(beta.W),
           CO2flux.veg=as.numeric(CO2flux.veg),CO2flux.ent=as.numeric(CO2flux.ent),CO2flux.tot=as.numeric(CO2flux.tot),
           dh.dt=as.numeric(dh.dt),E=as.numeric(E),Fhq=as.numeric(Fhq),deltaq=as.numeric(deltaq))
 
  return(list(c(DT,DTa,Dqa,DthetavM,DTsoil1,DWsoil1,Dh,DCO2),vars2))
  })
} # LAIM <-function(time,state,parms,SWdn_TIME,Ta_TIME){


########################################################
# Time integration: call LAIM function using ode()
if(atmrespondTF&ABLTF&t.day>1){
  print(paste("==========Multiple calls to ode:=========="))
  result <- NULL
  for(t.dd in 1:(ceiling(max((times)/(3600*24))))){
    print(paste("TIME INTEGRATION: DAY",t.dd))
    sel <- times>=(3600*24*(t.dd-1))&times<(3600*24*(t.dd))  #!!! still potentially missing last time step on last day !!!
    times.sub <- times[sel]
    # multiple calls to "ode", each time by 1 day, to allow for entrainment of residual layer
    if(t.dd>1){
      ilast <- nrow(result.tmp)
      #  find qa in ABL just before the ABL collapses, and use it as the humidity in residual layer that would be entrained into ABL following day
      qa.resid <- result.tmp[which(result.tmp$h==max(result.tmp$h)),"qa"]
      print(paste("specific humidity of residual layer [g/g]:",signif(qa.resid,5)))
      parms["qabove"] <- qa.resid   # assign residual layer humidity as humdity above ABL
      #  find [CO2] in ABL just before the ABL collapses, and use it as the [CO2] in residual layer that would be entrained into ABL following day
      Cair.resid <- result.tmp[which(result.tmp$h==max(result.tmp$h)),"CO2"]
      print(paste("CO2 of residual layer [ppm]:",signif(Cair.resid,5)))
      parms["Cabove"] <- Cair.resid   # assign residual layer [CO2] as [CO2] above ABL
      yini <- c(T=result.tmp$T[ilast], Ta=result.tmp$Ta[ilast], qa=result.tmp$qa[ilast], thetavM=result.tmp$thetavM[ilast],
                Tsoil1=result.tmp$Tsoil1[ilast], Wsoil1=result.tmp$Wsoil1[ilast], h=result.tmp$h[ilast], CO2=result.tmp$CO2[ilast])
      names(yini) <- c("T","Ta","qa","thetavM","Tsoil1","Wsoil1","h","CO2")
    } # if(t.dd>1){
    result.tmp <- ode(yini, times.sub, LAIM, parms, SWdn_DAY=SWdn_DAY, LWdn_DAY=LWdn_DAY,Ta.c_DAY=Ta.c_DAY, method = "lsoda")
    result.tmp <- data.frame(result.tmp)
    result <- rbind(result,result.tmp)
  } # for(t.dd in 1:(t.day+1)){
  filenm <- "result.csv"; write.csv(result,file=filenm)
  print(paste(filenm,"written out"))
} else {
  # single call "ode" to integrate LAIM model in time
  result <- ode(yini, times, LAIM, parms, SWdn_DAY=SWdn_DAY, LWdn_DAY=LWdn_DAY,Ta.c_DAY=Ta.c_DAY, method = "lsoda")
  result <- data.frame(result)
  filenm <- "result.csv"; write.csv(result,file=filenm)
  print(paste(filenm,"written out"))
} # if(atmrespondTF&ABLTF&t.day>1){


########################################################
# Plotting 
# text on plot 
xmain <- paste("atmrespondTF=",atmrespondTF)
xmain <- paste(xmain,"  ABLTF=",ABLTF)
xmain <- paste(xmain,"  cloudTF=",cloudTF)
xmain <- paste(xmain,"\nvegcontrolTF=",vegcontrolTF)
xmain <- paste(xmain,"  soilWTF=",soilWTF)
xmain <- paste(xmain,"\ndt=",dt,"[s]")
# regenerate VPD from qsat and qa 
e <- result[,"qa"]*Psurf/(Rd/Rv)      # vapor pressure [hPa]
esat <- result[,"qsat"]*Psurf*(Rv/Rd)  # saturation vapor pressure [hPa]
VPD <- 100*(esat-e)                     # vapor pressure deficit [Pa]

# generate 4 different plots in one window, with 2*2 configuration
dev.new(); par(mfrow=c(2,2),cex.main=0.7)   
ylims <- range(result[,c("T","Ta","Tsoil1")]-273.15)
plot(result[,"time"]/3600,result[,"T"]-273.15,type="l",xlab="Time [hour]",ylab="Temperature [deg-C]",
     cex.axis=1.3,cex.lab=1.3,lwd=3,ylim=ylims,main=xmain)
lines(result[,"time"]/3600,result[,"Ta"]-273.15,lty=3,lwd=1.5)
lines(result[,"time"]/3600,result[,"Tsoil1"]-273.15,lty=1,lwd=2,col="darkgray")
legend(x="topright",c("Ts","Ta","Tsoil"),lwd=c(3,2,2),lty=c(1,3,1),col=c("black","black","darkgray"))

plot(result[,"time"]/3600,VPD,type="l",xlab="Time [hour]",ylab="VPD [Pa]",lwd=2,
     cex.axis=1.3,cex.lab=1.3,main=xmain)

ylims <- range(result[,c("qsat","qa")]*1000)
plot(result[,"time"]/3600,result[,"qsat"]*1000,type="l",xlab="Time [hour]",ylab="Specific humidity [g/kg]",
     cex.axis=1.3,cex.lab=1.3,lwd=3,ylim=ylims,main=xmain)
lines(result[,"time"]/3600,result[,"qa"]*1000,lty=3,lwd=1.5)
legend(x="topright",c("qsat","qa"),lwd=c(3,2),lty=c(1,3))

ylims <- range(result[,c("raero","rveg")])
plot(result[,"time"]/3600,result[,"raero"],type="l",xlab="Time [hour]",ylab="Resistances [s/m]",
     cex.axis=1.3,cex.lab=1.3,lwd=1.5,lty=3,ylim=ylims,main=xmain,col="black")
lines(result[,"time"]/3600,result[,"rveg"],lty=1,lwd=3,col="black")
legend(x="topright",c("r_veg","r_aero"),lwd=c(3,2),lty=c(1,3),col=c("black","black"))
dev.copy(png,"T_q_r.png");dev.off();print("T_q_r.png written out")

# plot with energy fluxes 
dev.new()
matplot(result[,"time"]/3600,result[,c("Rn","LWdn","LWup","H","LE","G")],type="l",lty=c(1,3,1,1,1,1),lwd=c(3,2,3,2,2,2),
        cex.axis=1.5,cex.lab=1.5,col=c("black","black","darkgray","orange","blue","darkgreen"),xlab="Time [hr]",ylab="")
mtext(text=expression(paste("Energy Fluxes [W ",m^-2,"]",sep="")),line=2.3,cex=1.4,side=2)
legend(x="topright",c("Rn","LWdn","LWup","H","LE","G"),col=c("black","black","darkgray","orange","blue","darkgreen"),lwd=c(3,2,3,2,2,2),lty=c(1,3,1,1,1,1))
title(main=xmain)
dev.copy(png,"Energyfluxes.png");dev.off();print("Energyflux.png written out")

# plot soil water content 
if (soilWTF) {
  dev.new()
  plot(result[,"time"]/3600,result[,"Wsoil1"],type="l",xlab="Time [hour]",ylab="",cex.main=1.0,
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=paste("Soil type =",soiltype,"\n",xmain),col="black",ylim=c(Wwilt,Wfc))
  mtext(text=expression(paste("Soil Volumetric Water Content [",m^3,"/",m^3,"]",sep="")),line=2,cex=1.3,side=2)
  abline(h=Wsoil2,lty=3,lwd=2)
  legend(x="bottomright",c("Wsoil1","Wsoil2"),col=c("black","black"),lwd=2,lty=c(1,3))
  dev.copy(png,"Wsoil.png");dev.off();print("Wsoil.png written out")
} #if (soilWTF)

if (atmrespondTF) {
  # plot time series of ABL height 
  dev.new()
  plot(result[,"time"]/3600,result[,"h"],type="l",xlab="Time [hour]",ylab="ABL height  h(t) [m]",
       cex.axis=1.3,cex.lab=1.3,lwd=2)
  title(main=paste0(xmain,";  Beta=",Beta,";  gamma=",signif(gamma,4)," [K/m]"),cex.main=1.2)
  dev.copy(png,"ABLht.png");dev.off();print("ABLht.png written out") 
} #if(atmrespondTF){

if (vegcontrolTF) {
  # plot time series of photosynthetic uptake 
  dev.new()
  plot(result[,"time"]/3600,result[,"An"],type="l",xlab="Time [hour]",ylab="",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain)
  mtext(text=expression(paste("Net Photosynthesis (An) [",mu,"mole ",m^-2," ",s^-1,"]",sep="")),line=2,cex=1.3,side=2)
  dev.copy(png,"PSN.png");dev.off();print("PSN.png written out")
} #if(vegcontrolTF){

if (vegcontrolTF&atmrespondTF&co2budgetTF) {
  # plot time series of CO2 
  dev.new()
  plot(result[,"time"]/3600,result[,"CO2"],type="l",xlab="Time [hour]",ylab="CO2 [ppm]",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain)
  par(new=TRUE)
  ylims <- range(result[,c("CO2flux.ent","CO2flux.veg")],na.rm=TRUE)
  plot(result[,"time"]/3600,result[,"CO2flux.veg"],type="l",axes=F,xlab="",ylab="",col="darkgray",ylim=ylims,lty=1,lwd=2)
  lines(result[,"time"]/3600,result[,"CO2flux.ent"],col="darkgray",lty=3,lwd=2)
  abline(h=0,lty=1,lwd=0.5,col="darkgray")
  axis(4,cex.lab=1.3,cex.axis=1.3,col="darkgray",col.axis="darkgray")
  legend(x="topright",c("CO2tot","dCO2.veg","dCO2.ent"),lwd=2,lty=c(1,1,3),
         col=c("black","darkgray","darkgray"),text.col=c("black","darkgray","darkgray"))
  dev.copy(png,"CO2.png");dev.off();print("CO2.png written out")
} #if (vegcontrolTF&atmrespondTF) {

if (cloudTF) {
  # plot time series of cloud fraction, albedo, and relative humidity
  dev.new()
  ylims <- range(result[,c("cloud","albedo","RH")],na.rm=TRUE)
  plot(result[,"time"]/3600,result[,"cloud"],type="l",xlab="Time [hour]",ylab="Cloud Fraction/Albedo/RH",
       cex.axis=1.3,cex.lab=1.3,lwd=3,lty=1,main=xmain,ylim=ylims)
  lines(result[,"time"]/3600,result[,"albedo"],type="l",lwd=2,lty=3)
  lines(result[,"time"]/3600,result[,"RH"],type="l",lwd=3,lty=1,col="darkgray")
  legend(x="topright",c("cloud fraction","albedo","RH"),lwd=c(3,2,3),lty=c(1,3,1),
         col=c("black","black","darkgray"))
  dev.copy(png,"cloud_albedo_RH.png");dev.off();print("cloud_albedo_RH.png written out")
} #if(cloudTF){
