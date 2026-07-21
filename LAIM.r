# Land-Atmosphere Interactions Model (LAIM) 
# By John C. Lin (John.Lin@utah.edu)

require("deSolve")   #load deSolve package to access function "ode"

#################################################
# Flags to Turn On/Off Processes 
atmrespondTF <- FALSE    # does atmosphere respond to surface fluxes?
ABLTF <- FALSE           # does ABL grow or decay, according to surface heat fluxes?
cloudTF <- FALSE         # does cloud cover change as function of atmospheric humidity?
vegcontrolTF <- FALSE    # vegetation control?
soilWTF <- FALSE         # turn on soil moisture feedbacks?
co2budgetTF <- FALSE     # track atmospheric CO2, based on surface and entrainment fluxes? 
if (!atmrespondTF & ABLTF) stop ("atmrespondTF needs to be TRUE to allow ABL to grow and decay")
if (!vegcontrolTF & soilWTF) stop ("vegcontrolTF needs to be TRUE for soil moisture feedback to work")
if (co2budgetTF & !atmrespondTF) stop("for co2budgetTF to be TRUE, also requires atmrespondTF to be TRUE")
LWdnTF <- TRUE          # does LWdn respond dynamically?  
co2fluxprescTF <- FALSE # is CO2 flux (& ABL) prescribed, rather than simulated internally?
if (!co2budgetTF & co2fluxprescTF) stop ("co2budgetTF needs to be TRUE to prescribe CO2 flux")
if (!co2fluxprescTF) {if (!vegcontrolTF & co2budgetTF) stop ("vegcontrolTF needs to be TRUE to track CO2, if CO2 fluxes not prescribed")}
#################################################

#################################################
# Model timestep and duration
dt <- 60           # model timestep [s]
t.day <- 3         # run time in days
tmax <- t.day*24*3600  #maximum time [s]
times <- seq(0,tmax,dt) #vector of time steps [s]
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
k <- 0.4  # von Karman constant
#################################################

#################################################
# Load in basic functions
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
#################################################



#################################################
# Land surface characteristics
gvmax <- 1/50      # max vegetation conductance [m/s] (reciprocal of vegetation resistance) when vegcontrol is FALSE;  when TRUE, calculated by BBF function
albedo.surf <- 0.1    # surface albedo
albedo <- albedo.surf # surface albedo
z0 <- 0.5          # roughness length for momentum [m]
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
# b) heat capacity based on canopy + canopy air 
Hveg <- 10               # height of vegetation [m]
rho.veg <- 1.67          # bulk density of above-ground vegetation [kg/m3]; from Heidkamp et al. (2018): Geosci. Model Dev., 11, 3465–3479, https://doi.org/10.5194/gmd-11-3465-2018, 2018
Cp.veg <- 3000           # bulk heat capacity of above-ground vegetation [J/kg/K];  Sect. 7.2 of Bonan (2019)
Cs.veg <- Cp.veg*(rho.veg)*Hveg  # heat capacity of vegetation [J/K/m2]
rho.air <- 1.2           # air density (at sea level) [kg/m3]
Cs.air <- Cp*rho.air*Hveg
Cs <- Cs.veg + Cs.air

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
tau.soil <- 24*3600 # time constant of restoring Tsoil1 to Tsoil2 & Wsoil1 to Wsoil2 [s]
# Initialize two-layer (force-restore) soil model, from de Arellano et al. (2015)
Tsoil2 <- 286       # T of deep soil layer [K] that is constant
Tsoil1 <- Tsoil2    # T of top soil layer [K] that varies w/ time
Wsoil2 <- Wfc       # volumetric water content of deep soil layer [m3/m3]
Wsoil1 <- Wsoil2    # volumetric water content of top soil layer [m3/m3]
d1 <- 0.1           # soil depth to which diurnal variations in moisture penetrates [m], from Deardorff [1977]
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
#################################################

#################################################
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
if(!atmrespondTF)hmin <- 1000    # ABL ht matters implicitly even if atmrespondTF=FALSE, since it affects stability calculations.  Set to daytime value
thetavM0<-(Ta.c[1]+273.15)*(1+0.61*qa.presc) # initial virtual potential temperature [K]; Eq. 1.5.1b of Stull [1988]
Beta <- 0.2       # closure hypothesis:  fraction of surface virtual potential temperature flux that determines entrainment heat flux
gamma <- 5/1000   # slope of thetav above growing ABL [K/m]
qabove <- qa.presc/5  # specific humidity of air above ABL [g/g]
W <- 0            # subsidence rate [m/s]
Ur <- 1           # reference windspeed [m/s] at top of surface layer zsl
Cair <- 400       # atmospheric CO2 concentration [umole/mole, or ppm]; this is also the initial CO2 value within ABL if co2budgetTF = TRUE
Cfree <- 400      # CO2 concentration [ppm] in free troposphere (not modified by values in ABL)
Cabove <- Cfree   # CO2 concentration [ppm] above ABL (later modified by value in residual layer)
albedo.cloud <- 0.5 # albedo of cloud
#################################################

#################################################
# parameters determining CO2 greenhouse effect
CO2.SENSITIVITY <- 3.7  # CO2 doubling sensitivity [W/m2 per doubling of CO2] (IPCC 2007; Myhre et al. 1998)
CO2.baseline <- 280     # baseline to determine doubling (pre-industrial CO2 concentration [ppm])
# ave CO2 in atmospheric column, using scale height as weighting (i.e., density follows exponential decay)
CO2.colave <- Cair + (Cfree - Cair)*exp(-hmin/Hscale)           
GHG.FORCE <- CO2.SENSITIVITY*log(CO2.colave/CO2.baseline)/log(2) # GHG forcing--from CO2 elevated above CO2base.ppm [W/m2]
#################################################

#################################################
# prescribe ABL depth or CO2 fluxes
hini <- hmin  # initial ABL depth [m]
ABLdepth_DAY <- NULL
if(!ABLTF){
  print("Prescribing ABL depth...")
  # NOTE:  ABL depth needs to be prescribed with either hourly or 1-sec timestep
  ABLdepth_DAY <- SWdn_DAY
  ABLdepth_DAY[1:length(ABLdepth_DAY)] <- hmin         # prescribe daily cycle of ABL depth [m]
  hini <- ABLdepth_DAY[1]
} # if(!ABLTF){

CO2flux.veg_DAY <- NULL
if(co2fluxprescTF){
  print("Prescribing CO2 flux...")
  # NOTE:  CO2 flux needs to be prescribed with either hourly or 1-sec timestep
  CO2flux.veg_DAY <- SWdn_DAY
  CO2flux.veg_DAY[1:length(CO2flux.veg_DAY)] <- 5      # prescribe daily cycle of CO2 flux [umole/m2/s]
} # if(co2fluxprescTF){


#################################################
# Load in functions required for Monin-Obukhov stability calculations

# stability functions for momentum, based on CLASS model (https://github.com/classmodel/modelgui/blob/master/model.cpp)
psiM.f <- function(zeta){
  if(zeta <= 0){
    #unstable conditions:  from Paulson (1970) "The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer"
    x <- (1 - 16*(zeta))^0.25
    psiM <- pi/2 - 2*atan(x) + log(((1+x)^2)*(1+x^2) /8)
  } else {
    #stable conditions: from Beljaars & Holtslag (1991) “Flux Parameterization over Land Surfaces for Atmospheric Models”
    psiM <- (-2/3)*(zeta - 5/0.35)*exp(-0.35*zeta) - zeta - (10/3)/0.35
  } # if(zeta <= 0){
  return(psiM)
} # psiM.f <- function(zeta){

# stability functions for heat, based on CLASS model (https://github.com/classmodel/modelgui/blob/master/model.cpp)
psiH.f <- function(zeta){
  if(zeta <= 0){
    #unstable conditions:  from Paulson (1970) "The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer"
    x <- (1 - 16*(zeta))^0.25
    psiH <- 2*log((1+x^2)/2)
  } else {
    #stable conditions: from Beljaars & Holtslag (1991) “Flux Parameterization over Land Surfaces for Atmospheric Models”
    psiH <- -2/3 * (zeta - 5/0.35) * exp(-0.35 * zeta) -
      (1 + (2/3) * zeta)^1.5 - (10/3) / 0.35 + 1
  } # if(zeta <= 0){
  return(psiH)
} # psiH.f <- function(zeta){

clip <- function(x, xmin, xmax){
  return(pmin(pmax(x, xmin), xmax))
} # clip <- function(x, xmin, xmax){

finite_or <- function(x, fallback){
  if(!is.finite(x)) return(fallback)
  return(x)
} # finite_or <- function(x, fallback){

psiM.safe <- function(zeta, zeta_min = -5, zeta_max = 2){
  zeta <- finite_or(zeta, 0)
  zeta <- clip(zeta, zeta_min, zeta_max)
  return(psiM.f(zeta))
} # psiM.safe <- function(zeta, zeta_min = -5, zeta_max = 2){

psiH.safe <- function(zeta, zeta_min = -5, zeta_max = 2){
  zeta <- finite_or(zeta, 0)
  zeta <- clip(zeta, zeta_min, zeta_max)
  return(psiH.f(zeta))
} # psiH.safe <- function(zeta, zeta_min = -5, zeta_max = 2){

# calculate transfer-coefficients, based on Monin-Obukhov Similarity Theory (MOST)
most_transfer.f <- function(z0, Ur, zref, zeta, k = 0.4, Umin = 0.1, zeta_min = -5, zeta_max = 2,
                            raero_min = 2, raero_max = 5000){
  # roughness length for heat; from Garratt, J.R. (1978) Quart. J. Roy. Met. Soc. 104, 491-50
  z0H <- z0 * exp(-2.5)
  
  # Use an effective wind speed to avoid zero-conductance singularities.
  Ur.eff <- max(abs(Ur), Umin)
  
  # Ensure valid log-layer geometry.
  zref <- max(zref, 1.05 * max(z0, z0H))
  
  zeta <- finite_or(zeta, 0)
  zeta <- clip(zeta, zeta_min, zeta_max)
  
  # Because zeta = zref / L, roughness-level arguments are:
  # z0 / L  = (z0 / zref)  * zeta
  # z0H / L = (z0H / zref) * zeta
  zeta0  <- clip((z0  / zref) * zeta, zeta_min, zeta_max)
  zeta0H <- clip((z0H / zref) * zeta, zeta_min, zeta_max)
  
  denomM <- log(zref / z0) -
    psiM.safe(zeta, zeta_min, zeta_max) +
    psiM.safe(zeta0, zeta_min, zeta_max)
  
  denomH <- log(zref / z0H) -
    psiH.safe(zeta, zeta_min, zeta_max) +
    psiH.safe(zeta0H, zeta_min, zeta_max)
  
  # Prevent pathological denominator behavior.
  denomM <- max(denomM, 1e-6)
  denomH <- max(denomH, 1e-6)
  
  CD <- k^2 / denomM^2           # CD is drag coefficient for momentum
  CH <- k^2 / (denomM * denomH)  # CH is drag coefficient for heat (also relevant for other scalars like H2O & CO2)
  
  ustar <- k * Ur.eff / denomM   # friction velocity [m/s]
  
  raero <- 1 / (CH * Ur.eff)     # aerodynamic resistance [s/m]
  raero <- clip(raero, raero_min, raero_max)
  
  return(list(zeta = zeta,  zref = zref,
              Ur.eff = Ur.eff, z0H = z0H,  
              denomM = denomM, denomH = denomH,
              CD = CD, CH = CH, ustar = ustar, raero = raero))
} # most_transfer.f <- function(z0, Ur, zref, zeta, k = 0.4,

# calculate surface exchange, based on Monin-Obukhov Similarity Theory (MOST)
surface_exchange_most.f <- function(T, Ta, qa, qsat, rveg,
                                    Ur, zref, z0, rho, Cp, Lv,
                                    thetav, g = 9.80665,
                                    k = 0.4, Umin = 0.1,
                                    zeta_min = -5, zeta_max = 2,
                                    B0_neutral = 1e-7, allow_dew = FALSE){
  
  # Temperature and humidity differences driving surface fluxes
  dT <- T - Ta
  dq <- qsat - qa
  
  if(!allow_dew)  dq <- max(dq, 0)
  
  calc_at_zeta <- function(zeta){

    # use trial zeta to calculate transfer coefficients
    tr <- most_transfer.f(
      z0 = z0, Ur = Ur,
      zref = zref, zeta = zeta,
      k = k, Umin = Umin,
      zeta_min = zeta_min, zeta_max = zeta_max)
    
    raero <- tr$raero
    
    # sensible heat flux [W/m2]
    H <- rho * Cp * dT / raero
    
    # latent heat flux [W/m2]
    LE <- Lv * rho * dq / (raero + rveg)
    if(!allow_dew){
      LE <- max(LE, 0)
    } # if(!allow_dew){
    
    # Moisture flux in [kg/kg * m/s].
    wq <- LE / (Lv * rho)
    
    # Virtual potential temperature flux approximation
    wthetav <- H / (rho * Cp) + 0.61 * Ta * wq
    
    B0 <- g * wthetav / thetav   # surface buoyancy flux [m2/s3]
    
    if(!is.finite(B0) || abs(B0) < B0_neutral){
      L <- Inf
      zeta.new <- 0
    } else {
      L <- -tr$ustar^3 / (k * B0)
      zeta.new <- tr$zref / L
      zeta.new <- clip(zeta.new, zeta_min, zeta_max)
    } # if(!is.finite(B0) || abs(B0) < B0_neutral){
    
    # compare assumed zeta to implied zeta
    residual <- zeta - zeta.new
    
    return(list(zeta = zeta, zeta.new = zeta.new,
      residual = residual, L = L, B0 = B0,
      H = H, LE = LE, raero = raero,
      ustar = tr$ustar,
      CD = tr$CD, CH = tr$CH,
      denomM = tr$denomM, denomH = tr$denomH))
  } # calc_at_zeta <- function(zeta){
  
  neutral <- calc_at_zeta(0)
  
  if(abs(neutral$B0) < B0_neutral){
    neutral$converged <- TRUE
    neutral$method <- "neutral"
    return(neutral)
  } # if(abs(neutral$B0) < B0_neutral){
  
  # Choose a physically consistent bracket.
  # Positive buoyancy flux -> unstable -> zeta < 0.
  # Negative buoyancy flux -> stable   -> zeta > 0.
  if(neutral$B0 > 0){
    bracket <- c(zeta_min, 0)
  } else {
    bracket <- c(0, zeta_max)
  } # if(neutral$B0 > 0){
  
  residual.f <- function(zz){
    calc_at_zeta(zz)$residual
  } # residual.f <- function(zz){
  
  f1 <- residual.f(bracket[1])
  f2 <- residual.f(bracket[2])
  
  # find the stability parameter zeta that makes the MOST calculation self-consistent
  if(is.finite(f1) && is.finite(f2) && f1 * f2 <= 0){
    root <- uniroot(residual.f, lower = bracket[1], upper = bracket[2])
    out <- calc_at_zeta(root$root)
    out$converged <- TRUE
    out$method <- "uniroot"
  } else {
    # Fallback: do not crash the ODE solver.
    # Pick the zeta in the allowed range that minimizes the residual.
    opt <- optimize(f = function(zz) abs(residual.f(zz)), interval = bracket)
    out <- calc_at_zeta(opt$minimum)
    out$converged <- FALSE
    out$method <- "bounded_optimize_fallback"
  } # if(is.finite(f1) && is.finite(f2) && f1 * f2 <= 0){
  
  return(out)
} # surface_exchange_most.f <- function(T, Ta, qa, qsat, rveg,
#################################################


#################################################
# function to initialize T with equilibrium value (determined through "uniroot")
#################################################
f <- function(T, Ta, SWdn, LWdn, albedo.cloud, albedo.surf, epsilon.s, Tsoil1, Ur,
              zsl, z0, gvmax=gvmax, RH=RH, qa=qa.presc, CO2=Cair, Cfree=Cfree,
              Psurf=1000, Hscale=8000, h=hini,
              CO2.baseline=CO2.baseline, CO2.SENSITIVITY=CO2.SENSITIVITY,
              Wsoil1=Wsoil1, Wfc=Wfc, Wwilt=Wwilt,
              LAI=LAI, Kb=Kb,
              return.all=FALSE){  
  
  # --------------Physical constants--------#
  Cp <- 1005.7; Cv <- 719 # heat capacities @ constant pressure & volume [J/kg/K] (Appendix 2 of Emanuel [1994])
  g <- 9.80665 # standard surface gravity [m/s2]
  Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel [1994])
  Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson [1999])
  sigma <- 5.670373E-8 # Stefan-Boltzmann constant [W/m2/K4]
  Md <- 28.97  #molar mass of dry air [g/mole]
  k <- 0.4  # von Karman constant
  # --------------Physical constants--------#
  
  # calculate RH at ABLtop and near ground surface
  # NOTE:  this mirrors the dynamic model timestep
  e <- qa*Psurf/(Rd/Rv)      # vapor pressure [hPa]
  RH <- e/(satvap(Ta - 273.15)/100)
  P.h <- Psurf*exp(-h/Hscale)
  e.h <- qa*(Rv/Rd)*P.h # vapor pressure at ABL top [hPa]
  T.h <- Ta - (g/Cp)*h  # temperature at ABL top [K], where (g/Cp) is the adiabatic lapse rate
  esat.h <- satvap(T.h-273.15)/100 # saturation vapor pressure at ABL top [hPa]
  RH.h <- e.h/esat.h    # relative humidity at ABLtop
  
  if(cloudTF){
    # diagnose cloud fraction based on Eq. 3 of Slingo [1987]:  "The development and verification of a cloud prediction scheme for the ECMWF model"
    RHcrit <- 0.8
    # use RH at the ABLtop for the cloud scheme, as in the dynamic model timestep
    tmp <- (RH.h-RHcrit)/(1-RHcrit)
    tmp[tmp<0] <- 0
    cloud <- tmp^2
    if(cloud > 1.0)cloud <- 1.0
  } else { cloud <- 0 } # if(cloudTF){
  
  albedo <- (1-cloud)*albedo.surf + cloud*albedo.cloud
  SWup <- albedo*SWdn
  
  # ave CO2 in atmospheric column, using scale height as weighting (i.e., density follows exponential decay)
  #     NOTE: ignore the variation in CO2 within shallow residual layer (represented by updated Cabove) 
  CO2.colave <- CO2 + (Cfree - CO2)*exp(-h/Hscale)                 
  GHG.FORCE <- CO2.SENSITIVITY*log(CO2.colave/CO2.baseline)/log(2) # GHG forcing--from CO2 elevated above CO2.baseline
  
  LWup <- epsilon.s*sigma*T^4   # upward longwave radiation [W/m2]
  LWdn.t <- LWdn
  
  if (LWdnTF) {
    # empirical formula of downward longwave radiation based on Yang et al. (2023): https://doi.org/10.5194/acp-23-4419-2023
    epsilon.clr <- 0.532 + 0.808*((e/Ta)^(1/3))  # clear-sky emissivity
    epsilon.all <- epsilon.clr*(1-0.201*cloud^0.796) + 0.088*(cloud^1.038)*((RH*100)^0.221)  # all-sky emissivity
    LWdn.t <- epsilon.all * sigma * Ta^4 
    LWdn.t <- LWdn.t + GHG.FORCE  # add GHG forcing
  } # if (LWdnTF) {
  
  # determine net radiation
  Rn <- SWdn - SWup + LWdn.t - LWup
  
  # determine air density near surface
  rho.surf <- Psurf*100/(Rd*T)   # surface air density [kg/m3]
  
  # determine latent heat flux variables before MOST
  # NOTE:  LE affects buoyancy, so rveg and qsat need to be known before calling surface_exchange_most.f()
  beta.W <- 1   # water stress parameter (dependent on soil moisture)
  Lv <- 1000*latentheat(T-273.15)  # latent heat of vaporization [J/kg]
  esat <- satvap(T-273.15)/100     # saturation vapor pressure [hPa]
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
    hs <- e/esat  # RH at leaf surface [.]   
    if(!is.finite(hs)) hs <- 1.0
    if(hs<0.7) hs <- hs + 0.3   #!!! quick adjustment that ensures leaf surface is not too dry...accounts for higher humidity within canopy
    cs <- CO2    # CO2 concentration at leaf surface [umole/mole]
    
    BBFout <- BBF(SW=SWdn,Tleaf.C=T-273.15,hs=hs,beta.W=beta.W,cs=cs,Psurf=Psurf)  
    gsv <- BBFout["gsv"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]  
    rho.mole <- rho.surf*1000/Md # air density [kg/m3] => molar density [moles/m3]
    gsv <- gsv/rho.mole   # [mole/m2/s] => [m/s]
    
    # guard against zero or negative stomatal conductance
    gsv <- max(as.numeric(gsv), 1e-8)
    
    rveg <- 1/gsv         # vegetation resistance [s/m]
    An <- BBFout["An"]    # Net photosynthesis [umole/m2/s]
    ci <- BBFout["ci"]    # intercellular CO2 [umole/mole]
    
  } else {
    
    rveg <- 1/gvmax
    An <- NA
    ci <- NA
    
  } # if(vegcontrolTF){
  
  # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
  scale.canopy <- (1-exp(-Kb*LAI))/Kb
  An <- An*scale.canopy
  
  gv <- (1/rveg)*scale.canopy
  
  # guard against zero or negative canopy conductance
  gv <- max(as.numeric(gv), 1e-8)
  rveg <- 1/gv
  
  # diagnose initial virtual potential temperature
  thetavM <- Ta*(1+0.61*qa)   # virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]
  
  # determine sensible and latent heat fluxes using Monin-Obukhov Similarity Theory
  # NOTE:  surface_exchange_most.f solves for raero, H, LE, ustar, L, and zeta consistently
  MOST <- surface_exchange_most.f(
    T = T, Ta = Ta,
    qa = qa, qsat = qsat,
    rveg = rveg, Ur = Ur,
    zref = zsl, z0 = z0,
    rho = rho.surf,
    Cp = Cp, Lv = Lv,
    thetav = thetavM,
    g = g, k = k,
    Umin = 0.1,
    zeta_min = -5, zeta_max = 2,
    allow_dew = FALSE)
  
  raero <- MOST$raero
  H <- MOST$H
  LE <- MOST$LE
  ustar <- MOST$ustar
  L <- MOST$L
  zeta <- MOST$zeta
  
  # determine ground heat flux from two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
  G <- Lambda * (T - Tsoil1)
  
  # this should =0 when T is at equilibrium value
  imbalance <- Rn - H - LE - G
  
  if(return.all){
    return(list(
      imbalance = as.numeric(imbalance),
      Rn = as.numeric(Rn),
      H = as.numeric(H),
      LE = as.numeric(LE),
      G = as.numeric(G),
      LWup = as.numeric(LWup),
      LWdn = as.numeric(LWdn.t),
      SWup = as.numeric(SWup),
      cloud = as.numeric(cloud),
      albedo = as.numeric(albedo),
      RH = as.numeric(RH),
      RH.h = as.numeric(RH.h),
      qsat = as.numeric(qsat),
      VPD = as.numeric(VPD),
      rveg = as.numeric(rveg),
      raero = as.numeric(raero),
      beta.W = as.numeric(beta.W),
      ustar = as.numeric(ustar),
      L = as.numeric(L),
      zeta = as.numeric(zeta),
      MOST.converged = as.numeric(MOST$converged),
      MOST.method = MOST$method))
  } # if(return.all){
  
  return(imbalance)
  
} # f <- function(T, Ta, SWdn, LWdn, albedo.cloud, albedo.surf, epsilon.s, Tsoil1, Ur,

xinterv <- Ta.c[1]+273.15+c(-50,50)  # interval over which to search for equil temperature

# Reference height for MOST.
# limit range of zsl to avoid applying surface-layer similarity too high into the ABL or too near the roughness elements
zsl <- 0.1*hini  # surface layer height [m] assumed to be 10% of ABL height
zsl <- min(c(zsl, 100))
zsl <- max(c(zsl, 1.05 * z0))

# use initial radiation, temps to solve for initial equil. temperature
Tinit <- uniroot(
  f,  interval = xinterv,
  Ta = Ta.c[1]+273.15,
  SWdn = SWdn[1], LWdn = LWdn[1],
  Tsoil1 = Tsoil1,
  albedo.cloud = albedo.cloud,
  albedo.surf = albedo.surf,
  epsilon.s = epsilon.s,
  Ur = Ur, zsl = zsl,
  z0 = z0, gvmax = gvmax,
  RH = RH,  qa = qa.presc,
  CO2 = Cair,  Cfree = Cfree,
  Psurf = Psurf,  Hscale = Hscale,
  h = hini,
  CO2.baseline = CO2.baseline,
  CO2.SENSITIVITY = CO2.SENSITIVITY,
  Wsoil1 = Wsoil1,  Wfc = Wfc,  Wwilt = Wwilt,
  LAI = LAI,  Kb = Kb)$root

init.out <- f(
  T = Tinit,
  Ta = Ta.c[1]+273.15,
  SWdn = SWdn[1],  LWdn = LWdn[1],
  Tsoil1 = Tsoil1,
  albedo.cloud = albedo.cloud,
  albedo.surf = albedo.surf,
  epsilon.s = epsilon.s,
  Ur = Ur,
  zsl = zsl, z0 = z0,
  gvmax = gvmax,
  RH = RH,  qa = qa.presc,
  CO2 = Cair,  Cfree = Cfree,
  Psurf = Psurf,  Hscale = Hscale,
  h = hini,
  CO2.baseline = CO2.baseline,
  CO2.SENSITIVITY = CO2.SENSITIVITY,
  Wsoil1 = Wsoil1,  Wfc = Wfc,  Wwilt = Wwilt,  
  LAI = LAI,  Kb = Kb,
  return.all = TRUE)

imbalance <- init.out$imbalance

print(paste("Tinit [oC]:",signif(Tinit-273.15,5),
            ";   (Rn-H-LE-G) =",signif(imbalance,4),"[W/m2]",
            ";   raero =",signif(init.out$raero,4),"[s/m]",
            ";   ustar =",signif(init.out$ustar,4),"[m/s]",
            ";   zeta =",signif(init.out$zeta,4),
            ";   MOST.converged =",init.out$MOST.converged))

#############################################################################################################
#--------------------------------- ODE ODE ODE ODE ODE ODE ODE ODE ODE ODE ---------------------------------#
#############################################################################################################
# initialize state variables
thetaM <- Ta.c[1]+273.15
qa <- qa.presc   # initialize with prescribed specific humidity [g/g]
thetavM <- thetaM*(1+0.61*qa)   # virtual potential temperature [K];  Eq. 1.5.1b of Stull [1988]

zeta <- 0  #initialize zsl/L to 0 (neutral conditions)
yini <- c(T=Tinit, Ta=Ta.c[1]+273.15, qa=qa, thetavM=thetavM,
          Tsoil1=Tsoil1, Wsoil1=Wsoil1, h=hini, CO2=Cair) 
names(yini) <- c("T","Ta","qa","thetavM","Tsoil1","Wsoil1","h","CO2")

########################################################
# initialize parameters
# 0. numerical parameters
parms <- c(dt=dt)
# 1.  flags
parms <- c(parms,vegcontrolTF=vegcontrolTF,atmrespondTF=atmrespondTF,ABLTF=ABLTF,
           soilWTF=soilWTF,co2budgetTF=co2budgetTF,cloudTF=cloudTF,LWdnTF=LWdnTF,co2fluxprescTF=co2fluxprescTF)
# 2.  atmospheric conditions 
parms <- c(parms,Psurf=Psurf,qa.presc=qa.presc,Hscale=Hscale,hmin=hmin,Beta=Beta,
           gamma=gamma,qabove=qabove,W=W,Ur=Ur,Cabove=Cabove,Cfree=Cfree,albedo.cloud=albedo.cloud,
           CO2.baseline=CO2.baseline,CO2.SENSITIVITY=CO2.SENSITIVITY)
# 3.  land surface characteristics
parms <- c(parms,gvmax=gvmax,albedo.surf=albedo.surf,z0=z0,epsilon.s=epsilon.s,
           LAI=LAI,Kb=Kb,Hveg=Hveg,rho.veg=rho.veg,Cp.veg=Cp.veg,Cs=Cs,Resp25=Resp25,Q10=Q10)
# 4.  soil characteristics
parms <- c(parms,Wsat=Wsat,Wfc=Wfc,Wwilt=Wwilt,aa=aa,bb=bb,pp=pp,rTsoil.sat=rTsoil.sat,
           C1sat=C1sat,C2ref=C2ref,Lambda=Lambda,Tsoil2=Tsoil2,Wsoil2=Wsoil2,d1=d1,tau.soil=tau.soil)

########################################################
# define LAIM model function (what happens each time step)
LAIM <-function(time,state,parms,SWdn_DAY,LWdn_DAY,Ta.c_DAY,ABLdepth_DAY=NULL,CO2flux.veg_DAY=NULL) {
  #------------------#
  # Physical constants
  Cp <- 1005.7;Cv <- 719 # specific heat capacities of dry air @ constant pressure & volume [J/kg/K] 
  g <- 9.80665 # standard surface gravity [m/s2]
  Rd <- 287.04 # Ideal Gas Constant of DRY air [J/kg/K] (Appendix 2 of Emanuel (1994))
  Rv <- 461.40 # Ideal Gas Constant of water vapor [J/kg/K] (Appendix A.1.4 of Jacobson (1999)
  sigma <- 5.670373E-8    # Stefan-Boltzmann constant [W/m2/K4]
  Md <- 28.97  #molar mass of dry air [g/mole]
  rho.W <- 1000 # density of water [kg/m3]
  k <- 0.4  # von Karman constant
  #------------------#
  
  if (atmrespondTF & !ABLTF & is.null(ABLdepth_DAY)) {
    stop("ABLdepth_DAY must be supplied when atmrespondTF = TRUE and ABLTF = FALSE")
  } # if (atmrespondTF & !ABLTF & is.null(ABLdepth_DAY)) {
  
  # print model time, with code to deal with fact that numerical method could call LAIM() multiple instances to evaluate derivatives
  hr <- round(time / 3600, 10)
  if (abs(hr - round(hr)) < 1e-8 && round(hr) > last_printed_hour) {
    print(paste("Running model: time=", round(hr), "[hr]"))
    last_printed_hour <<- round(hr)  # "<<-" is a super-assignment: updated as a global variable
  } # if (abs(hr - round(hr)) < 1e-8 && round(hr) > last_printed_hour) {
  
  with(as.list(c(state,parms)),{
    if (!atmrespondTF) {
      Ta <- approx(x=as.numeric(names(Ta.c_DAY))*3600,y=Ta.c_DAY,xout=time%%(24*3600))$y+273.15  #use prescribed value
      qa <- qa.presc
    } # if(atmrespondTF){
    
    # calculate RH at ABLtop and near ground surface
    e <- qa*Psurf/(Rd/Rv)      # vapor pressure [hPa]
    RH <- e/(satvap(Ta - 273.15)/100)
    P.h <- Psurf*exp(-h/Hscale)
    e.h <- qa*(Rv/Rd)*P.h # vapor pressure at ABL top [hPa]
    T.h <- Ta - (g/Cp)*h  # temperature at ABL top [K], where (g/Cp) is the adiabatic lapse rate
    esat.h <- satvap(T.h-273.15)/100 # saturation vapor pressure at ABL top [hPa]
    RH.h <- e.h/esat.h    # relative humidity at ABLtop
    if(cloudTF){
      # diagnose cloud fraction based on Eq. 3 of Slingo [1987]:  "The development and verification of a cloud prediction scheme for the ECMWF model"
      RHcrit <- 0.8
      tmp <- (RH.h-RHcrit)/(1-RHcrit)
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
    
    # Reference height for MOST.
    # limit range of zsl to avoid applying surface-layer similarity too high into the ABL or too near the roughness elements
    zsl <- 0.1*h  # surface layer height [m] assumed to be 10% of ABL height
    zsl <- min(c(zsl, 100))
    zsl <- max(c(zsl, 1.05 * z0))
    
    rho.surf <- Psurf * 100 / (Rd * T)
    
    # Prepare humidity quantities before solving MOST.
    beta.W <- 1
    Lv <- 1000 * latentheat(T - 273.15)
    esat <- satvap(T - 273.15) / 100
    e <- qa * Psurf / (Rd / Rv)
    VPD <- 100 * (esat - e)
    qsat <- (Rd / Rv) * esat / Psurf
    
    if (vegcontrolTF) {
      
      if (soilWTF) {
        # Eq. (12.56) of Bonan (2019)
        beta.W <- (Wsoil1 - Wwilt) / (Wfc - Wwilt)   
        if (Wsoil1 >= Wfc) beta.W <- 1.0
        if (Wsoil1 <= Wwilt) beta.W <- 0
      } # if (soilWTF) {
      
      hs <- e / esat
      if(hs < 0.7) hs <- hs + 0.3   #!!! quick adjustment that ensures leaf surface is not too dry...accounts for higher humidity within canopy
      cs <- CO2   # CO2 concentration at leaf surface [umole/mole]
      
      # Ball-Berry + Farquhar coupled stomatal conductance & photosynthesis model for vegetation resistance [s/m]
      BBFout <- BBF(
        SW = SWdn.t, Tleaf.C = T - 273.15,
        hs = hs, beta.W = beta.W,
        cs = cs,Psurf = Psurf)
      
      gsv <- BBFout["gsv"]  # stomatal conductance with respect to water vapor [mole H2O/m2/s]
      rho.mole <- rho.surf * 1000 / Md  # air density [kg/m3] => molar density [moles/m3]
      gsv <- gsv / rho.mole # [mole/m2/s] => [m/s]
      
      # Guard against zero or negative stomatal conductance.
      gsv <- max(gsv, 1e-8)
      
      rveg <- 1 / gsv     # vegetation resistance [s/m]
      An <- BBFout["An"]  # Net photosynthesis [umole/m2/s]
      ci <- BBFout["ci"]  # intercellular CO2 [umole/mole]
      
    } else {
      rveg <- 1 / gvmax
      An <- NA
      ci <- NA
    } # if (vegcontrolTF) {
    
    # scale up photosynthesis and stomatal conductance to CANOPY values using Big-Leaf Model, based on Eq. (15.5) of Bonan (2019)
    scale.canopy <- (1 - exp(-Kb * LAI)) / Kb
    An <- An * scale.canopy
    gv <- (1 / rveg) * scale.canopy
    gv <- max(gv, 1e-8)
    rveg <- 1 / gv
    
    # Solve MOST diagnostically
    MOST <- surface_exchange_most.f(
      T = T, Ta = Ta, qa = qa, qsat = qsat,
      rveg = rveg, Ur = Ur, zref = zsl, z0 = z0,
      rho = rho.surf, Cp = Cp,  Lv = Lv,
      thetav = thetavM, g = g, k = k,
      Umin = 0.1, zeta_min = -5, zeta_max = 2)
    
    raero <- MOST$raero
    H <- MOST$H
    LE <- MOST$LE
    ustar <- MOST$ustar
    L <- MOST$L
    zeta <- MOST$zeta
    
    # determine respiration flux of CO2 to atmosphere
    Resp <- Resp25*(Q10^((T-298.15)/10))  # respiration flux based on Q10 formulation [umole CO2/m2/s]
    
    # determine ground heat flux 
    # use two-layer (force-restore) soil model to calculate ground heat flux and soil moisture
    G <- Lambda * (T - Tsoil1)
    
    Storage <- Rn - LE - H - G
    
    dT.dt <- Storage/Cs  # rate of change of T, making use of bulk heat capacity [K/s]
    
    # heat transport between surface and deep soil layer to update Tsoil1 from CLASS model (https://github.com/classmodel/modelgui/blob/master/model.cpp)
    rTsoil <- rTsoil.sat * (Wsat/Wsoil2)^(bb/(2*log(10)))
    dTsoil1.dt <- (rTsoil*G - (2*pi/tau.soil)*(Tsoil1 - Tsoil2)) #Eq. (9.32) of de Arellano et al. (2015)
    
    if (soilWTF) {
      # update soil water content, based on CLASS model (https://github.com/classmodel/modelgui/blob/master/model.cpp)
      C1 <- C1sat*(Wsat/Wsoil1)^(bb/2 + 1)  #Eq. (9.35) of de Arellano et al. (2015)
      Wsmall <- 1E-3
      C2 <- C2ref*(Wsoil2/(Wsat - Wsoil2 + Wsmall))  #Eq. (9.36) of de Arellano et al. (2015)
      Wsoil1eq <- Wsoil2 - aa*Wsat*((Wsoil2/Wsat)^pp)*(1-(Wsoil2/Wsat)^(8*pp))  #Eq. (9.37) of de Arellano et al. (2015)
      # Eq. (9.34) of de Arellano et al. (2015); NOTE:  use LE instead of LEsoil as in (9.34), and -1 multiplied by C1 that is missing in (9.34)
      dWsoil1.dt <- ((-C1/(rho.W*d1))*(LE/Lv) - (C2/tau.soil)*(Wsoil1 - Wsoil1eq))
      # make sure that Wsoil1 does not dip below Wwilt;  NOTE:  this does NOT conserve water (since could stll have residual E from minimum gv)
      if(Wsoil1 < Wwilt){dWsoil1.dt <- (Wwilt-Wsoil1)/dt;Wsoil1 <- Wwilt}  
      if (Wsoil1 < 0) {dWsoil1.dt <- (0-Wsoil1)/dt;Wsoil1 <- 0}
    } else {
      dWsoil1.dt <- 0
    } #if (soilWTF) {
    
    # if want atmosphere to respond
    # Based on "zero-order jump" or "slab" model of convective boundary layer, described in Pg. 151~155 of Garratt [1992]
    CO2flux.veg <- NA; CO2flux.ent <- NA; CO2flux.tot <- NA
    Fhthetav <- 0
    if (atmrespondTF) {
      #calculate surface virtual heat flux
      Lv <- latentheat(T-273.15)  # latent heat of vaporization [J/g]
      E <- LE/Lv   # surface moisture flux [g/m^2/s] 
      F0theta <- H/Cp  # potential heat flux [K-kg/m^2/s]
      F0thetav <- F0theta+0.073*Lv*E/Cp # virtual heat flux [K-kg/m^2/s]
      if (ABLTF) {
        Fhthetav <- -1*Beta*F0thetav   # closure hypothesis (Eq. 6.15 of Garratt [1992])
        # calculate ABL growth rate [m/s]
        dh.dt<-(1+2*Beta)*F0thetav/(rho.surf*gamma*h)  # Eq. (6.18) of Garratt [1992]
        if (F0thetav<=0.00){dh.dt <- (hmin - h)/dt;Fhthetav <- 0}   # override value:  ABL collapses
      } else {
        STEP <- 3600  # time stamp in prescribed object [s]--default is hourly
        if(max(as.numeric(names(ABLdepth_DAY)))>86000) STEP <- 1  # time stamp is in [s]
        h.t <- approx(x=as.numeric(names(ABLdepth_DAY))*STEP,y=ABLdepth_DAY,xout=time%%(24*3600))$y  
        h <- h.t
        h.tnext <- approx(x=as.numeric(names(ABLdepth_DAY))*STEP,y=ABLdepth_DAY,xout=(time+dt)%%(24*3600))$y  
        dh.dt <- (h.tnext-h.t)/dt
        if (F0thetav<=0.00|dh.dt==0){Fhthetav <- 0} # override value:  ABL collapses
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
        if(co2fluxprescTF){
          STEP <- 3600  # time step in prescribed object [s]
          if(max(as.numeric(names(CO2flux.veg_DAY)))>86000) STEP <- 1  # time stamp is in [s]
          CO2flux.veg.t <- approx(x=as.numeric(names(CO2flux.veg_DAY))*STEP,y=CO2flux.veg_DAY,xout=time%%(24*3600))$y  
          CO2flux.veg <- CO2flux.veg.t
        } # if(co2fluxprescTF){
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
    DT <- dT.dt
    DTa <- (dthetavM.dt - 0.61*Ta*dq.dt)/(1+0.61*qa)
    Dqa <- dq.dt
    DthetavM <- dthetavM.dt
    DTsoil1 <- dTsoil1.dt
    DWsoil1 <- dWsoil1.dt
    Dh <- dh.dt
    DCO2 <- dC.dt
    
    #variables that aren't integrated with time and aren't returned as derivatives
    vars2 <- c(SWdn=SWdn.t,LWdn=LWdn.t,GHG.FORCE=GHG.FORCE,Rn=as.numeric(Rn),LWup=as.numeric(LWup),H=as.numeric(H),LE=as.numeric(LE),G=as.numeric(G),
               RH=as.numeric(RH),RH.h=as.numeric(RH.h),cloud=as.numeric(cloud),albedo=as.numeric(albedo),qsat=as.numeric(qsat),
               An=as.numeric(An),Resp=as.numeric(Resp),rveg=as.numeric(rveg),raero=as.numeric(raero),beta.W=as.numeric(beta.W),
               CO2flux.veg=as.numeric(CO2flux.veg),CO2flux.ent=as.numeric(CO2flux.ent),CO2flux.tot=as.numeric(CO2flux.tot),
               dh.dt=as.numeric(dh.dt),E=as.numeric(E),Fhq=as.numeric(Fhq),deltaq=as.numeric(deltaq),Fhthetav=as.numeric(Fhthetav),
               zeta=as.numeric(zeta),ustar=as.numeric(ustar),L=as.numeric(L))
    return(list(c(DT,DTa,Dqa,DthetavM,DTsoil1,DWsoil1,Dh,DCO2),vars2))
  })
} # LAIM <-function(time,state,parms,SWdn_DAY,LWdn_DAY,Ta.c_DAY,ABLdepth_DAY=NULL,CO2flux.veg_DAY=NULL) {


########################################################
# Time integration: call LAIM function using ode()
last_printed_hour <- -Inf  # "last_printed_hour" gets updated within ode as a global variable--helps with printing out time when numerical method calls LAIM repeated times
if(atmrespondTF&ABLTF&t.day>1){
  print(paste("==========Multiple calls to ode:=========="))
  result <- NULL
  
  for (t.dd in seq_len(ceiling(max(times) / (3600 * 24)))) {
    print(paste("TIME INTEGRATION: DAY",t.dd))
    day.start <- 3600 * 24 * (t.dd - 1)
    day.end <- min(3600 * 24 * t.dd, max(times))
    times.sub <- times[times >= day.start & times <= day.end]
    
    if (t.dd > 1) {
      ilast <- nrow(result.tmp)
      #  find qa in ABL just before the ABL collapses, and use it as the humidity in residual layer that would be entrained into ABL following day
      imax <- tail(which(result.tmp$h == max(result.tmp$h, na.rm = TRUE)), 1)
      qa.resid <- result.tmp$qa[imax]
      print(paste("specific humidity of residual layer [g/g]:",signif(qa.resid,5)))
      parms["qabove"] <- qa.resid   # assign residual layer humidity as humdity above ABL
      #  find [CO2] in ABL just before the ABL collapses, and use it as the [CO2] in residual layer that would be entrained into ABL following day
      Cair.resid <- result.tmp$CO2[imax]
      print(paste("CO2 of residual layer [ppm]:",signif(Cair.resid,5)))
      parms["Cabove"] <- Cair.resid   # assign residual layer [CO2] as [CO2] above ABL
      yini <- c(
        T = result.tmp$T[ilast],
        Ta = result.tmp$Ta[ilast],
        qa = result.tmp$qa[ilast],
        thetavM = result.tmp$thetavM[ilast],
        Tsoil1 = result.tmp$Tsoil1[ilast],
        Wsoil1 = result.tmp$Wsoil1[ilast],
        h = result.tmp$h[ilast],
        CO2 = result.tmp$CO2[ilast])
    } #   if (t.dd > 1) {
    
    result.tmp <- data.frame(
      ode(yini, times.sub, LAIM, parms,
          SWdn_DAY=SWdn_DAY, LWdn_DAY=LWdn_DAY, Ta.c_DAY=Ta.c_DAY,
          ABLdepth_DAY=ABLdepth_DAY, CO2flux.veg_DAY=CO2flux.veg_DAY,
          method = "rk4")
    )
    
    result <- if (is.null(result)) result.tmp else rbind(result, result.tmp[-1, ])
  } # for (t.dd in seq_len(ceiling(max(times) / (3600 * 24)))) {
  filenm <- "result.csv"; write.csv(result,file=filenm)
  print(paste(filenm,"written out"))
} else {
  # single call "ode" to integrate LAIM model in time
  result <- ode(yini, times, LAIM, parms, SWdn_DAY=SWdn_DAY, LWdn_DAY=LWdn_DAY,Ta.c_DAY=Ta.c_DAY, 
                ABLdepth_DAY=ABLdepth_DAY, CO2flux.veg_DAY=CO2flux.veg_DAY, method = "rk4")
  result <- data.frame(result)
  filenm <- "result.csv"; write.csv(result,file=filenm)
  print(paste(filenm,"written out"))
} # if(atmrespondTF&ABLTF&t.day>1){


########################################################
# Plotting 
colorplotsTF <- TRUE  # generate colored plots? (esp. for surface energy fluxes)

# text on plot 
xmain <- paste("atmrespondTF=",atmrespondTF)
xmain <- paste(xmain,"  ABLTF=",ABLTF)
xmain <- paste(xmain,"  cloudTF=",cloudTF)
xmain <- paste(xmain,"\nvegcontrolTF=",vegcontrolTF)
xmain <- paste(xmain,"  soilWTF=",soilWTF)
if(co2fluxprescTF) xmain <- paste(xmain,"  co2fluxprescTF=",co2fluxprescTF)
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
if(colorplotsTF){
  matplot(result[,"time"]/3600,result[,c("Rn","LWdn","LWup","H","LE","G")],type="l",lty=c(1,3,1,1,1,1),lwd=c(3,2,3,2,2,2),
          cex.axis=1.5,cex.lab=1.5,col=c("black","black","darkgray","orange","blue","darkgreen"),xlab="Time [hr]",ylab="")
  legend(x="topright",c("Rn","LWdn","LWup","H","LE","G"),col=c("black","black","darkgray","orange","blue","darkgreen"),lty=c(1,3,1,1,1,1),lwd=c(3,2,3,2,2,2))
} else {
  cols <- c("black","darkgray","darkgray","black","black","darkgray")
  matplot(result[,"time"]/3600,result[,c("Rn","LWdn","LWup","H","LE","G")],type="l",lty=c(1,1,4,3,1,3),lwd=c(4,2,3,3,2,2),
          cex.axis=1.5,cex.lab=1.5,col=cols,xlab="Time [hr]",ylab="")
  legend(x="topright",c("Rn","LWdn","LWup","H","LE","G"),col=cols,lty=c(1,1,4,3,1,3),lwd=c(4,2,3,3,2,2),text.col=cols,cex=1.1,ncol=1,bty="o")
} # if(colorplotsTF)
mtext(text=expression(paste("Energy Fluxes [W ",m^-2,"]",sep="")),line=2.3,cex=1.4,side=2)
title(main=xmain)
dev.copy(png,"Energyfluxes.png");dev.off();print("Energyflux.png written out")

# plot soil water content 
if (soilWTF) {
  dev.new()
  plot(result[,"time"]/3600,result[,"Wsoil1"],type="l",xlab="Time [hour]",ylab="",cex.main=1.0,
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=paste("Soil type =",soiltype,"\n",xmain),col="black",ylim=c(Wwilt,Wfc))
  mtext(text=expression(paste("Soil Volumetric Water Content [",m^3,"/",m^3,"]",sep="")),line=2,cex=1.3,side=2)
  abline(h=Wsoil2,lty=3,lwd=2)
  legend(x="topright",c("Wsoil1","Wsoil2"),col=c("black","black"),lwd=2,lty=c(1,3))
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
       cex.axis=1.3,cex.lab=1.3,lwd=3,main=xmain)
  mtext(text=expression(paste("Photosynthesis or Respiration [",mu,"mole ",m^-2," ",s^-1,"]",sep="")),line=2,cex=1.3,side=2)
  lines(result[,"time"]/3600,result[,"Resp"],lty=3,lwd=2)
  legend(x="topright",c("Photosynthesis (An)","Respiration"),lwd=c(3,2),lty=c(1,3))
  dev.copy(png,"An_Resp.png");dev.off();print("An_Resp.png written out")
} #if(vegcontrolTF){

if ((vegcontrolTF|co2fluxprescTF)&atmrespondTF&co2budgetTF) {
  # plot time series of CO2 
  dev.new()
  plot(result[,"time"]/3600,result[,"CO2"],type="l",xlab="Time [hour]",ylab="",
       cex.axis=1.3,cex.lab=1.3,lwd=2,main=xmain)
  mtext(text=expression(paste("CO"[2]," [ppm]"),sep=""),line=2.5,cex=1.3,side=2)
  par(new=TRUE)
  ylims <- range(result[,c("CO2flux.ent","CO2flux.veg")],na.rm=TRUE)
  plot(result[,"time"]/3600,result[,"CO2flux.veg"],type="l",axes=F,xlab="",ylab="",col="darkgray",ylim=ylims,lty=1,lwd=2)
  lines(result[,"time"]/3600,result[,"CO2flux.ent"],col="darkgray",lty=3,lwd=2)
  abline(h=0,lty=1,lwd=0.5,col="darkgray")
  axis(4,cex.lab=1.3,cex.axis=1.3,col="darkgray",col.axis="darkgray")
  legend(x="topright",c("CO2 in ABL","CO2flux.veg","CO2flux.ent"),lwd=2,lty=c(1,1,3),
         col=c("black","darkgray","darkgray"),text.col=c("black","darkgray","darkgray"))
  mtext(text=expression(paste("CO"[2]," Flux [",mu,"mole ",m^-2," ",s^-1,"]",sep="")),line=-1,cex=1.3,side=4,col="darkgray")
  dev.copy(png,"CO2.png");dev.off();print("CO2.png written out")
} #if ((vegcontrolTF|co2fluxprescTF)&atmrespondTF&co2budgetTF) {

if (cloudTF) {
  # plot time series of cloud fraction, albedo, and relative humidity
  dev.new()
  ylims <- range(result[,c("cloud","albedo","RH")],na.rm=TRUE)
  plot(result[,"time"]/3600,result[,"cloud"],type="l",xlab="Time [hour]",ylab="Cloud Fraction/Albedo/RH",
       cex.axis=1.3,cex.lab=1.3,lwd=3,lty=1,main=xmain,ylim=ylims)
  lines(result[,"time"]/3600,result[,"albedo"],type="l",lwd=2,lty=3)
  lines(result[,"time"]/3600,result[,"RH.h"],type="l",lwd=3,lty=1,col="darkgray")
  lines(result[,"time"]/3600,result[,"RH"],type="l",lwd=3,lty=3,col="darkgray")
  legend(x="topright",c("cloud fraction","albedo","RH@ABLtop","RH near surf"),lwd=c(3,2,3,3),lty=c(1,3,1,3),
         col=c("black","black","darkgray","darkgray"))
  dev.copy(png,"cloud_albedo_RH.png");dev.off();print("cloud_albedo_RH.png written out")
} #if(cloudTF){
