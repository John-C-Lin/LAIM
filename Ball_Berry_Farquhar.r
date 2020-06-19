# Solves simultaneously for both the stomatal conductance and the photosynthetic rate
# based on the Ball-Berry (1987) and Farquhar (1980) models, respectively.
# The solution is described in pages 196~197 of Bonan (2019):  "Climate change and terrestrial ecosystem modeling"
# 6/1/2020 by John C. Lin (John.Lin@utah.edu)

BBF <- function(SW,Tleaf.C,hs,beta=1.0,cs=400,Psurf=1000){
# Arguments:  Atmospheric variables 
# SW is shortwave radiation [W/m2]
# Psurf is surface pressure [hPa]
# hs is fractional humidity (=RH/100) at leaf surface [.]
# cs is CO2 concentration at leaf surface [umole/mole]

Psurf <- Psurf*100  # [hPa]=>[Pa]

#################################################
# Stomatal conductance variables 
g1 <- 9    #slope of Ball-Berry relatinnship
g0 <- 0.01 # minimum stomatal conductance [mole H2O/m2/s]
#################################################

#################################################
# Plant physiological parameters for C3 plants 
# Bernacchi et al. (2001); From Table 11.1 of Bonan (2019)
Kc25 <- 404.9*Psurf/1E6  # Michaelis-Menten constant for CO2 @ 25oC [Pa]
Ko25 <- 278.4*Psurf/1E3  # Michaelis-Menten constant for O2 @ 25oC [Pa]
Gamma.star25 <- 42.75*Psurf/1E6 # CO2 compensation point @ 25oC [Pa]
Vcmax25 <- 60 # max carboxylation rate @ 25oC [umole/m2/s];  Table 11.4 of Bonan (2019); value of 60 is typical among different kinds of trees

Rd25 <- 0.015*Vcmax25 # leaf respiration [umole/m2/s];  Table 11.5 of Bonan (2019)
  
# electron transport parameters
THETAJ <- 0.7    # curvature parameter (von Caemmerer et al., 2009) [.]
PHI.PSII <- 0.85 # quantum yield of photosystem II (von Caemmerer et al., 2009) [mole/mole]
alpha.l <- 0.8   # leaf absorptance (section 11.2 of Bonan (2019) [.]
Jmax25 <- Vcmax25*2.0  # max electron transport rate @ 25oC (Leuning 2002) [umole/m2/s]
#################################################

#################################################
# Temperature dependence of physiological parameters 
# Arrhenius parameters from Bernacchi et al. (2001, 2003); from Table 11.2 of Bonan (2019)
DHa.Kc <- 79430  # [J/mole]
DHa.Ko <- 36380  # [J/mole]
DHa.Gamma.star <- 37830  # [J/mole]
DHa.Vcmax <- 65330  # [J/mole]
DHa.Rd <- 46390     # [J/mole]
DHa.Jmax <- 43540   # [J/mole]

# High temperature inhibition:  thermal breakdown parameters from Table 11.3 of Bonan (2019)
DHd <- 150000  # [J/mole]
DS <- 490      # [J/K/mole]
#################################################

# Apply Arrhenius function to modify physiological parameters, based on leaf T (Eq. 11.34 of Bonan (2019))
Rg <- 8.314 # universal gas constant [J/K/mole]
Tleaf.K<-Tleaf.C + 273.15
Vcmax <- Vcmax25*exp((DHa.Vcmax*(1-298.15/Tleaf.K)/(298.15*Rg)))   # [umole/m2/s]
Kc <- Kc25*exp((DHa.Kc*(1-298.15/Tleaf.K)/(298.15*Rg))) # [Pa]
Ko <- Ko25*exp((DHa.Ko*(1-298.15/Tleaf.K)/(298.15*Rg))) # [Pa]
Gamma.star <- Gamma.star25*exp((DHa.Gamma.star*(1-298.15/Tleaf.K)/(298.15*Rg)))  #[Pa]
Jmax <- Jmax25*exp((DHa.Jmax*(1-298.15/Tleaf.K)/(298.15*Rg))) # [umole/m2/s]
Rd <- Rd25*exp((DHa.Rd*(1-298.15/Tleaf.K)/(298.15*Rg)))       # [umole/m2/s]

# Further apply thermal breakdown function to represent high T inhibition of a subset of parameters (Eq. 11.36 of Bonan (2019))
Vcmax <- Vcmax*(1+exp((298.15*DS-DHd)/(298.15*Rg)))/(1+exp((DS*Tleaf.K-DHd)/(Rg*Tleaf.K)))
Jmax <- Jmax*(1+exp((298.15*DS-DHd)/(298.15*Rg)))/(1+exp((DS*Tleaf.K-DHd)/(Rg*Tleaf.K)))
Rd <- Rd*(1+exp((298.15*DS-DHd)/(298.15*Rg)))/(1+exp((DS*Tleaf.K-DHd)/(Rg*Tleaf.K)))

# Scale Vcmax and Jmax with water stress parameter "beta" (dependent on soil moisture)
Vcmax <- beta*Vcmax
Jmax <- beta*Jmax

# Calculate J, the rate of electron transport [umole/m2/s]
PAR <- 1.98*SW  # conversion from shortwave radiation [W/m2] to PAR [umole-photon/m2/s];  Lin et al. (2011)

I.PSII<-(PHI.PSII/2)*alpha.l*PAR   #amt of light used by photosystem II [umole/m2/s];  (Eq. 11.23 of Bonan (2019))
J <- ((I.PSII+Jmax)-sqrt((I.PSII+Jmax)^2-4*THETAJ*I.PSII*Jmax))/(2*THETAJ)

oi <- 209           # intercellular O2 [mmole/mole]
Kc <- Kc*1E6/Psurf  # [Pa]=>[umole/mole]
Ko <- Ko*1E3/Psurf  # [Pa]=>[mmole/mole]
Gamma.star <- Gamma.star*1E6/Psurf  # [Pa]=>[umole/mole]

# Solve for ci, the intercellular CO2 [umole/mole], subject to constraint that the demand function for CO2 (photosynthesis) = 
# supply function for CO2 (diffusion/stomatal conductance);  Eq. (12.26) of Bonan (2019): 
vv <- g1*hs/(1.6*cs)
# 1) solve quadratic eqn for ci under Ac, the Rubisco-limited assimilation rate [umole/m2/s]
aa <- Vcmax
bb <- Kc*(1 + oi/Ko)
AA <- (g0/1.6) + vv*(aa - Rd)
BB <- ((1 - vv*cs)*(aa - Rd) + (g0/1.6)*(bb - cs) - vv*(aa*Gamma.star + bb*Rd))
CC <- ((vv*cs - 1)*(aa*Gamma.star + bb*Rd) - bb*cs*g0/1.6)
ci.Ac <- (-BB + sqrt(BB^2 - 4*AA*CC))/(2*AA)
# 2) solve quadratic eqn for ci under Aj, the light-limited assimilation rate [umole/m2/s]
aa <- J/4
bb <- 2*Gamma.star
AA <- (g0/1.6) + vv*(aa - Rd)
BB <- ((1 - vv*cs)*(aa - Rd) + (g0/1.6)*(bb - cs) - vv*(aa*Gamma.star + bb*Rd))
CC <- ((vv*cs - 1)*(aa*Gamma.star + bb*Rd) - bb*cs*g0/1.6)
ci.Aj <- (-BB + sqrt(BB^2 - 4*AA*CC))/(2*AA)

# Rubisco-limited assimilation rate [umole/m2/s];  Table 11.5 of Bonan (2019)
Ac <- Vcmax*(ci.Ac - Gamma.star)/(ci.Ac + Kc*(1 + oi/Ko))  
# Light-limited assimilation rate [umole/m2/s];  Table 11.5 of Bonan (2019)
Aj <- (J/4)*((ci.Aj - Gamma.star)/(ci.Aj + 2*Gamma.star))
if (is.nan(Ac)) Ac <- 0
if (is.nan(Aj)) Aj <- 0

# tmp<-c(cs,ci.Ac,ci.Aj,Gamma.star,J);names(tmp)<-c("cs","ci.Ac","ci.Aj","Gamma.star","J")
# print(signif(tmp,5))
# print(paste("Ac:",signif(Ac,5),"; Aj:",signif(Aj,5)))

# photosynthesis rate is which ever of Ac or Aj that is limiting
# A <- apply(cbind(Ac,Aj),1,min)
if (Ac <= Aj) {A <- Ac; ci <- ci.Ac}
if (Aj < Ac)  {A <- Aj; ci <- ci.Aj}

# Net photosynthesis [umole/m2/s]
An <- A - Rd

# stomatal conductance with respect to water vapor, following Ball-Berry model [mole H2O/m2/s];  Eq. (12.14) of Bonan (2019)
gsw <- g0 + g1*(An/cs)*hs

result <- c(gsw, An, ci)
names(result) <- c("gsw","An","ci")
return(result)

} #BBF <- function(SW,Tleaf.C,hs,cs=400,Psurf=1000){
  
  

