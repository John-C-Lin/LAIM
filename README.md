# Land-Atmosphere Interactions Model (LAIM)
Land-Atmosphere Interactions Model (LAIM)
by John Chun-Han Lin (John.Lin@utah.edu)

## Simple 1-D model of the land and the atmosphere in order to study the coupling and interactions.
The LAIM model adopts numerous simplifications/assumptions as part of its “spherical cow” approach:
- No horizontal fluxes and horizontal advection
- Variations within mixed-layer is minimal
- Variations within thin surface layer neglected
- Inversion layer in the entrainment zone represented by a “jump”

## Here are some specific steps to run LAIM in the <a href="https://www.r-project.org/">R programming language</a>:
1) Download the LAIM model from Github: https://github.com/John-C-Lin/LAIM.git <br>
2) Run the coupled land-atmosphere surface model, coded up in LAIM.r. This can be done using
the source command in R (see Chapt. 4 of "Intro to Land-Atmosphere Interactions"):
source(“PATH/LAIM.r”)  <br>
3) Make sure that the R working directory is set to the directory where both LAIM.r and
Ball_Berry_Farquhar.r are found. 
You can use the 'getwd' command in R to get the current working directory, and the 'setwd' command to set the working directory. <br>
Why?  This is because LAIM.r also reads in the Ball_Berry_Farquhar.r script, which is the model of vegetation resistance, or stomatal opening.
Ball_Berry_Farquhar.r creates a function named BBF, which is the Ball-Berry-Farquhar stomatal model that is called within LAIM.
