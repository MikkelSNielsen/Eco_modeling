library(deSolve)
library(fields)

# Creates an empty list to add the parameters
param <- list()

# List of parameters 
param$I0 <- 350 # incoming light radiation at surface level (micro mol photons / (m2*s); micro-Einstein)
param$kw <- 0.0375 # light absorbed by water (1/m)
param$kc <- 0.05 # light absorbed by phyto plankton (m2/mol * N)
param$gmax <- 1.5 # maximum growth (1/day) 
param$alpha <- 0.1 / (param$I0 / 200)# Slope of PI-curve (1 / micro-Einstein*d)
param$HN <- 0.3 # Half saturation constant for nutrients (mmol *  N /m^3)
param$m <- 0.03 # specific loss rate (d^-1)
param$gamma <- 0.1 # grazing (m^3*mmol^-1*day^-1) 
param$eps <- 0.1 # remenirilazation (d^-1)  
param$Dv <-  5 # Diffusivity (m^2/d)
param$up <- 0.5 # Settling velocity phytoplankton (m/d) 
param$ud <- 15 # Settling velcity detritus (m/d)
param$Nb <- 30 # nutrient content at the bottom (mmol * N/m^3) 
param$dz <- 2.5 # grid space (m)
param$z <- seq(param$dz/2, 100, by = param$dz) # depth (m)
param$n <- length(param$z) # number of grid cells

# Define initial conditions - concentrations of phytoplankton, nutrients and detritus
PND <- c(rep(10,param$n),rep(param$Nb, param$n), rep(0, param$n))


# Create a function
CalLight <- function(P, D, t, param) {
  
  seasons <- cos(2*pi/365*t)+1
  
  Q <- param$kc * param$dz * ((cumsum(P) - P/2) + (cumsum(D) - D/2))
  
  I <- param$I0 * exp(-param$kw * param$z - Q) * seasons
  
  return(I)
}


# Create function that creates the differential equation at each step with initial conditions
FuncPND <- function(t, PND, param) {
  P <- PND[1:param$n]
  N <- PND[(1+param$n):(2*param$n)]
  D <- PND[(2*param$n+1):(3*param$n)]
  
  # Phytoplankton flux
  Jap <- rep(0,param$n+1)
  Jdp <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jap[i] <- param$up* P[i-1]
    Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jap[1] = 0
  Jap[param$n+1] = 0
  
  # Diffusive flux boundary
  Jdp[1] = 0
  Jdp[param$n+1] = 0
  
  Jp = Jap + Jdp
  
  # Nutrient flux
  Jdn <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
  }

  
  # Diffusive flux boundary
  Jdn[1] = 0
  Jdn[param$n+1] = -param$Dv * (param$Nb - N[param$n]) / param$dz
  
  Jn = Jdn
  
  
  # Detritus flux
  Jad <- rep(0,param$n+1)
  Jdd <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jad[i] <- param$ud * D[i-1]
    Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jad[1] = 0
  Jad[param$n+1] =  param$ud * D[param$n]
  
  # Diffusive flux boundary
  Jdd[1] = 0
  Jdd[param$n+1] = 0
  
  Jd = Jad + Jdd
  
  # Light function
  I <- CalLight(P, D, t, param)
  
  #Growth rquation
  g <- param$gmax * pmin(param$alpha * I / sqrt((param$alpha * I)**2 + param$gmax**2), N / (N + param$HN))
  
  #Phytoplankton equation
  dPdt <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g - param$m - param$gamma * P) * P
  
  #Nutrients equation
  dNdt <- -g * P + param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
  
  #Detritus equation
  dDdt <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P - param$eps * D
  
  
  return(list(c(dPdt,dNdt, dDdt)))
  }

# Time 
time <- seq(0, 365*4, by = 1)

# Solve the differential equations
res <- ode(PND, time, FuncPND, param)
P <- res[,2:(param$n+1)]
N <- res[,(param$n+2):(param$n*2+1)]
D <- res[,(param$n*2+2):(param$n*3+1)]

 # Surface plot of phytoplankton, yearly change
 par(mfrow=c(1,3))
 image.plot(time/365, param$z, P, ylim = rev(range(param$z)), col = rainbow(400),
           xlab = "Time [Years]", ylab = "Depth [m]",
           main = "Phytoplankton concentration [cells/m^3]")
 box()
 # Add a text box to the first plot
 #text(max(par("usr")[2]), max(par("usr")[2]), "diffusivity = 9", pos=2, cex=1.2)

  # Surface plot of Nutrients,  yearly change
 image.plot(time/365, param$z, N, ylim = rev(range(param$z)), col = rainbow(400),
           xlab = "Time [Years]", ylab = "",
           main = "Nutrient concentration [mmol/m^3]")
 box()
 # Add a text box to the first plot
 #text(max(par("usr")[2]), max(par("usr")[2]), "diffusivity = 9", pos=2, cex=1.2)

  # Surface plot of Detritus,  yearly change
 image.plot(time/365, param$z, D, ylim = rev(range(param$z)), col = rainbow(400),
           xlab = "Time [Years]", ylab = "",
           main = "Detritus [mmol/m^3]")
 box()
 # Add a text box to the first plot
 #text(max(par("usr")[2]), max(par("usr")[2]), "diffusivity = 9", pos=2, cex=1.2)

 # Vertical plot of phytoplankton
 #par(mfrow=c(1,2))
 #plot(P[1,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 0 days", ylab = "Depth", xlab = "Concentration [cells/m^3]", type ="l", lwd = 2.5)
 #plot(P[365,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 1 year", ylab = "", xlab = "Concentration [cells/m^3]",type ="l", lwd = 2.5)
 #plot(P[365*2,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 2 years", ylab = "", xlab = "Concentration [cells/m^3]",type ="l", lwd = 2.5)
 #plot(P[365*3,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 3 years, Light absorption = 0.2 [1/m]", ylab = "Depth", xlab = "Concentration [cells/m^3]",type ="l", lwd = 2.5, cex.main = 0.9)
 #plot(P[365*4,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 4 years", ylab = "", xlab = "Concentration [cells/m^3]",type ="l", lwd = 2.5)
 #plot(P[365*5,], param$z, ylim = rev(range(param$z)), col = "darkolivegreen", main = "Phytoplankton: 5 yars", ylab = "", xlab = "Concentration [cells/m^3]",type ="l", lwd = 2.5)

 # Vertical plot of nutrients
 # plot(N[1,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 0 days", ylab = "Depth", xlab = "Concentration [mmol/m^3]", type ="l", lwd = 2.5)
 # plot(N[365,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 1 year", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(N[365*2,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 2 years", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(N[365*3,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 3 years", ylab = "Depth", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(N[365*4,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 4 years", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(N[365*5,], param$z, ylim = rev(range(param$z)), pch = 19, col = "blue", main = "Nutrients: 5 yars", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # 
 # # Vertical plot of detritus
 # plot(D[1,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 0 days", ylab = "Depth", xlab = "Concentration [mmol/m^3]", type ="l", lwd = 2.5)
 # plot(D[365,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 1 year", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(D[365*2,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 2 years", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(D[365*3,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 3 years", ylab = "Depth", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(D[365*4,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 4 years", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # plot(D[365*5,], param$z, ylim = rev(range(param$z)), pch = 19, col = "chocolate4", main = "Detritus: 5 yars", ylab = "", xlab = "Concentration [mmol/m^3]",type ="l", lwd = 2.5)
 # 

 
