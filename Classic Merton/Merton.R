library(readxl)
library(pracma)
library(xlsx)
library(scales)

################################################################################
################################################################################
## This code is a routine used to compute default probabilities of an issuer ###
####### using his debt and equity value in the Merton Model's framework ########
### We use Vassalou (2003) algorithm to estimate the parameters of the model ###
################################################################################
################################################################################

# Import the data
data_equity = read_excel("Data issuers.xlsx",sheet = "Mod Market Cap")
data_debt = read_excel("Data issuers.xlsx",sheet = "Gross Debt")

######### Définition des variables #############
r <- 0 # Risk-free rate
Time <- 0.5 # Time
convergence_threshold <- 0.00001 # Convergence threshold used for Vassalou algorithm
default_probabilities <- matrix(1,1,ncol(data_debt)) # Receives default probability
distances_to_default <- matrix(1,1,ncol(data_debt)) # Receives distance to default
parameters_issuer <- matrix(1,2,ncol(data_debt)) # Receives calibrated parameters


Theoretical_equity_value  <- function(x){
  ###### This function computes the theoretical value of the equity using asset value and asset value volatility #####
  d1 = (log(x/K)+(r+(vol_actif^2)/2)*t)/(vol_actif*sqrt(t))
  d2 = d1 - vol_actif*sqrt(t)
  Theoretical_equity_value = x*pnorm(d1)-K*exp(-r*t)*pnorm(d2)-equity_value
}


Distance_to_default  <- function(K,actif_value,vol_actif,r,t){
  ###### This function computes the theoretical distance to default in Merton's framework #####
  Distance_to_default = (log(actif_value/K)+(r-(vol_actif^2)/2)*t)/(vol_actif*sqrt(t))
}

Default_probability  <- function(K,actif_value,vol_actif,r,t){
  ###### This function computes the theoretical risk neutral default probability in Merton's framework. Input : distance to default of the issuer #####
  Default_probability = 1-pnorm((log(actif_value/K)+(r-(vol_actif^2)/2)*t)/(vol_actif*sqrt(t)))
}

for(Time in 1:15){
# For each issuer we apply Vassalou's algorithm to estimate the parameters (asset value and asset value volatility)
for (k in 1:length(data_debt)){
  K <- as.double(data_debt[1,k]) # debt of the current issuer
  equity <- as.double(as.matrix(data_equity[(nrow(data_equity)-251):nrow(data_equity),k+1])) # Equity of the issuer at each time
  
  actif_value <- matrix(0,length(equity),1) # The matrix receiving asset value at each time
  vol_equity <- sd(diff(log(equity)))*sqrt(252) # equity volatility
  hist_vol_actif <- append(-100,vol_equity) # our first guess for the asset volatility is the equity volatility
  
  print(k)
  cpt <- 0
  
  # Vassalou (2004) algorithm to find the asset volatility and the corresponding asset value at each time
  while (abs(hist_vol_actif[cpt+2]-hist_vol_actif[cpt+1])>convergence_threshold){
    for (i in 1:length(equity)){
      equity_value <- as.double(equity[i])
      vol_actif <- tail(hist_vol_actif,1)
      t <- Time + (252-i)/252
      # Being given the asset volatility, we solve the equation of the theoretical equity value to find the asset value of the issuer
      actif_value[i,1] = newtonRaphson(Theoretical_equity_value,K)$root
    }
    hist_vol_actif[cpt+3] <- sd(diff(log(actif_value)))*sqrt(252)
    cpt <- cpt+1
  }
  
  # We print the result for nb iterations, default probability and distance to default
  print(paste("Number of iterations for Vassalou's algorithm :",cpt,"iter",sep=" "))
  default_probabilities[1,k] <- Default_probability(K,actif_value[length(actif_value),1],tail(hist_vol_actif,1),r,t)
  distances_to_default[1,k] <- Distance_to_default(K,actif_value[length(actif_value),1],tail(hist_vol_actif,1),r,t)
  parameters_issuer[1,k] <- tail(hist_vol_actif,1)
  parameters_issuer[2,k] <- tail(actif_value,1)
  print(paste("Default probability of the issuer =",default_probabilities[1,k],sep=" "))
  print(paste("Distance to default of the issuer =",distances_to_default[1,k],sep=" "))
  
  
  # Export to Excel
  data_to_export <- as.data.frame(rbind(format(default_probabilities,scientific=FALSE),as.numeric(distances_to_default),parameters_issuer))
  colnames(data_to_export) <- colnames(data_debt)
  rownames(data_to_export) <- c(paste("Probabilité de défaut ",t,"Y",sep=""),paste("Distance to default ",t,"Y",sep=""),"Calibrated volatility", "Calibrated Asset Value")
  setwd("my_folder")
  write.xlsx(data_to_export,file=paste("DD_DP_ModMarketCap_Debt_",t,"Y.xlsx",sep=""),append=FALSE)
}
}