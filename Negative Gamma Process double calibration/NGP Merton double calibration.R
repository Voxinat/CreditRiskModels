ptm <- proc.time()

library(readxl)
library(pracma)
library(xlsx)
library(scales)
library(GeneralizedHyperbolic)
library(moments)

################################################################################
################################################################################
## This code is a routine used to compute default probabilities of an issuer ###
### using his debt and equity value in the case of an extended Merton Model ####
############# where the log returns of the equity value follows ################
####################### a Negative Gamma Process ###############################
##### We fixed one of the parameter of the NGP and  use Vassalou (2003)  #######
########### algorithm to estimate the 2 other parameters of the model ##########
################################################################################
################################################################################

par_rho <- function(kurt){
  ###### This function computes rho as a function of  kurtosis
  par_rho = 6/kurt
}

par_lambda <- function(variance,rho){
  ###### This function computes lambda as a function of variance and kurtosis
  par_lambda = sqrt(rho/variance)
}

Theoretical_equity_value  <- function(x){
  ###### This function computes the theoretical value of the equity using a Negative Gamma Process #####
  k_G <- log(x/K)+(r+rho*log(1+1/lambda))*t
  if (k_G>0){
    # If k_G > 0 we apply equation (29)....
    Theoretical_equity_value = x*as.numeric(gammainc((lambda+1)*k_G,rho*t)[1])/gamma(rho*t)-K*exp(-r*t)*as.numeric(gammainc(lambda*k_G,rho*t)[1])/gamma(rho*t)-equity_value
  }else{
    # ...Else, V_E = 0
    Theoretical_equity_value = -equity_value
  }
}

Distance_to_default  <- function(K,actif_value,rho,lambda,r,t){
  ###### This function computes the theoretical distance to default in the Negative Gamma Process framework #####
  Distance_to_default = log(actif_value/K)+(r+rho*log(1+1/lambda))*t
}

Default_probability  <- function(K,actif_value,rho,lambda,r,t){
  ###### This function computes the theoretical risk neutral default probability in the Negative Gamma Process framework. Input : distance to default of the issuer ####
  k_G <- log(actif_value/K)+(r+rho*log(1+1/lambda))*t
  
  if (k_G>0){
    Default_probability = as.numeric(gammainc(lambda*k_G,rho*t)[2])/gamma(rho*t)
  }else{
    Default_probability = 1
  }
}

##################################################################
########################### Main Code ############################
##################################################################


# Import the data
data_equity = read_excel("Data issuers.xlsx",sheet = "Mod Market Cap")
data_debt = read_excel("Data issuers.xlsx",sheet = "Gross Debt")

######### Définition des variables #############
r <- 0 # Risk-free rate
Time <- 1 # Time
convergence_threshold <- 0.0001 # Convergence threshold used for Vassalou algorithm
default_probabilities <- matrix(1,1,ncol(data_debt)) # Receives default probability
distances_to_default <- matrix(1,1,ncol(data_debt)) # Receives distance to default
parameters_issuer <- matrix(1,3,ncol(data_debt)) # Receives the calibrated parameters

for(Time in 2:15){
# For each issuer we apply Vassalou's algorithm to estimate the parameters (asset value and asset value volatility)
for (k in 1:length(data_debt)){
  K <- as.double(data_debt[1,k]) # debt of the current issuer
  equity <- as.double(as.matrix(data_equity[(nrow(data_equity)-251):nrow(data_equity),k+1])) # Equity of the issuer at each time
  rV_E <- diff(log(equity)) # Equity returns
  rho_0 <- par_rho(kurtosis(rV_E)) # First guess for rho = equity rho
  lambda_0 <- par_lambda(var(rV_E)*252,rho_0) # First guess for lambda = equity lambda
  actif_value <- matrix(0,length(equity),1) # The matrix receiving the asset value at each time
  hist_rho_actif <- append(-100,rho_0) # Receives rho at each iteration
  hist_lambda_actif <- append(-100,lambda_0) # Receives lambda at each iteration
  
  print(k)
  cpt <- 0
  
  # Vassalou (2004) algorithm to find the asset volatility and the corresponding asset value at each time, we stop when convergence
  while ((abs(hist_rho_actif[cpt+2]-hist_rho_actif[cpt+1])>convergence_threshold)|(abs(hist_lambda_actif[cpt+2]-hist_lambda_actif[cpt+1])>convergence_threshold)){
    rho <- tail(hist_rho_actif,1)
    lambda <- tail(hist_lambda_actif,1)
    for (i in 1:length(equity)){
      equity_value <- as.double(equity[i])
      t <- Time + (252-i)/252
      # Being given the asset rho and lambda, we solve the equation of the theoretical equity value to find the asset value of the issuer
      actif_value[i,1] <- newtonRaphson(Theoretical_equity_value,K)$root
    }
    
    rV_A <- diff(log(actif_value)) # Asset log returns
    variance <- var(rV_A)*252 # Annualized variance of asset
    kurt <- kurtosis(rV_A) # Kurtosis of asset
    print(paste("Itération numéro",cpt+1,sep=" "))
    print(paste("rho:",hist_rho_actif[cpt+2],sep=" "))
    print(paste("lambda:",hist_lambda_actif[cpt+2],sep=" "))
    print(paste("V_A:",tail(actif_value,1),sep=" "))
    print(paste("Variance:",variance,sep=" "))
    
    hist_rho_actif[cpt+3] <- par_rho(kurt) # Computes the new rho
    hist_lambda_actif[cpt+3] <- par_lambda(variance,hist_rho_actif[cpt+3]) # Computes the new lambda
    
    cpt <- cpt+1
  }
  
  # We print the result for nb iterations, default probability and distance to default
  default_probabilities[1,k] <- Default_probability(K,tail(actif_value,1),tail(hist_rho_actif,1),tail(hist_lambda_actif,1),r,t)
  distances_to_default[1,k] <- Distance_to_default(K,tail(actif_value,1),tail(hist_rho_actif,1),tail(hist_lambda_actif,1),r,t)
  parameters_issuer[1,k] <- tail(hist_rho_actif,1)
  parameters_issuer[2,k] <- tail(hist_lambda_actif,1)
  parameters_issuer[3,k] <- tail(actif_value,1)
  print(paste("Number of iterations for Vassalou's algorithm :",cpt,"iter",sep=" "))
  print(paste("Default probability of the issuer =",default_probabilities[1,k],sep=" "))
  print(paste("Distance to default of the issuer =",distances_to_default[1,k],sep=" "))
  print(paste("Parameters = rho :",tail(hist_rho_actif,1),"lambda :",tail(hist_lambda_actif,1),"V_A:",tail(actif_value,1),sep=" "))
}

# Export to Excel
data_to_export <- as.data.frame(rbind(format(default_probabilities,scientific=FALSE),as.numeric(distances_to_default),parameters_issuer))
colnames(data_to_export) <- colnames(data_debt)
rownames(data_to_export) <- c(paste("Probabilité de défaut ",Time,"Y",sep=""),paste("Distance to default ",Time,"Y",sep=""),"Calibrated Rho","Calibrated Lambda","Calibrated Asset Value")
setwd("my_folder")
write.xlsx(data_to_export,file=paste("DD_DP_ModMarketCap_Debt_",Time,"Y.xlsx",sep=""),append=FALSE)
}

print(proc.time()-ptm)


