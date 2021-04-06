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

par_mu <- function(variance,kurt){
  ###### This function computes rho as a function of  kurtosis
  par_mu = sqrt(15*variance/kurt)
}

par_lambda <- function(variance,kurt){
  ###### This function computes lambda as a function of variance and kurtosis
  par_lambda = 15/kurt*sqrt(15*variance/kurt)
}

phi <- function(x,t,lambda,mu){
  phi = pnorm(sqrt((lambda*t^2)/x)*(x/(mu*t)-1))+exp(2*lambda*t/mu)*pnorm(-sqrt((lambda*t^2)/x)*(x/(mu*t)+1))
}

Theoretical_equity_value  <- function(x){
  ###### This function computes the theoretical value of the equity using a Negative Gamma Process #####
  k_I <- log(x/K)+(r+(lambda/mu)*(sqrt(1+2*(mu^2/lambda))-1))*t
  if (k_I>0){
    # If k_G > 0 we apply equation (29)....
    Theoretical_equity_value = x*phi(k_I*sqrt(1+(2*mu^2)/lambda),t,lambda*sqrt(1+(2*mu^2)/lambda),mu)-K*exp(-r*t)*phi(k_I,t,lambda,mu)-equity_value
    # ...Else, V_E = 0
  }else{
    Theoretical_equity_value = -equity_value
  }
}

Distance_to_default  <- function(K,actif_value,lambda,mu,r,t){
  ###### This function computes the theoretical distance to default in the Negative Gamma Process framework #####
  Distance_to_default = log(actif_value/K)+(r+(lambda/mu)*(sqrt(1+2*(mu^2/lambda))-1))*t
}

Default_probability  <- function(K,actif_value,lambda,mu,r,t){
  ###### This function computes the theoretical risk neutral default probability in the Negative Gamma Process framework. Input : distance to default of the issuer ####
  k_I <- log(actif_value/K)+(r+(lambda/mu)*(sqrt(1+2*(mu^2/lambda))-1))*t
  
  if (k_I>0){
    Default_probability = pnorm(-sqrt((lambda*t^2)/k_I)*(k_I/(mu*t)-1))-exp(2*lambda*t/mu)*pnorm(-sqrt((lambda*t^2)/k_I)*(k_I/(mu*t)+1))
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
Time <- 0.5 # Time
convergence_threshold <- 0.0001 # Convergence threshold used for Vassalou algorithm
default_probabilities <- matrix(1,1,ncol(data_debt)) # Receives default probability
distances_to_default <- matrix(1,1,ncol(data_debt)) # Receives distance to default
parameters_issuer <- matrix(1,3,ncol(data_debt)) # Receives the calibrated parameters


# For each issuer we apply Vassalou's algorithm to estimate the parameters (asset value and asset value volatility)
for (k in 1:length(data_debt)){
  K <- as.double(data_debt[1,k]) # debt of the current issuer
  equity <- as.double(as.matrix(data_equity[(nrow(data_equity)-251):nrow(data_equity),k+1])) # Equity of the issuer at each time
  rV_E <- diff(log(equity)) # Equity returns
  mu_0 <- par_mu(var(rV_E)*252,kurtosis(rV_E)) # First guess for mu = equity mu
  lambda_0 <- par_lambda(var(rV_E)*252,kurtosis(rV_E))  # First guess for lambda = equity lambda
  actif_value <- matrix(0,length(equity),1) # The matrix receiving the asset value at each time
  hist_mu_actif <- append(-100,mu_0) # Receives rho at each iteration
  hist_lambda_actif <- append(-100,lambda_0) # Receives lambda at each iteration
  
  print(k)
  cpt <- 0
  
  # Vassalou (2004) algorithm to find the asset volatility and the corresponding asset value at each time, we stop when convergence
  while ((abs(hist_mu_actif[cpt+2]-hist_mu_actif[cpt+1])>convergence_threshold)|(abs(hist_lambda_actif[cpt+2]-hist_lambda_actif[cpt+1])>convergence_threshold)){
    mu <- tail(hist_mu_actif,1)
    lambda <- tail(hist_lambda_actif,1)
    
    for (i in 1:length(equity)){
      equity_value <- as.double(equity[i])
      t <- Time + (252-i)/252
      # Being given the asset mu and lambda, we solve the equation of the theoretical equity value to find the asset value of the issuer
      actif_value[i,1] <- newtonRaphson(Theoretical_equity_value,K)$root
    }
    
    rV_A <- diff(log(actif_value)) # Asset log returns
    variance <- var(rV_A)*252 # Annualized variance of asset
    kurt <- kurtosis(rV_A) # Excess Kurtosis of asset
    print(paste("Itération numéro",cpt+1,sep=" "))
    print(paste("mu:",hist_mu_actif[cpt+2],sep=" "))
    print(paste("lambda:",hist_lambda_actif[cpt+2],sep=" "))
    print(paste("V_A:",tail(actif_value,1),sep=" "))
    print(paste("Variance:",variance,sep=" "))
    
    hist_mu_actif[cpt+3] <- par_mu(variance,kurt) # Computes the new mu
    hist_lambda_actif[cpt+3] <- par_lambda(variance,kurt) # Computes the new lambda
    
    cpt <- cpt+1
  }
  
  # We print the result for nb iterations, default probability and distance to default
  default_probabilities[1,k] <- Default_probability(K,tail(actif_value,1),tail(hist_lambda_actif,1),tail(hist_mu_actif,1),r,t)
  distances_to_default[1,k] <- Distance_to_default(K,tail(actif_value,1),tail(hist_lambda_actif,1),tail(hist_mu_actif,1),r,t)
  parameters_issuer[1,k] <- tail(hist_mu_actif,1)
  parameters_issuer[2,k] <- tail(hist_lambda_actif,1)
  parameters_issuer[3,k] <- tail(actif_value,1)
  print(paste("Number of iterations for Vassalou's algorithm :",cpt,"iter",sep=" "))
  print(paste("Default probability of the issuer =",default_probabilities[1,k],sep=" "))
  print(paste("Distance to default of the issuer =",distances_to_default[1,k],sep=" "))
  print(paste("Parameters = mu :",tail(hist_mu_actif,1),"lambda :",tail(hist_lambda_actif,1),"V_A:",tail(actif_value,1),sep=" "))
}

# Export to Excel
data_to_export <- as.data.frame(rbind(format(default_probabilities,scientific=FALSE),as.numeric(distances_to_default),parameters_issuer))
colnames(data_to_export) <- colnames(data_debt)
rownames(data_to_export) <- c(paste("Probabilité de défaut ",Time,"Y",sep=""),paste("Distance to default ",Time,"Y",sep=""),"Calibrated Mu","Calibrated Lambda","Calibrated Asset Value")
setwd("my_folder")
write.xlsx(data_to_export,file=paste("DD_DP_ModMarketCap_Debt_",Time,"Y.xlsx",sep=""),append=FALSE)


print(proc.time()-ptm)


