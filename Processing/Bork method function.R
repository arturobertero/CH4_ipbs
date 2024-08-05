#The script below uses adapted code originally based on Appendix A of van Bork et al. (2019) 
#https://doi.org/10.1080/00273171.2019.1672515
#That paper was shared via an Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) licence 
#https://creativecommons.org/licenses/by-nc-nd/4.0/
#However, Dr van Bork provided emailed approval on 13 July 2021 to share this adapted script.

#By default this script applies our version of cross-validation, although
#It can be applied in such a way as not to use cross-validation.

#Input arguments
#n is the number of simulations to apply
#train is the training data
#test is the test data
#cfa is a unidimensional factor model fit through lavaan
#ggm_mat is an omega matrix of partial correlations for a fitted network model

Bork_fun = function(n, train, test, cfa, ggm_mat){
  
  #implied correlation matrices
  diag(ggm_mat) <- 1 #Set main diagonal of partial correlation matrix from network analysis to 1s
  nV <- ncol(train) #define the number of variables
  cov_FA <- fitted(cfa)$cov #obtain the implied covariance matrix for the UFM
  cor_FA <- cov2cor(cov_FA) #obtain the implied correlation matrix for the UFM
  cor_NW <- pcor2cor(ggm_mat) #implied correlation matrix from network analysis
  cor_data <- cor(test) #observed correlation matrix (based on test data)
  pcor_data <- cor2pcor(cor_data) #observed partial correlation matrix
  
  U_data_cor <- cor_data[lower.tri(cor_data, diag = F)]
  #only the unique elements in the observed (polychoric) correlation matrix
  
  U_data_pcor <- pcor_data[lower.tri(pcor_data, diag = F)] 
  #only the unique elements in the observed partial correlation matrix
  
  signswitch_data <- which(U_data_pcor * U_data_cor < 0) 
  #identifying indices of zero-order correlations that have different sign than their corresponding  partial correlations
  
  stronger_data <- which((U_data_pcor^2) > (U_data_cor^2)) 
  #identifying partial correlations of absolute value greater than zero-order correlations
  
  together_data <- union(signswitch_data, stronger_data)
  total_data <- length(together_data) #cases that don't match what's expected with a UFM
  
  propData <- total_data/length(U_data_cor) #proportion of anomalous partial correlations in matrix based on test data
  
  #Preliminary steps for simulation
  propFA <- rep(NA, n)
  propNW <- rep(NA, n)
  nobs = nrow(test)
  
  #In the simulations below, multivariate normal data is still assumed
  #I.e., we take the categorical nature of the data into account when simulating data for our power analyses
  #And also when estimating the CFA and network models
  #But not in van Bork's simulations for checking how often particular numbers of partial correlations are observed
  #Interested readers could develop modifications of van Bork's method that don't simulate multivariate normal data.

  for(i in 1:n) {
    DataNW <- data.frame(mvrnorm(n = nobs, mu = rep(0,nV), Sigma = cor_NW, empirical=F))
    DataFA <- data.frame(mvrnorm(n=nobs, mu = rep(0,nV), Sigma = cor_FA, empirical=F))
    
    corFA <- cor(DataFA) 
    pcorFA <- cor2pcor(corFA)
    corNW <- cor(DataNW) 
    pcorNW <- cor2pcor(corNW)
    U_FA_cor <- corFA[lower.tri(corFA, diag = F)]
    U_FA_pcor <- pcorFA[lower.tri(pcorFA, diag = F)]
    U_NW_cor <- corNW[lower.tri(corNW, diag = F)]
    U_NW_pcor <- pcorNW[lower.tri(pcorNW, diag = F)]
    signswitch_FA <- which(U_FA_pcor * U_FA_cor < 0)
    signswitch_NW <- which(U_NW_pcor * U_NW_cor < 0)
    stronger_FA <- which((U_FA_pcor ^ 2) >(U_FA_cor^2))
    stronger_NW <- which((U_NW_pcor ^ 2) > (U_NW_cor^2))
    together_FA <- union(signswitch_FA, stronger_FA)
    together_NW <- union(signswitch_NW, stronger_NW)
    total_FA <- length(together_FA)
    total_NW <- length(together_NW)
    propFA[i] <- total_FA/length(U_FA_cor)
    propNW[i] <- total_NW/length(U_NW_cor)
  }
  P_net = sum(propNW == propData)/n
  P_FA = sum(propFA == propData)/n
  LR = (sum(propFA == propData)/n) / (sum(propNW == propData)/n )
  if(P_net>P_FA) pick = "network model" else pick = "unidimensional model"
  out = list(
    "Proportion of anomalous correlations" = propData, 
    "Probability of this many anomalous correlations given network model" = P_net, 
    "Probability of this many anomalous correlations given factor model" = P_FA, 
    "Likelihood ratio, network model in denominator" = LR,
    "Favoured model" = pick)
  return(out)
}

