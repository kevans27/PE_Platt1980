###This script will output a one row 'data frame' called modelFitParameters which will contain the model output, plus Ik and Im which are derived from the output.
###It will also plot your data to show how the curve performed
###For a single sample, depending on your computer, I'd expect this to take somewhere between five minutes and an hour. If it goes much beyond that, something has probably gone wrong. Check your data format first, but I'm happy to help check it out too.

library(rstan)

##Putting in fake data. Replace with your own
##Light of each incubation
light <- c(533, 1446, 21.1, 424, 681, 86.3, 132, 240, 2665, 16.3, 163, 731, 61, 
           900, 1520, 108, 49.4, 363, 12.5, 150, 504, 52.8, 217, 32.1)
##Productivity response rate (in your units of interest)
phot <- c(1.0852575, 0.9827610, 0.5004243, 1.2601046, 1.1696665, 0.9707026, 
          1.0370239, 1.2179001, 0.6692421, 0.6390961, 1.2721630, 1.1455496, 
          0.9104105, 1.0490823, 0.8018847, 1.1998125, 0.7837971, 1.1395204, 
          0.4099862, 1.1998125, 1.1636372, 0.9284981, 1.1636372, 0.7898263)

##These are some extractable parameters that don't emerge directly from the model fits, but emerge from the model output.
IMFun <- function(pbsVar, alphVar, betVar){
  x = pbsVar/alphVar * log((alphVar + betVar)/betVar)
  return(x)
}
PBMFun <- function(pbsVar, alphVar, betVar){
  x = pbsVar*(alphVar/(alphVar+betVar))*(betVar/(alphVar+betVar))**(betVar/alphVar)
  return(x)
}
IKFun <- function(pbmVar, alphVar){
  x = pbmVar/alphVar
  return(x)
}



sink("nonHierarchicalPlatt1980.stan")
cat("
      data {
      int <lower=1> K; //number of samples 
      vector [K] Light; //input model vector
      vector [K] Phot; //response model vector
      }
      
      parameters {
      real <lower=0> alph;
      real <lower=0> bet;
      real <lower=0> PBS;
      real <lower=0> sigm;
      }
      
      transformed parameters {
      real PBM;
        PBM = PBS*(alph/(alph+bet))*(bet/(alph+bet))^(bet/alph);
      }
      
      model { 
      //priors. 
      alph ~ normal(0.03, 0.02);
      bet ~ normal(0.0003, 0.0003);
      PBS ~ normal(1.5, 1);

      //likelihood    	
      for (j in 1:K) {
        Phot[j] ~ normal(PBS*(1-exp(-(alph*Light[j])/PBS))*(exp(-(bet*Light[j])/PBS)), sigm);
      }
      }",
    fill=TRUE)
sink()

##Now to stan it
stanDat <- list(K=length(light), Light=light, Phot=phot)

stanOutput <- stan(file="nonHierarchicalPlatt1980.stan", data = stanDat, 
                  iter = 10000, chains = 8, control = list(adapt_delta = 0.99, 
                                                           max_treedepth = 15))


##This is me extracting the output. There's probably a cleaner way to do this
allVals <- extract(stanOutput, permuted = TRUE)
alphaVal <- mean(allVals$alph)
alphaValSd <- sd(allVals$alph)
betaVal <- mean(allVals$bet)
betaValSd <- sd(allVals$bet)
PBSVal <- mean(allVals$PBS)
PBSValSd <- sd(allVals$PBS)
tronSigm <- mean(allVals$sigm)

IMAll <- IMFun(allVals$PBS, allVals$alph, allVals$bet)
IMVal <- mean(IMAll)
IMSd <- sd(IMAll)
PBMAll <- PBMFun(allVals$PBS, allVals$alph, allVals$bet)
PBMVal <- mean(PBMAll)
PBMSd <- sd(PBMAll)
IKAll <- IKFun(PBMAll, allVals$alph)
IKVal <- mean(IKAll)
IKSd <- sd(IKAll)


#Turns the model fit parameters into a 'dataframe'
modelFitParameters <- data.frame("Alph" = alphaVal,
                                 "Alph_Sd" = alphaValSd,
                                 "Bet" = betaVal, 
                                 "Bet_Sd" = betaValSd,
                                 "PBS" = PBSVal, 
                                 "PBS_Sd" = PBSValSd,
                                 "PBM" = PBMVal,
                                 "PBM_Sd" = PBMSd,
                                 "IM" = IMVal,
                                 "IM_Sd" = IMSd,
                                 "IK" = IKVal,
                                 "IK_Sd" = IKSd)


#Plots the sample fit
plot(light, phot, pch = 18, xlab = "Light", ylab = "Photosynthetic production")
forCurve <- seq(0, max(light), length.out = 1000)
curveOutput <- modelFitParameters$PBS*
  (1-exp(-(modelFitParameters$Alph*forCurve)/modelFitParameters$PBS))*
  exp(-(modelFitParameters$Bet*forCurve)/modelFitParameters$PBS)
lines(forCurve, curveOutput, col = 'red')
