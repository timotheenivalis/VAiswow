#### For variances back-tranformation ####

## function for expectation of a zero-inflated Poisson model
Ey<-function(l1,l2){
  plogis(-l2)*exp(l1)
}#end Ey()

## population mean given MC samples of latent variables
barW<-function(MC){
  l1<-MC[,1]
  l2<-MC[,2]
  return(mean(Ey(l1,l2)))
}#end barW()

# function to give the genetic variance of relative
# fitness given ZIP model M for MCMC sample s.  This
# is itself done by MC integration, with nMC samples
Va_w_zip<-function(m, s=1, nMC=50000, h=0.02,
                   fitnessfile=NULL, fixedpred=c("inbreeding", "Qgg"),
                   ranpred = c("id", "dam", "cohort"), 
                   focalvar = "id"
)
{
  
  muP <- m$Sol[s,"traitLBS"]
  muZI <- m$Sol[s,"traitzi_LBS"]
  if(!is.null(fixedpred))
  {
    if(is.null(fitnessfile)){
      stop("No fitnessfile provided but need to integrate over fixed effects")
    }else{ fit <- read.csv(fitnessfile)
    }
    for(i in 1:length(fixedpred))
    {
      muP <- muP + fit[,fixedpred[i]] * m$Sol[s,paste0("traitLBS:",fixedpred[i])]
      muZI <- muZI + fit[,fixedpred[i]] * m$Sol[s,paste0("traitzi_LBS:",fixedpred[i])]
    }
  }
  
  mu<-cbind(muP, muZI)
  
  PVP <- m$VCV[s,"traitLBS.units"] + 
    var(mu[,1])
  for (i in 1:length(ranpred))
  {
    PVP <- PVP + m$VCV[s,paste0("traitLBS:traitLBS.", ranpred[i])]
  }
  
  PVZI <- m$VCV[s,"traitzi_LBS.units"] + 
    var(mu[,2])
  for (i in 1:length(ranpred))
  {
    PVZI <- PVZI + m$VCV[s,paste0("traitzi_LBS:traitzi_LBS.", ranpred[i])]
  }
  
  PcovZIP <- cov(mu[,1], mu[,2])
  for (i in 1:length(ranpred))
  {
    PcovZIP <- PcovZIP + m$VCV[s,paste0("traitzi_LBS:traitLBS.", ranpred[i])]
  }
  
  
  sigma <-matrix(c(PVP,
                   PcovZIP,PcovZIP,
                   PVZI),2,2)
  
  MC_latent_smaples<-rmvnorm(nMC,colMeans(mu),sigma)
  
  # mean fitness
  muW<-barW(MC_latent_smaples)
  
  # values of mean fitness with perturbations of mean latent
  # values for use in getting derivatives by finite differences
  
  muW1<-barW(MC_latent_smaples+cbind(rep(h,nMC),0))
  muW2<-barW(MC_latent_smaples+cbind(0,rep(h,nMC)))
  
  # derivatives of mean relative fitness w.r.t. mean latent values
  d1<-(muW1/muW-1)/h
  d2<-(muW2/muW-1)/h
  
  # genetic variance in relative fitness from the linear 
  # projection of the distribution of latent values for the 
  #two components of the zip model
  if(focalvar=="units")
  {
    sigma_a<-matrix(data = c(m$VCV[s,paste0("traitLBS.",focalvar)],
                             0,0,
                             m$VCV[s,paste0("traitzi_LBS.",focalvar)]),
                    2,2)
    
  }else{
    sigma_a<-matrix(m$VCV[s,c(paste0(c("traitLBS:traitLBS.", 
                                       "traitzi_LBS:traitLBS.",
                                       "traitLBS:traitzi_LBS.",
                                       "traitzi_LBS:traitzi_LBS."), focalvar))],2,2)
  }
  d<-matrix(c(d1,d2),2,1)
  
  return(as.numeric(t(d)%*%sigma_a%*%d))
} #end Va_w_zip()

Va_w_zip_posterior <- function(m, nMC=50000, h=0.02,
                               fitnessfile, fixedpred=c("inbreeding", "Qgg"),
                               ranpred = c("id", "dam", "cohort"))
{
  length_posterior <- length(m$Sol[,1])
  ranpredunits <- c(ranpred, "units")
  # create a MCMC matrix for the transformed ZIP values
  btzip <- as.mcmc(matrix(data = NA, nrow = length_posterior,
                          ncol = length(ranpredunits)+1, 
                          dimnames = list(1:length_posterior, c(ranpredunits, "residual") ) ) )
  #copy MCMC attributes to the new matrix
  attr(btzip, "mcpar") <- attr(m$VCV, "mcpar")
  
  #compute the back-transformed effects for each variance component
  for(va in 1:length(ranpredunits))
  {
    btva <- vector(length = length_posterior)
    for(i in 1:length_posterior)
    {
      btva[i] <- Va_w_zip(m = m, s = i, fitnessfile = fitnessfile, 
                          fixedpred=fixedpred, nMC = nMC, h = h,
                          ranpred = ranpred, focalvar = ranpredunits[va])
    }
    btzip[, ranpredunits[va]] <- btva
  }
  fit <- read.csv(fitnessfile)
  btzip[,"residual"] <- var(fit$LBS) - rowSums(btzip[,-ncol(btzip)])
  m$btzip <- btzip
  return(m)
}#end Va_w_zip_posterior()


#### For model prediction ####
zip_pred<-function(m, s=1, nMC=50,
                   fitnessfile=NULL, fixedpred=c("inbreeding", "Qgg"),
                   ranpred = c("id", "dam", "cohort"))
{
  
  muP <- m$Sol[s,"traitLBS"]
  muZI <- m$Sol[s,"traitzi_LBS"]
  if(!is.null(fixedpred))
  {
    if(is.null(fitnessfile)){
      stop("No fitnessfile provided but need to integrate over fixed effects")
    }else{ fit <- read.csv(fitnessfile)
    }
    for(i in 1:length(fixedpred))
    {
      muP <- muP + fit[,fixedpred[i]] * m$Sol[s,paste0("traitLBS:",fixedpred[i])]
      muZI <- muZI + fit[,fixedpred[i]] * m$Sol[s,paste0("traitzi_LBS:",fixedpred[i])]
    }
  }
  mu<-cbind(muP, muZI)
  
  PVP <- m$VCV[s,"traitLBS.units"] + 
    var(mu[,1])
  for (i in 1:length(ranpred))
  {
    PVP <- PVP + m$VCV[s,paste0("traitLBS:traitLBS.", ranpred[i])]
  }
  
  PVZI <- m$VCV[s,"traitzi_LBS.units"] + 
    var(mu[,2])
  for (i in 1:length(ranpred))
  {
    PVZI <- PVZI + m$VCV[s,paste0("traitzi_LBS:traitzi_LBS.", ranpred[i])]
  }
  
  PcovZIP <- cov(mu[,1], mu[,2])
  for (i in 1:length(ranpred))
  {
    PcovZIP <- PcovZIP + m$VCV[s,paste0("traitzi_LBS:traitLBS.", ranpred[i])]
  }
  
  sigma <-matrix(c(PVP,
                   PcovZIP,PcovZIP,
                   PVZI),2,2)
  pred <- vector(length = nrow(mu)*nMC)
  for (i in 1:nrow(mu))
  {
    MC_latent_smaples<-rmvnorm(nMC,mu[i,],sigma)
    pred0 <- sapply(MC_latent_smaples[,2], function(x) {rbinom(n = 1, size = 1, prob = 1/(1+exp(x)))})
    predcond <- sapply(MC_latent_smaples[,1], function(x) {rpois(n = 1,lambda = exp(x))})
    pred[(1+(i-1)*(nMC)):(nMC*i)] <- pred0*predcond
  }
  return(pred)
}#end zip_pred()


zip_pred_posterior <- function(m, nMC=10,
                               fitnessfile=NULL, fixedpred=c("inbreeding", "Qgg"),
                               ranpred = c("id", "dam", "cohort"))
{
  length_posterior <- length(m$Sol[,1])
  fit <- read.csv(fitnessfile)
  limdis <- max(fit$LBS, na.rm = TRUE)+2
  predmat <- matrix(0, nrow = length_posterior, ncol = limdis+1)
  for (i in 1:length_posterior)
  {
    pred <- zip_pred(m = m, s = i, nMC = nMC, fitnessfile = fitnessfile,
                     fixedpred = fixedpred, ranpred = ranpred)
    pred[pred>limdis] <- limdis
    tabp <- table(pred)
    tabp <- tabp/sum(tabp)
    tabpnames <- as.numeric(names(tabp))
    for (j in 1:length(tabpnames))
    {
      predmat[i,tabpnames[j]+1] <- tabp[j]
    }
  }
  return(predmat)
}#end zip_pred_posterior()
