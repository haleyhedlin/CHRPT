#######
## code to obtain predictions from model estimates for new cohorts
## April 17, 2018
## Haley Hedlin
#######


#### read in required packages
#### run the install.packages lines if you have not already installed these packages
# install.packages("mice")
library(mice)
# install.packages("cmprsk")
library(cmprsk)


## contact Haley Hedlin at hedlin@stanford.edu to receive the model fit file
# load("WHIfit.Rdata")

#### read in data - the .csv should be saved in the same directory as this R code file
# assume that aspirin and statins are 0, 1
dat <- read.csv("val.csv")

# looking at a 10 year event horizon
hor <- 10
hordy <- hor*365.25

# define outcomes for the competing risk model
tmp.out.all <- cbind(dat$MI, dat$strk, dat$lc, dat$bc, dat$crc, 
                     dat$hip, dat$death)
tmp.time.all <- cbind(dat$MItime, dat$strktime, dat$lctime, dat$bctime, 
                      dat$crctime, dat$hiptime, dat$deathtime)
dat$tmall <- apply(tmp.time.all,1,min,na.rm=TRUE)
dat$statall <- rep(NA,nrow(dat))  
#0 = censored prior to any event, 1 = MI, 2 = stroke, 3 = LC, 4 = BC, 5 = CRC, 6 = hip fx, 7 = death
dat$statall[rowSums(tmp.out.all)==0] <- 0 # cens
dat$statall[dat$tmall >= hordy] <- 0
# loop through each woman who has an NA for statall
notcens <- which(is.na(dat$statall))
for(j in notcens){
  tmpstat <- which(tmp.time.all[j,]==min(tmp.time.all[j,]))
  if(length(tmpstat)==1){ 
    dat$statall[j] <- tmpstat
  } else if(sum(tmp.out.all[j,])==1){  
    dat$statall[j] <- which(tmp.out.all[j,]==1)
  } else if(sum(tmp.out.all[j,])==2 & tmp.out.all[j,7]==1){ 
    dat$statall[j] <- min(which(tmp.out.all[j,]==1))
  } 
}

## define vectors to use for individual outcomes
## for women who do not fall into one of the categories defined above
## tiebreakers needed for these women, to be used in the appropriate models
dat$statMI <- dat$statstrk <- dat$statlc <- dat$statbc <- dat$statcrc <- 
  dat$stathip <- dat$statdeath <- dat$statall

statNA <- which(is.na(dat$statall))
# loop through the outcomes
for(k in statNA){
  tmpstat <- which(tmp.time.all[k,]==min(tmp.time.all[k,]))
  if(length(tmpstat)==7) tmpstat <- tmpstat[tmp.out.all[k,]==1] 
  dat$statMI[k] <- ifelse(1 %in% tmpstat, 1, sample(tmpstat)[1])
  dat$statstrk[k] <- ifelse(2 %in% tmpstat, 2, sample(tmpstat)[1])
  dat$statlc[k] <- ifelse(3 %in% tmpstat, 3, sample(tmpstat)[1])
  dat$statbc[k] <- ifelse(4 %in% tmpstat, 4, sample(tmpstat)[1])
  dat$statcrc[k] <- ifelse(5 %in% tmpstat, 5, sample(tmpstat)[1])
  dat$stathip[k] <- ifelse(6 %in% tmpstat, 6, sample(tmpstat)[1])
  dat$statdeath[k] <- ifelse(7 %in% tmpstat, 7, sample(tmpstat)[1])
}

dat$statMI1 <- (dat$statMI==1)
dat$statstrk1 <- (dat$statstrk==2)
dat$statlc1 <- (dat$statlc==3)
dat$statbc1 <- (dat$statbc==4)
dat$statcrc1 <- (dat$statcrc==5)
dat$stathip1 <- (dat$stathip==6)
dat$statdeath1 <- (dat$statdeath==7)

rm(tmp.out.all, tmp.time.all)

dat$statMI1 <- (dat$statMI==1)
dat$statstrk1 <- (dat$statstrk==2)
dat$statlc1 <- (dat$statlc==3)
dat$statbc1 <- (dat$statbc==4)
dat$statcrc1 <- (dat$statcrc==5)
dat$stathip1 <- (dat$stathip==6)
dat$statdeath1 <- (dat$statdeath==7)

# calculate nelson aalen estimates here
dat$naimpMI <- nelsonaalen(dat,MItime,statMI1)
dat$naimpstrk <- nelsonaalen(dat,strktime,statstrk1)
dat$naimplc <- nelsonaalen(dat,lctime,statlc1)
dat$naimpbc <- nelsonaalen(dat,bctime,statbc1)
dat$naimpcrc <- nelsonaalen(dat,crctime,statcrc1)
dat$naimphip <- nelsonaalen(dat,hiptime,stathip1)
dat$naimpdeath <- nelsonaalen(dat,deathtime,statdeath1)


predmat <- matrix(1,nrow=length(c(1:89,91:93)),ncol=length(c(1:89,91:93)))
predmat[,90] <- 0  # don't use time var because we are using the nelson aalen estimate
diag(predmat) <- 0
impMI <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("MItime", "statMI", "naimpMI")))], 
               m=1, method='pmm',predictorMatrix=predmat)
impMI <- cbind(complete(impMI, action=1)) 
impstrk <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("strktime", "statstrk", "naimpstrk")))], 
               m=1, method='pmm',predictorMatrix=predmat)
impstrk <- cbind(complete(impstrk, action=1)) 
implc <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("lctime", "statlc", "naimplc")))], 
               m=1, method='pmm',predictorMatrix=predmat)
implc <- cbind(complete(implc, action=1)) 
impbc <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("bctime", "statbc", "naimpbc")))], 
               m=1, method='pmm',predictorMatrix=predmat)
impbc <- cbind(complete(impbc, action=1)) 
impcrc <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("crctime", "statcrc", "naimpcrc")))], 
               m=1, method='pmm',predictorMatrix=predmat)
impcrc <- cbind(complete(impcrc, action=1)) 
imphip <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                        c("hiptime", "stathip", "naimphip")))], 
               m=1, method='pmm',predictorMatrix=predmat)
imphip <- cbind(complete(imphip, action=1)) 
impdeath <- mice(data=dat[,c(1:89,which(names(dat) %in% 
                                          c("deathtime", "statdeath", 
                                            "naimpdeath")))], 
                 m=1, method='pmm',predictorMatrix=predmat)
impdeath <- cbind(complete(impdeath, action=1)) 

# define interactions
impMI$diabXpkyrs <- impMI$DIAB*impMI$PACKYRS
impMI$diabXhtntrt2 <- impMI$DIAB*impMI$HTNTRT_2
imphip$fh552Xalc6 <- imphip$fracthip55_2*imphip$ALCOHOL_6
impdeath$diabXpkyrs <- impdeath$DIAB*impdeath$PACKYRS
impdeath$diabXhtntrt2 <- impdeath$DIAB*impdeath$HTNTRT_2
impdeath$fh552Xalc6 <- impdeath$fracthip55_2*impdeath$ALCOHOL_6



# here is the code for predictions --
ndat <- nrow(dat)
MIval10 <- strkval10 <- lcval10 <- bcval10 <- crcval10 <- hipval10 <- 
  deathval10 <- numeric(ndat)
indMI <- findInterval(hordy,MIfit$uftime)
indstrk <- findInterval(hordy,strkfit$uftime)
indlc <- findInterval(hordy,lcfit$uftime)
indbc <- findInterval(hordy,bcfit$uftime)
indcrc <- findInterval(hordy,crcfit$uftime)
indhip <- findInterval(hordy,hipfit$uftime)
inddeath <- findInterval(hordy,deathfit$uftime)
for(i in 1:ndat){
  MIval10[i] <- predict(MIfit, cov1=impMI[i,c(MIvars,93:94)])[indMI,-1]
  strkval10[i] <- predict(strkfit, cov1=impstrk[i,strkvars])[indstrk,-1]
  lcval10[i] <- predict(lcfit, cov1=implc[i,lcvars])[indlc,-1]
  bcval10[i] <- predict(bcfit, cov1=impbc[i,bcvars])[indbc,-1]
  crcval10[i] <- predict(crcfit, cov1=impcrc[i,crcvars])[indcrc,-1]
  hipval10[i] <- predict(hipfit, cov1=imphip[i,c(hipvars,93)])[indhip,-1]
  deathval10[i] <- predict(deathfit, cov1=impdeath[i,c(deathvars,93:95)])[
    inddeath,-1]
  }

##### output the individual women's predictions in a .csv
predvalmat <- cbind(MIval10,strkval10,lcval10,bcval10,crcval10,hipval10,
                    deathval10)
write.csv(predvalmat,file="ValPred.csv")

covmatvec <- expression(impMI,impstrk,implc,impbc,impcrc,imphip,impdeath)


decl5 <- apply(predvalmat,2,quantile,seq(0,1,by=0.2))
qtlcat <- matrix(NA, nrow=ndat, ncol=ncol(decl5))
KMmat <- matrix(NA, nrow=ndat, ncol=ncol(decl5)*3)  
# the first column contains survival estimate
# the second column contains the time
# the third column indicates stratum
KMlist <- NULL
for(i in 1:ncol(decl5)){
  covmat <- eval(covmatvec[i])
  qtlcat[,i] <- cut(predvalmat[,i],decl5[,i],include.lowest=TRUE)
  ## need to get this so I can tell the quintiles
  tmp <- survfit(Surv(covmat[,90], covmat[,91]==i)~qtlcat[,i])
  ## strata contains the n for each stratum
  KMmat[1:length(tmp$surv),(3*i-2)] <- tmp$surv
  KMmat[1:length(tmp$surv),(3*i-1)] <- tmp$time
  tmpind = cumsum(tmp$strata)
  KMmat[1:tmpind[1],(3*i)] <- 1
  KMmat[(tmpind[1]+1):tmpind[2],(3*i)] <- 2
  KMmat[(tmpind[2]+1):tmpind[3],(3*i)] <- 3
  KMmat[(tmpind[3]+1):tmpind[4],(3*i)] <- 4
  KMmat[(tmpind[4]+1):tmpind[5],(3*i)] <- 5
  
  ## also output a list of the KM objects in case I want to plot CI or censoring
  KMlist[[i]] = tmp
}

write.csv(KMmat, file="KMmatval.csv")
save(KMlist, file="KMlistval.Rdata")



## methods of validation

# plot the histogram of the distribution of predicted probabilities

# # create the plot of predictions based on training data
# # then add in the predictions from the test data
# MIdecileT <- quantile(MItest10,seq(0,1,by=0.1))
# MIdeccatT <- cut(MItest10,MIdecileT,include.lowest=TRUE)
# MIdecrateT <- aggregate(test.impMI$statMI==1, by=list(MIdeccatT), FUN=mean)
# points(seq(0.05, 0.95, by=0.1),MIdecrateT[,2], col="green", pch=19)

# look at KM plot separating deciles (or quintiles)




