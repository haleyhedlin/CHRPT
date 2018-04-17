
#######
## model building and internal validation code for event-first model
## April 17, 2018
## Haley Hedlin
#######


library(Hmisc)
library(riskRegression)
library(cmprsk)
library(pec)
library(CoxBoost)
library(mice)
library(crrstep)
library(ggplot2)


# factor2ind function from http://www.nature.com/bmt/journal/v45/n9/full/bmt2009359a.html#app1
# x is the variable to be converted to a factor
# baseline is the reference categoryg
factor2ind <- function(x,baseline){  
  baseline <- as.character(baseline)
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  ## updated on 5/20 to allow for missing values to remain missing
  ## added the following two lines
  x[x=="NaN"] <- NA
  x <- factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep= "_"))  
  return(X[,-1,drop=FALSE])
}


## read in your data file
## use the same naming conventions as in val.csv
#dat <- 


test = 1   # 1 for test, 0 for training
hor <- 10   # number of years for prediction - 5, 10, or 15 years
hordy <- hor*365.25

 
# define outcome day as minimum of event day and the end of ext1 or horizon
dat$MItime <- pmin(dat$MIDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$strktime <- pmin(dat$STROKEDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$lctime <- pmin(dat$LUNGDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$bctime <- pmin(dat$BREASTDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$crctime <- pmin(dat$COLORECTALDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$hiptime <- pmin(dat$BKHIPDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)
dat$deathtime <- pmin(dat$DEATHDY, dat$ENDEXT1DY, hordy, na.rm=TRUE)

tmp.out.all <- cbind(dat$MI, dat$STROKE, dat$LUNG, dat$BREAST, dat$COLORECTAL, 
                     dat$BKHIP.y, dat$DEATH)
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

dat$statMI <- dat$statstrk <- dat$statlc <- dat$statbc <- dat$statcrc <- 
  dat$stathip <- dat$statdeath <- dat$statall

# subdistribution hazard models (eg Fine and Gray) have separate models for each outcome of interest 
# (issues with this discussed in the Beyersmann ebook on competing risks)
# I think this means it's ok to define ties differently in each outcome

statNA <- which(is.na(dat$statall))
# loop through the eight outcomes
# this code defines the outcomes as event of interest, competing event, or censored
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

## covariate coding specific to WHI data

# note - 0 and -1 defined the same for agefbir2 and parity 
# I'm redefining agefbir2 to collapse those levels
dat$agefbir3 <- dat$agefbir2
dat$agefbir3[dat$agefbir3==-1] <- 0

# collapsing the levels in MENARCHE and qsmokage
dat$menarche <- dat$MENARCHE
dat$menarche[dat$menarche==1] <- 2
dat$menarche[dat$menarche > 6] <- 7
dat$qsmokage2 <- dat$qsmokage # making the categories intervals of 10 instead of 5
dat$qsmokage2[dat$qsmokage2==1] <- 4  # further combined qsmokeage2 to avoid very small numbers
dat$qsmokage2[dat$qsmokage2==2] <- 4 
dat$qsmokage2[dat$qsmokage2==3] <- 4
dat$qsmokage2[dat$qsmokage2==5] <- 6
dat$qsmokage2[dat$qsmokage2==7] <- 8
dat$qsmokage2[dat$qsmokage2==9] <- 10
dat$qsmokage2[dat$qsmokage2==11] <- 10  # further combined qsmokeage2 to avoid very small numbers

# combined ooph so that part of an ovary and an unknown number of ovaries combined with one ovary removed
dat$ooph <- dat$OOPH
dat$ooph[dat$ooph>2] <- 1

# mimomage, midadage combining so that "yes, don't know age" grouped with "yes, age 65 or older"
dat$mimomage2 <- dat$mimomage
dat$mimomage2[dat$mimomage2 > 3] <- 3
dat$midadage2 <- dat$midadage 
dat$midadage2[dat$midadage2 > 3] <- 3

# collapsing bkhip55 and fract55 into a single variable
# 0 = never broke bone, 1 & 2 never broke hip, 3 & 4 broke hip
# 1 = broke bone <= age 55 or not sure age, 2 = first broke bone at > age 55
# 3 = broke hip 
# combined broke hip categories due to small numbers
dat$fracthip55 <- dat$fract55
dat$fracthip55[!is.na(dat$bkhip55) & dat$bkhip55!=0] <- dat$bkhip55[
  !is.na(dat$bkhip55) & dat$bkhip55!=0] + 2
dat$fracthip55[dat$fracthip55==4] <- 3

# create a broke hip or not for mom and dad
dat$BKHIPMOM[dat$BKBONMOM==0] <- 0
dat$bkhipmom <- dat$BKHIPMOM
dat$bkhipmom[dat$bkhipmom == 9] <- NA
dat$bkhipmom[dat$bkhipmom >0] <- 1
dat$BKHIPDAD[dat$BKBONDAD==0] <- 0
dat$bkhipdad <- dat$BKHIPDAD
dat$bkhipdad[dat$bkhipdad == 9] <- NA
dat$bkhipdad[dat$bkhipdad >0] <- 1


# combine ethnic categories
# keep NH black, white, make AI/AN, API, Hispanic, and Other into a broader "Other"
dat$ethnic2 <- dat$ETHNIC
dat$ethnic2[dat$ethnic2 < 3 | dat$ethnic==4] <- 8

# combine GENHEL
dat$genhel <- dat$GENHEL
dat$genhel[dat$genhel==5] <- 4

# collapse strkreln
dat$strkreln2 <- dat$strkreln
dat$strkreln2[dat$strkreln2 > 2] <- 2

# collapse oldest ages in momdiedage and daddiedage
dat$momdiedage2 <- dat$momdiedage
dat$momdiedage2[dat$momdiedage2 > 7] <- 7
dat$daddiedage2 <- dat$daddiedage
dat$daddiedage2[dat$daddiedage2 > 7] <- 7

# create nelson-aalen estimates for imputing time to event variables
# need one for each outcome

dat$statMI1 <- (dat$statMI==1)
dat$statstrk1 <- (dat$statstrk==2)
dat$statlc1 <- (dat$statlc==3)
dat$statbc1 <- (dat$statbc==4)
dat$statcrc1 <- (dat$statcrc==5)
dat$stathip1 <- (dat$stathip==6)
dat$statdeath1 <- (dat$statdeath==7)

dat$naimpMI <- nelsonaalen(dat,MItime,statMI1)
dat$naimpstrk <- nelsonaalen(dat,strktime,statstrk1)
dat$naimplc <- nelsonaalen(dat,lctime,statlc1)
dat$naimpbc <- nelsonaalen(dat,bctime,statbc1)
dat$naimpcrc <- nelsonaalen(dat,crctime,statcrc1)
dat$naimphip <- nelsonaalen(dat,hiptime,stathip1)
dat$naimpdeath <- nelsonaalen(dat,deathtime,statdeath1)
  

dat.fac <- cbind(AGE=dat$AGE,factor2ind(dat$ethnic2,5), DIAB=dat$DIAB, # white
                 HICHOLRP=dat$HICHOLRP, MIGRAINE=dat$MIGRAINE, 
                 ATRIALFB=dat$ATRIALFB, undthy=dat$undthy, 
                 ovrthy=dat$ovrthy, factor2ind(dat$fracthip55,0), # never broke bone
                 factor2ind(dat$HTNTRT,0), factor2ind(dat$menarche,5), # no HTN
                 BRSTFEED=dat$BRSTFEED, factor2ind(dat$ooph,0), # none removed
                 BRSTBIOP=dat$BRSTBIOP, factor2ind(dat$PARITY,-1), # never preg
                 factor2ind(dat$agefbir3,0), factor2ind(dat$momdiedage2,0), # never preg
                 factor2ind(dat$daddiedage2,0), MIREL=dat$MIREL, BKBONREL=dat$BKBONREL,
                 bkhipmom=dat$bkhipmom, bkhipdad=dat$bkhipdad,
                 brcafrel2=dat$brcafrel2, colofrel2=dat$colofrel2,
                 colomrel2=dat$colomrel2, factor2ind(dat$mimomage2,0), # no MI
                 factor2ind(dat$strkreln2,0),
                 CANCFREL=dat$CANCFREL, CANCMREL=dat$CANCMREL,
                 factor2ind(dat$midadage2,0), LACTDIET=dat$LACTDIET, # no MI,
                 TEPIWK=dat$TEPIWK, TMINWK=dat$TMINWK, 
                 factor2ind(dat$ALCOHOL,1), PACKYRS=dat$PACKYRS, 
                 factor2ind(dat$qsmokage2,0), PULSE30=dat$PULSE30, # non-drinker
                 HEIGHT=dat$HEIGHT, WEIGHT=dat$WEIGHT, WAIST=dat$WAIST, 
                 HIP=dat$HIP, BMI=dat$BMI, WHR=dat$WHR, SYST=dat$SYST,
                 genhel=factor2ind(dat$genhel,1), # excellent
                 aspirin=as.numeric(dat$aspirin),statin=as.numeric(dat$statin),
                 MItime=dat$MItime, statMI=dat$statMI1, naMI=dat$naimpMI,
                 strktime=dat$strktime, statstrk=dat$statstrk1, nastrk=dat$naimpstrk,
                 lctime=dat$lctime, statlc=dat$statlc1, nalc=dat$naimplc,
                 bctime=dat$bctime, statbc=dat$statbc1, nabc=dat$naimpbc,
                 crctime=dat$crctime, statcrc=dat$statcrc1, nacrc=dat$naimpcrc,
                 hiptime=dat$hiptime, stathip=dat$stathip1, nahip=dat$naimphip,
                 deathtime=dat$deathtime, statdeath=dat$statdeath1, nadeath=dat$naimpdeath,
                 REGION=dat$REGION, HRTARM=dat$HRTARM)
dat.fac <- as.data.frame(dat.fac)
names(dat.fac) <- gsub("dat\\$", "", names(dat.fac))

# check missingness
#colSums(is.na(dat.fac))

# fit the competing risk model 

MIvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                c("AGE","DIAB","HICHOLRP","ATRIALFB","TEPIWK",
                                  "TMINWK","PACKYRS","PULSE30","HEIGHT",
                                  "MIREL","WEIGHT","WAIST","HIP",# "BMI", "WHR",
                                  "SYST","aspirin")),
            apply(as.matrix(c("ethnic","HTNTRT","ooph","mimomage","midadage",
                              "ALCOHOL","qsmokage","genhel")), 1, grep, 
                  x=names(dat.fac), ignore.case=TRUE))))
strkvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                c("AGE","DIAB","HICHOLRP","MIGRAINE",
                                  "ATRIALFB","PACKYRS","PULSE30",
                                  "WEIGHT","HEIGHT","WAIST",#"BMI","WHR",
                                  "SYST", "aspirin")),
                        apply(as.matrix(c("ethnic","HTNTRT","qsmokage",
                                          "genhel","strkreln2")), 1, grep, 
                              x=names(dat.fac), ignore.case=TRUE))))
lcvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                c("AGE","DIAB","PACKYRS","WEIGHT","WAIST")),
                        apply(as.matrix(c("ethnic","qsmokage","genhel")), 1, grep, 
                              x=names(dat.fac), ignore.case=TRUE))))
bcvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                c("AGE","DIAB","BRSTFEED","BRSTBIOP","PACKYRS",
                                  "WEIGHT","WAIST")),
                        apply(as.matrix(c("ethnic","menarche","ooph","PARITY",
                                          "agefbir","brcafrel","ALCOHOL",
                                          "qsmokage","genhel")), 1, grep, 
                              x=names(dat.fac), ignore.case=TRUE))))
crcvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                c("AGE","DIAB","PACKYRS",
                                  "WEIGHT","WAIST","aspirin" )),
                        apply(as.matrix(c("ethnic","colofrel","colomrel",
                                          "qsmokage","genhel")), 1, grep, 
                              x=names(dat.fac), ignore.case=TRUE))))
hipvars <- sort(unlist(c(which(names(dat.fac) %in% 
                                 c("AGE","DIAB","ATRIALFB","undthy","ovrthy",
                                   "BRSTFEED","BKBONREL","LACTDIET","TEPIWK",
                                   "TMINWK","PACKYRS","HEIGHT","WEIGHT",
                                   "WAIST","HIP","bkhipmom",#"BMI","WHR",
                                   "bkhipdad")),
                         apply(as.matrix(c("ethnic","fracthip","ooph","ALCOHOL",
                                           "qsmokage","genhel")), 1, grep, 
                               x=names(dat.fac), ignore.case=TRUE))))
deathvars <- 1:89



#### set aside a test set to be used after validation

set.seed(161)
# randomly select three regions to be the training dataset
testREGION <- sample(1:4,1)

if(test==0){
  tdat <- dat.fac[dat.fac$REGION!=testREGION,!names(dat.fac)=="REGION"]  
}
if(test==1){
  tdat <- dat.fac[dat.fac$REGION==testREGION,!names(dat.fac)=="REGION"]
}


## need a separate imputation for each outcome type
predmat <- matrix(1,nrow=length(c(1:89,91:93)),ncol=length(c(1:89,91:93)))
predmat[,90] <- 0  # don't use time var because we are using the nelson aalen estimate
diag(predmat) <- 0
dat.fac.impMI <- mice(data=tdat[,c(1:89,91:93)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.impMI <- cbind(complete(dat.fac.impMI, action=1))
dat.fac.impstrk <- mice(data=tdat[,c(1:89,94:96)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.impstrk <- cbind(complete(dat.fac.impstrk, action=1))
dat.fac.implc <- mice(data=tdat[,c(1:89,97:99)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.implc <- cbind(complete(dat.fac.implc, action=1))
dat.fac.impbc <- mice(data=tdat[,c(1:89,100:102)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.impbc <- cbind(complete(dat.fac.impbc, action=1))
dat.fac.impcrc <- mice(data=tdat[,c(1:89,103:105)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.impcrc <- cbind(complete(dat.fac.impcrc, action=1))
dat.fac.imphip <- mice(data=tdat[,c(1:89,106:108)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.imphip <- cbind(complete(dat.fac.imphip, action=1))
dat.fac.impdeath <- mice(data=tdat[,c(1:89,109:111)], m=1, method='pmm',
                      predictorMatrix=predmat)
cov.fac.impdeath <- cbind(complete(dat.fac.impdeath, action=1))


# run all models once with a small subset 
subfit <- sample(1:nrow(tdat), size=50000)
subfit01 <- (1:nrow(tdat) %in% subfit)
#subfit2 <- sample(1:nrow(trdat)[-subfit], size=50000)
#subfit201 <- (1:nrow(trdat)[-subfit] %in% subfit2)

sink(paste0(pathname,"MIC",hor,".txt"))
MIform <- as.formula(paste0("MItime~",paste(names(cov.fac.impMI)[MIvars],
                                            collapse="+")))
MIformfull <- as.formula(paste0("MItime~",paste(names(cov.fac.impMI)[MIvars],
                                            collapse="+"),
  "+BMI+WHR+DIAB*PACKYRS+DIAB*HICHOLRP+DIAB*HTNTRT_1+DIAB*HTNTRT_2+
  HICHOLRP*PACKYRS+HICHOLRP*HTNTRT_1+HICHOLRP*HTNTRT_1+HTNTRT_1*PACKYRS+
  HTNTRT_2*PACKYRS+MIREL*DIAB+MIREL*PACKYRS+MIREL*HICHOLRP+MIREL*HTNTRT_1+
  MIREL*HTNTRT_2+qsmokage2_4*DIAB+qsmokage2_6*DIAB+qsmokage2_8*DIAB+
  qsmokage2_10*DIAB+qsmokage2_12*DIAB+qsmokage2_4*HICHOLRP+
  qsmokage2_6*HICHOLRP+qsmokage2_8*HICHOLRP+qsmokage2_10*HICHOLRP+
  qsmokage2_12*HICHOLRP+qsmokage2_4*HTNTRT_1+qsmokage2_4*HTNTRT_2+
  qsmokage2_6*HTNTRT_2+qsmokage2_6*HTNTRT_1+qsmokage2_8*HTNTRT_2+
  qsmokage2_8*HTNTRT_1+qsmokage2_10*HTNTRT_2+qsmokage2_10*HTNTRT_1+
  qsmokage2_12*HTNTRT_2+qsmokage2_12*HTNTRT_1+qsmokage2_4*MIREL+
  qsmokage2_6*MIREL+qsmokage2_8*MIREL+qsmokage2_10*MIREL+qsmokage2_12*MIREL+
  SYST*DIAB+SYST*HICHOLRP+SYST*PACKYRS+SYST*qsmokage2_4+SYST*qsmokage2_6+
  SYST*qsmokage2_8+SYST*qsmokage2_10+SYST*qsmokage2_12+aspirin*DIAB+
  aspirin*HICHOLRP+aspirin*PACKYRS+aspirin*qsmokage2_4+aspirin*qsmokage2_6
  +aspirin*qsmokage2_8+aspirin*qsmokage2_10+aspirin*qsmokage2_12"))
mics <- crrstep(MIformfull, scope.min=MIform, etype=statMI,
                failcode=1, data=cov.fac.impMI, direction="forward",
		subset=subfit01,
                criterion="BICcr")
write.csv(mics, file=paste0(pathname,"MIC",hor,".csv"))
sink()
unlink(paste0(pathname,"MIC",hor,".txt"))

sink(paste0(pathname,"strkC",hor,".txt"))
strkform <- as.formula(paste0("strktime~",paste(names(cov.fac.impstrk)[strkvars],
                                            collapse="+")))
strkformfull <- as.formula(paste0("strktime~",paste(names(cov.fac.impstrk
                                                          )[strkvars],
                                                collapse="+"),
  "+BMI+WHR+DIAB*PACKYRS+DIAB*HICHOLRP+DIAB*HTNTRT_1+DIAB*HTNTRT_2+
  HICHOLRP*PACKYRS+HICHOLRP*HTNTRT_1+HICHOLRP*HTNTRT_2+ATRIALFB*DIAB+
  ATRIALFB*HTNTRT_1+ATRIALFB*HTNTRT_2+ATRIALFB*PACKYRS+ATRIALFB*qsmokage2_4+
  ATRIALFB*qsmokage2_6+ATRIALFB*qsmokage2_8+ATRIALFB*qsmokage2_10+
  ATRIALFB*qsmokage2_12+HTNTRT_1*PACKYRS+HTNTRT_2*PACKYRS+strkreln2_1*DIAB+
  strkreln2_2*DIAB+strkreln2_1*PACKYRS+
  strkreln2_2*PACKYRS+qsmokage2_4*DIAB+qsmokage2_6*DIAB+qsmokage2_8*DIAB+
  qsmokage2_10*DIAB+qsmokage2_12*DIAB+qsmokage2_4*HICHOLRP+
  qsmokage2_6*HICHOLRP+qsmokage2_8*HICHOLRP+qsmokage2_10*HICHOLRP+
  qsmokage2_12*HICHOLRP+qsmokage2_4*HTNTRT_1+qsmokage2_4*HTNTRT_2+
  qsmokage2_6*HTNTRT_1+qsmokage2_6*HTNTRT_2+qsmokage2_8*HTNTRT_1+
  qsmokage2_8*HTNTRT_2+qsmokage2_10*HTNTRT_1+qsmokage2_10*HTNTRT_2+
  qsmokage2_12*HTNTRT_1+qsmokage2_12*HTNTRT_2+qsmokage2_4*strkreln2_1+
  qsmokage2_6*strkreln2_1+qsmokage2_8*strkreln2_1+qsmokage2_10*strkreln2_1+
  qsmokage2_12*strkreln2_1+qsmokage2_4*strkreln2_2+qsmokage2_6*strkreln2_2+
  qsmokage2_8*strkreln2_2+qsmokage2_10*strkreln2_2+qsmokage2_12*strkreln2_2+
  SYST*DIAB+SYST*PACKYRS+SYST*qsmokage2_4+SYST*qsmokage2_6+SYST*qsmokage2_8+
  SYST*qsmokage2_10+SYST*qsmokage2_12"))
strkcs <- crrstep(strkformfull, scope.min=strkform,
                  etype=statstrk,failcode=1, data=cov.fac.impstrk,
		  subset=subfit01,
                  direction="forward",criterion="BICcr")
write.csv(strkcs, file=paste0(pathname,"strkC",hor,".csv"))
sink()
unlink(paste0(pathname,"strkC",hor,".txt"))

sink(paste0(pathname,"bcC",hor,".txt"))
bcform <- as.formula(paste0("bctime~",paste(names(cov.fac.impbc)[bcvars],
                                            collapse="+")))
bcformfull <- as.formula(paste0("bctime~",paste(names(cov.fac.impbc)[bcvars],
                                            collapse="+"),
  "+brcafrel2*menarche_2+brcafrel2*menarche_3+brcafrel2*menarche_4+
  brcafrel2*menarche_6+brcafrel2*menarche_7+brcafrel2*BRSTFEED+
  brcafrel2*ooph_1+brcafrel2*ooph_2+brcafrel2*BRSTBIOP+brcafrel2*PARITY_1+
  brcafrel2*PARITY_2+brcafrel2*PARITY_3+brcafrel2*PARITY_4+brcafrel2*PARITY_5+
  brcafrel2*agefbir3_1+brcafrel2*agefbir3_2+brcafrel2*agefbir3_3+
  brcafrel2*ALCOHOL_2+brcafrel2*ALCOHOL_3+brcafrel2*ALCOHOL_4+
  brcafrel2*ALCOHOL_5+brcafrel2*ALCOHOL_6+brcafrel2*PACKYRS+
  brcafrel2*qsmokage2_4+brcafrel2*qsmokage2_6+brcafrel2*qsmokage2_8+
  brcafrel2*qsmokage2_10+brcafrel2*qsmokage2_12"))
bccs <- crrstep(bcformfull, scope.min=bcform, etype=statbc,
                failcode=1, data=cov.fac.impbc, direction="forward",
		subset=subfit01,
                criterion="BICcr")
write.csv(bccs, file=paste0(pathname,"bcC",hor,".csv"))
sink()
unlink(paste0(pathname,"bcC",hor,".txt"))

sink(paste0(pathname,"crcC",hor,".txt"))
crcform <- as.formula(paste0("crctime~",paste(names(cov.fac.impcrc)[crcvars],
                                            collapse="+")))
crcformfull <- as.formula(paste0("crctime~",
                                 paste(names(cov.fac.impcrc)[crcvars],
                                                  collapse="+"),
                                 "+AGE*colofrel2+AGE*colomrel2+
                                 PACKYRS*colofrel2+PACKYRS*colomrel2"))
crccs <- crrstep(crcformfull, scope.min=crcform, etype=statcrc,
                 failcode=1, data=cov.fac.impcrc, direction="forward",
		 subset=subfit01,
                 criterion="BICcr")
write.csv(crccs, file=paste0(pathname,"crcC",hor,".csv"))
sink()
unlink(paste0(pathname,"crcC",hor,".txt"))

sink(paste0(pathname,"hipC",hor,".txt"))
hipform <- as.formula(paste0("hiptime~",paste(names(cov.fac.imphip)[hipvars],
                                              collapse="+")))
hipformfull <- as.formula(paste0("hiptime~",paste(names(cov.fac.imphip
                                                        )[hipvars],
                                              collapse="+"),
  "+BMI+WHR+BMI*WHR+BKBONREL*LACTDIET+fracthip55_1*LACTDIET+fracthip55_2*LACTDIET+
  fracthip55_3*LACTDIET+BKBONREL*TEPIWK+fracthip55_1*TEPIWK+
  fracthip55_2*TEPIWK+BKBONREL*TMINWK+fracthip55_1*TMINWK+fracthip55_2*TMINWK+
  BKBONREL*ALCOHOL_2+BKBONREL*ALCOHOL_3+BKBONREL*ALCOHOL_4+BKBONREL*ALCOHOL_5+
  BKBONREL*ALCOHOL_6+fracthip55_1*ALCOHOL_2+fracthip55_2*ALCOHOL_2+
  fracthip55_1*ALCOHOL_3+fracthip55_2*ALCOHOL_3+fracthip55_1*ALCOHOL_4+
  fracthip55_2*ALCOHOL_4+fracthip55_1*ALCOHOL_5+fracthip55_2*ALCOHOL_5+
  fracthip55_1*ALCOHOL_6+fracthip55_2*ALCOHOL_6+BKBONREL*PACKYRS+
  fracthip55_1*PACKYRS+fracthip55_2*PACKYRS+BKBONREL*qsmokage2_4+
  BKBONREL*qsmokage2_6+BKBONREL*qsmokage2_8+BKBONREL*qsmokage2_10+
  BKBONREL*qsmokage2_12+fracthip55_1*qsmokage2_4+fracthip55_2*qsmokage2_4+
  fracthip55_1*qsmokage2_6+fracthip55_2*qsmokage2_6+fracthip55_1*qsmokage2_8+
  fracthip55_2*qsmokage2_8+fracthip55_1*qsmokage2_10+fracthip55_2*qsmokage2_10+
  fracthip55_1*qsmokage2_12+fracthip55_2*qsmokage2_12"))
hipcs <- crrstep(hipformfull, scope.min=hipform, etype=stathip,
                 failcode=1, data=cov.fac.imphip, direction="forward",
		  subset=subfit01,
                 criterion="BICcr")
write.csv(hipcs, file=paste0(pathname,"hipC",hor,".csv"))
sink()
unlink(paste0(pathname,"hipC",hor,".txt"))


# # fit the models with all data in the training data
## use the parameters chosen in crrstep

## interaction for 15 year model
if(hor == 15){ 
  cov.fac.impMI$diabXhtntrt2 <- cov.fac.impMI$DIAB*cov.fac.impMI$HTNTRT_2
  MIfit <- crr(cov.fac.impMI$MItime, cov.fac.impMI$statMI, 
               cov.fac.impMI[,c(MIvars,93)],failcode=1) 
  
}
if(hor == 10){ 
  cov.fac.impMI$hicholrpXsyst <- cov.fac.impMI$HICHOLRP*cov.fac.impMI$SYST 
  cov.fac.impMI$diabXhtntrt2 <- cov.fac.impMI$DIAB*cov.fac.impMI$HTNTRT_2
  MIfit <- crr(cov.fac.impMI$MItime, cov.fac.impMI$statMI, 
               cov.fac.impMI[,c(MIvars,93:94)],failcode=1) 
  
}
if(hor == 5){  
  cov.fac.impMI$htntrt1Xqsmk12 <- cov.fac.impMI$HTNTRT_1 *cov.fac.impMI$qsmokage2_12 
  cov.fac.impMI$diabXsyst <- cov.fac.impMI$DIAB*cov.fac.impMI$SYST
  MIfit <- crr(cov.fac.impMI$MItime, cov.fac.impMI$statMI, 
               cov.fac.impMI[,c(MIvars,93:94)],failcode=1) 
}

if(hor == 15){ 
  cov.fac.impstrk$hicholXhtntrt2 <- cov.fac.impstrk$HTNTRT_2 *cov.fac.impstrk$HICHOLRP 
  strkfit <- crr(cov.fac.impstrk$strktime, cov.fac.impstrk$statstrk,  
                 cov.fac.impstrk[,c(strkvars,93)],failcode=1) 
}

if(hor == 10 ){  
  cov.fac.impstrk$atrialfbXhtntrt1 <- cov.fac.impstrk$ATRIALFB*cov.fac.impstrk$HTNTRT_1
  strkfit <- crr(cov.fac.impstrk$strktime, cov.fac.impstrk$statstrk,  # should this include ooph?
                 cov.fac.impstrk[,c(strkvars,93)],failcode=1) 
}

if(hor == 5){   
  cov.fac.impstrk$HTNTRT2Xqsmk4 <- cov.fac.impstrk$HTNTRT_2*cov.fac.impstrk$qsmokage2_4
  strkfit <- crr(cov.fac.impstrk$strktime, cov.fac.impstrk$statstrk,  # should this include ooph?
             cov.fac.impstrk[,c(strkvars,93)],failcode=1) 
}

lcfit <- crr(cov.fac.implc$lctime, cov.fac.implc$statlc, 
               cov.fac.implc[,lcvars],failcode=1) 
bcfit <- crr(cov.fac.impbc$bctime, cov.fac.impbc$statbc, 
               cov.fac.impbc[,bcvars],failcode=1) 
crcfit <- crr(cov.fac.impcrc$crctime, cov.fac.impcrc$statcrc, 
               cov.fac.impcrc[,crcvars],failcode=1) 
hipfit <- crr(cov.fac.imphip$hiptime, cov.fac.imphip$stathip, 
                cov.fac.imphip[,hipvars],failcode=1) 

if(hor == 15){  
  cov.fac.impdeath$hicholXhtntrt2 <- cov.fac.impdeath$HTNTRT_2 *cov.fac.impdeath$HICHOLRP 
  cov.fac.impdeath$diabXhtntrt2 <- cov.fac.impdeath$DIAB*cov.fac.impdeath$HTNTRT_2
  deathform15 = as.formula(paste0("Surv(deathtime, statdeath) ~ ",
                                 paste(names(cov.fac.impdeath)[
                                   c(deathvars,93:94)], collapse="+")))
  deathfit <- coxph(deathform15, data = cov.fac.impdeath)
}

if(hor == 10){
  cov.fac.impdeath$hicholrpXsyst <- cov.fac.impdeath$HICHOLRP*cov.fac.impdeath$SYST 
  cov.fac.impdeath$diabXhtntrt2 <- cov.fac.impdeath$DIAB*cov.fac.impdeath$HTNTRT_2
  cov.fac.impdeath$atrialfbXhtntrt1 <- cov.fac.impdeath$ATRIALFB*cov.fac.impdeath$HTNTRT_1
  deathform10 = as.formula(paste0("Surv(deathtime, statdeath) ~ ",
                                  paste(names(cov.fac.impdeath)[
                                    c(deathvars,93:95)], collapse="+")))
  deathfit <- coxph(deathform10, data = cov.fac.impdeath)
}

if(hor == 5){
  cov.fac.impdeath$htntrt1Xqsmk12 <- cov.fac.impdeath$HTNTRT_1 *cov.fac.impdeath$qsmokage2_12 
  cov.fac.impdeath$diabXsyst <- cov.fac.impdeath$DIAB*cov.fac.impdeath$SYST
  cov.fac.impdeath$HTNTRT2Xqsmk4 <- cov.fac.impdeath$HTNTRT_2*cov.fac.impdeath$qsmokage2_4
  deathform5 = as.formula(paste0("Surv(deathtime, statdeath) ~ ",
                                 paste(names(cov.fac.impdeath)[
                                   c(deathvars,93:95)], collapse="+")))
  deathfit <- coxph(deathform5, data = cov.fac.impdeath)
}



# save this output to a .csv
# 5 cols for each of 7 output
# coef, se(coef), 95% CI lb, ub, spacer
fonr <- length(deathvars)+3
fitout <- as.data.frame(matrix(NA, nrow = fonr, ncol=5*7))
rownames(fitout) <- names(cov.fac.impdeath)[c(deathvars,93:95)]
colnames(fitout) <- rep(c("HR", "se(HR)", "lb", "ub", ""),7)
fitout[c(MIvars,fonr-2,fonr-1), 1:4] <- cbind(summary(MIfit)$coef[,2:3], summary(MIfit)$conf.int[,3:4])
fitout[c(strkvars,fonr), 6:9] <- cbind(summary(strkfit)$coef[,2:3], summary(strkfit)$conf.int[,3:4])
fitout[lcvars, 11:14] <- cbind(summary(lcfit)$coef[,2:3], summary(lcfit)$conf.int[,3:4])
fitout[bcvars, 16:19] <- cbind(summary(bcfit)$coef[,2:3], summary(bcfit)$conf.int[,3:4])
fitout[crcvars, 21:24] <- cbind(summary(crcfit)$coef[,2:3], summary(crcfit)$conf.int[,3:4])
fitout[hipvars, 26:29] <- cbind(summary(hipfit)$coef[,2:3], summary(hipfit)$conf.int[,3:4])
fitout[c(deathvars,(fonr-2):fonr), 31:34] <- cbind(summary(deathfit)$coef[,2:3], summary(deathfit)$conf.int[,3:4])
## added code to convert the coefs for age, etc so interp
# tepiwk, tminwk, pulse30-syst
# take exp(C*log(HR) +/- 1.96*C*se); C = 1 SD
sdconv <- apply(tdat[c(65:66,78:83)],2,sd,na.rm=TRUE)
for(var in 1:7){
  coefnew <- log(fitout[c(65:66,78:83),(var-1)*5+1])*sdconv
  senew <- fitout[c(65:66,78:83),(var-1)*5+2]*sdconv
  fitout[c(65:66,78:83),(var-1)*5+1] <- exp(coefnew)
  fitout[c(65:66,78:83),(var-1)*5+3] <- exp(coefnew-1.96*senew)
  fitout[c(65:66,78:83),(var-1)*5+4] <- exp(coefnew+1.96*senew)
  rm(coefnew, senew)
}


###
# get each woman's predicted risk

# this gives the probability of an event at just the 10 year horizon
# for each woman
MIpred10 <- numeric(nrow(cov.fac.impMI))
indMI <- findInterval(hordy,MIfit$uftime)
for(i in 1:nrow(cov.fac.impMI)){  
  MIpred10[i] <- predict(MIfit,
    cov1=cov.fac.impMI[i,c(MIvars,93:94)])[indMI,-1] }
   
strkpred10 <- numeric(nrow(cov.fac.impstrk))
indstrk <- findInterval(hordy,strkfit$uftime)
for(i in 1:nrow(cov.fac.impstrk)){  
  strkpred10[i] <- predict(strkfit,
    cov1=cov.fac.impstrk[i,c(strkvars,93)])[indstrk,-1] }

lcpred10 <- numeric(nrow(cov.fac.implc))
indlc <- findInterval(hordy,lcfit$uftime)
for(i in 1:nrow(cov.fac.implc)){  
  lcpred10[i] <- predict(lcfit,
    cov1=cov.fac.implc[i,lcvars])[indlc,-1] }

bcpred10 <- numeric(nrow(cov.fac.impbc))
indbc <- findInterval(hordy,bcfit$uftime)
for(i in 1:nrow(cov.fac.impbc)){  
  bcpred10[i] <- predict(bcfit,
    cov1=cov.fac.impbc[i,bcvars])[indbc,-1] }

crcpred10 <- numeric(nrow(cov.fac.impcrc))
indcrc <- findInterval(hordy,crcfit$uftime)
for(i in 1:nrow(cov.fac.impcrc)){  
  crcpred10[i] <- predict(crcfit,
    cov1=cov.fac.impcrc[i,crcvars])[indcrc,-1] }

hippred10 <- numeric(nrow(cov.fac.imphip))
indhip <- findInterval(hordy,hipfit$uftime)
for(i in 1:nrow(cov.fac.imphip)){  
  hippred10[i] <- predict(hipfit,
    cov1=cov.fac.imphip[i,hipvars])[indhip,-1] }

deathpred10 <- numeric(nrow(cov.fac.impdeath))

inddeath <- findInterval(hordy,survfit(deathfit)$time)
cov.fac.impdeath = cov.fac.impdeathTE
for(i in 1:nrow(cov.fac.impdeathTE)){ 
  deathpred10[i] <- 1 - survfit(deathfit,newdata=cov.fac.impdeathTE[i,c(deathvars,93:95)])$surv[inddeath] }



rm(cov.fac.impdeath)
cov.fac.impdeath = cov.fac.impdeathTE

endpoint <- c("MI", "strk", "lc", "bc", "crc", "hip", "death")
predvec <- cbind(MIpred10, strkpred10, lcpred10,  bcpred10, crcpred10,
                hippred10, deathpred10)  # only need the 1-deathpred15 if I ran old code above

covmatvec <- expression(cov.fac.impMI,cov.fac.impstrk,cov.fac.implc,
                          cov.fac.impbc,cov.fac.impcrc,cov.fac.imphip,
                          cov.fac.impdeath)

kmlb <- c(0.85, 0.85, 0.85, 0.85, 0.95, 0.85, 0.5)

for(ep in 1:length(endpoint)){
  ept <- endpoint[ep]
  predtmp <- predvec[,ep]
  sorted <- sort(predtmp)
  pctl <- ecdf(sorted)(sorted)
  covmat <- eval(covmatvec[ep])

  pdf(paste0("calib10Ctest_",ept,".pdf"))
  plot(pctl, sorted, ylab="Predicted risk", xlab="Risk percentile", type="l",
     main=paste0("Calibration of ", ept, " model"))
  decl <- quantile(predtmp,seq(0,1,by=0.1))
  deccat <- cut(predtmp,decl,include.lowest=TRUE)
  decrate <- aggregate(covmat[,91]==1, by=list(deccat), FUN=mean)
  points(seq(0.05, 0.95, by=0.1),decrate[,2], col="blue", pch=19)
  legend("topleft", legend=c("Predicted risk", "Observed event rate by decile of predicted risk"),
         lty=c(1,NA), pch=c(NA,19), col=c("black", "blue"))
  # label the points as the actual event rate for women in each decile of predicted risk,
      # the line as the predicted risk
  dev.off()

  pdf(paste0("hist10Ctest_",ept,".pdf"))
  hist(predtmp, xlab="Predicted risk", breaks=40,
       main=paste0("Histogram of predicted risk for ", ept))
  dev.off()
  
  ## add in a KM plot showing the separation of deciles/quintiles

  decl10 <- quantile(predtmp,seq(0,1,by=0.2))
  qtlcat <- cut(predtmp,decl10,include.lowest=TRUE)
  
  tmp <- survfit(Surv(covmat[,90], covmat[,91]==1)~qtlcat)
  KMmat <- as.data.frame(matrix(NA, nrow=nrow(covmat), ncol=3))
  colnames(KMmat) <- c("surv", "time", "strat")
  KMmat[1:length(tmp$surv),1] <- tmp$surv
  KMmat[1:length(tmp$surv),2] <- tmp$time
  tmpind = cumsum(tmp$strata)
  KMmat[1:tmpind[1],3] <- 1
  KMmat[(tmpind[1]+1):tmpind[2],3] <- 2
  KMmat[(tmpind[2]+1):tmpind[3],3] <- 3
  KMmat[(tmpind[3]+1):tmpind[4],3] <- 4
  KMmat[(tmpind[4]+1):tmpind[5],3] <- 5
  KMmat$strat=factor(KMmat$strat)
  
  pdf(paste0("KM10Ctest_",ept,".pdf"))
  print(ggplot(KMmat, aes(x=time, y=surv, group=strat)) +
    geom_line(aes(color=strat)) + xlab("Time (years)") + ylab("Survival") + 
    theme_bw() + ylim(c(kmlb[ep],1)) + 
    labs(title=paste0("Survival by predicted risk quintile for ",ept)) +
    scale_x_continuous(breaks=seq(0,hordy,by=365.25), labels=0:10) +
    scale_color_discrete(guide = guide_legend(title = "Risk\nQuintile")) )
  dev.off()
  
 
}


CindMI <- cindex(as.matrix(MIpred10), 
               formula=as.formula(paste0("Hist(MItime, statMI)~",
                                         paste(names(cov.fac.impMI)[c(MIvars,93:94)],
                                               collapse="+"))),  
               cens.model="cox", data=cov.fac.impMI,eval.times=hordy,cause=1)
print("MI")
CindMI$AppCindex    

Cindstrk <- cindex(as.matrix(strkpred10), 
                 formula=as.formula(paste0("Hist(strktime, statstrk)~",
                                           paste(names(cov.fac.impstrk)[c(strkvars,93)],
                                                 collapse="+"))),  
                 cens.model="cox",  
                 data=cov.fac.impstrk,eval.times=hordy,cause=1)  
                
print("stroke")
Cindstrk$AppCindex    


Cindlc <- cindex(as.matrix(lcpred10), 
                 formula=as.formula(paste0("Hist(lctime, statlc)~",
                                           paste(names(cov.fac.implc)[lcvars],
                                                 collapse="+"))),  
                 cens.model="cox",  
                 data=cov.fac.implc,eval.times=hordy,cause=1)
                
print("lc")
Cindlc$AppCindex   

Cindbc <- cindex(as.matrix(bcpred10),  
                   formula=as.formula(paste0("Hist(bctime, statbc)~",
                                             paste(names(cov.fac.impbc)[bcvars],
                                                   collapse="+"))),  
                   cens.model="cox",  
                   data=cov.fac.impbc,eval.times=hordy,cause=1)
                  
print("bc")
Cindbc$AppCindex   


Cindcrc <- cindex(as.matrix(crcpred10),  
                     formula=as.formula(paste0("Hist(crctime, statcrc)~",
                                               paste(names(cov.fac.impcrc)[crcvars],
                                                     collapse="+"))),  
                     cens.model="cox",  
                     data=cov.fac.impcrc,eval.times=hordy,cause=1)
                     
print("crc")
Cindcrc$AppCindex   


Cindhip <- cindex(as.matrix(hippred10),  
                     formula=as.formula(paste0("Hist(hiptime, stathip)~",
                                               paste(names(cov.fac.imphip)[hipvars],
                                                     collapse="+"))),  
                     cens.model="cox",  
                     data=cov.fac.imphip,eval.times=hordy,cause=1)
print("hip")
Cindhip$AppCindex    


Cinddeath <- rcorr.cens(deathpred10, Surv(cov.fac.impdeath$deathtime, 
                                          cov.fac.impdeath$statdeath))
print("death")
Cinddeath[1]   
1- Cinddeath[1]   




#### set up for shiny app

###############
## creating histogram for the figure - will output the bins

load("testpred10C.Rdata")
bc10 = bcpred10
crc10 = crcpred10
death10 = deathpred10
hip10 = hippred10
lc10 = lcpred10
MI10 = MIpred10
strk10 = strkpred10

load("trpred10C.Rdata")  
bc10 = c(bc10,bcpred10)
crc10 = c(crc10,crcpred10) 
death10 = c(death10,deathpred10)
hip10 = c(hip10,hippred10)
lc10 = c(lc10,lcpred10)
MI10 = c(MI10,MIpred10)
strk10 = c(strk10,strkpred10)

bc10Ch = hist(bc10,breaks=seq(0,1,by=0.001))
crc10Ch = hist(crc10,breaks=seq(0,1,by=0.001))
death10Ch = hist(death10,breaks=seq(0,1,by=0.001))
hip10Ch = hist(hip10,breaks=seq(0,1,by=0.001))
lc10Ch = hist(lc10,breaks=seq(0,1,by=0.001))
MI10Ch = hist(MI10,breaks=seq(0,1,by=0.001))
strk10Ch = hist(strk10,breaks=seq(0,1,by=0.001))

save(bc10Ch,crc10Ch,death10Ch,hip10Ch,lc10Ch,MI10Ch,strk10Ch,file="pred10Chist.Rdata")


load("testpred5C.Rdata")
bc5 = bcpred5
crc5 = crcpred5
death5 = deathpred5
hip5 = hippred5
lc5 = lcpred5
MI5 = MIpred5
strk5 = strkpred5

load("trpred5C.Rdata")  
bc5 = c(bc5,bcpred5)
crc5 = c(crc5,crcpred5) 
death5 = c(death5,deathpred5)
hip5 = c(hip5,hippred5)
lc5 = c(lc5,lcpred5)
MI5 = c(MI5,MIpred5)
strk5 = c(strk5,strkpred5)

bc5Ch = hist(bc5,breaks=seq(0,1,by=0.001))
crc5Ch = hist(crc5,breaks=seq(0,1,by=0.001))
death5Ch = hist(death5,breaks=seq(0,1,by=0.001))
hip5Ch = hist(hip5,breaks=seq(0,1,by=0.001))
lc5Ch = hist(lc5,breaks=seq(0,1,by=0.001))
MI5Ch = hist(MI5,breaks=seq(0,1,by=0.001))
strk5Ch = hist(strk5,breaks=seq(0,1,by=0.001))

save(bc5Ch,crc5Ch,death5Ch,hip5Ch,lc5Ch,MI5Ch,strk5Ch,file="pred5Chist.Rdata")



load("testpred15C.Rdata")
bc15 = bcpred15
crc15 = crcpred15
death15 = deathpred15
hip15 = hippred15
lc15 = lcpred15
MI15 = MIpred15
strk15 = strkpred15

load("trpred15C.Rdata")  
bc15 = c(bc15,bcpred15)
crc15 = c(crc15,crcpred15) 
death15 = c(death15,1-deathpred15)
hip15 = c(hip15,hippred15)
lc15 = c(lc15,lcpred15)
MI15 = c(MI15,MIpred15)
strk15 = c(strk15,strkpred15)

bc15Ch = hist(bc15,breaks=seq(0,1,by=0.001))
crc15Ch = hist(crc15,breaks=seq(0,1,by=0.001))
death15Ch = hist(death15,breaks=seq(0,1,by=0.001))
hip15Ch = hist(hip15,breaks=seq(0,1,by=0.001))
lc15Ch = hist(lc15,breaks=seq(0,1,by=0.001))
MI15Ch = hist(MI15,breaks=seq(0,1,by=0.001))
strk15Ch = hist(strk15,breaks=seq(0,1,by=0.001))

save(bc15Ch,crc15Ch,death15Ch,hip15Ch,lc15Ch,MI15Ch,strk15Ch,file="pred15Chist.Rdata")


###############




