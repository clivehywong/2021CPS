#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(lme4)))
suppressMessages(suppressWarnings(library(lmerTest)))
suppressMessages(suppressWarnings(library(glmmLasso)))
suppressMessages(suppressWarnings(library(progress)))

args<-commandArgs(TRUE)
optMod <- args[1]
scriptstartedat<-Sys.time()

#### Timing ####
dhms <- function(label=""){
  if (!exists("scriptstartedat")) {
    scriptstartedat<<-Sys.time()
  } else {
    t<-as.integer(Sys.time()-scriptstartedat)
    cat(paste(label, t %/% (60*60*24) 
              ,paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
                     ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
                     ,formatC(t %% 60, width = 2, format = "d", flag = "0")
                     ,sep = ":"
              ),"\n"
    )
    )
  }
}


#### Read Data ####
initrtco <- function(m,g=NULL) {
  print("Lasso.R: Prepare RT and COPE data")
  suppressMessages(roibeta <- read_delim("roi/PeakCopeAV23.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
  suppressMessages(roi <- read_delim(sprintf("param.%s.roi",m), " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
  colnames(roibeta)<-c("Group","Did","Mod","Cont","ROI","Beta")
  roibeta<-roibeta[,1:6]
  roi<-roi$X1
  roibeta<-roibeta[roibeta$ROI%in%roi,]
  
  roibeta$Group <- factor(roibeta$Group)
  if (!is.null(g) && g!=0) {
    roibeta<-roibeta[roibeta$Group %in% g,]
    ngroup<<-1
  } else {
    ngroup<<-nlevels(roibeta$Group)
  }
  roibeta$Did <- factor(roibeta$Did)
  roibeta$Mod <- factor(roibeta$Mod)
  roibeta$Cont <- factor(roibeta$Cont)
  roibeta$ROI <- factor(roibeta$ROI)
  load("/bcs2/projects/cw/ps15/data/RT.o1boot.RData")

  S1.rt <- RT[RT$Did %in% levels(roibeta$Did),]
  S1.rt <- S1.rt[order(S1.rt$Did),]
  
  # cast #
  S1.rt.23 <- melt(S1.rt[,c("Did","Group","ACm","AIm","VCm","VIm")],id=c("Did","Group"))
  colnames(S1.rt.23) <- c("Did","Group","RTL","RT")
  S1.rt.11 <- melt(S1.rt[,c("Did","Group","ANm","VNm")],id=c("Did","Group"))
  colnames(S1.rt.11) <- c("Did","Group","CTL","CT")
  S1.rt.11 <- rbind(S1.rt.11[S1.rt.11$CTL=="ANm",],S1.rt.11[S1.rt.11$CTL=="ANm",],S1.rt.11[S1.rt.11$CTL=="VNm",],S1.rt.11[S1.rt.11$CTL=="VNm",])
  S1.rt.23 <- data.frame(S1.rt.23,S1.rt.11$CT)
  colnames(S1.rt.23) <- c("Did","Group","RTL","RT","CT")
  rtco <- dcast(droplevels(roibeta[roibeta$Cont %in% c(2,3),]), Mod + Cont + Did + Group ~ ROI , mean, value.var = "Beta")
  rtco <- data.frame(S1.rt.23,rtco)
  rtco$Did <- as.factor(rtco$Did)
  rtco$Group <- as.factor(rtco$Group)
  rtco$RT <- scale(rtco$RT)
  rtco$RCT <- scale(residuals(lmer(RT~CT+(1|RTL),data = rtco)))
  rtco <- rtco %>% group_by(RTL) %>% mutate(RRT = scale(residuals(lm(RT~CT))))
  rtco <- as.data.frame(rtco)
  
  # Normalize data
  rtco <- cbind(rtco[,-grep("_",colnames(rtco))],scale(rtco[,grep("_",colnames(rtco))]))

  #https://stats.stackexchange.com/questions/17336/how-exactly-does-one-control-for-other-variables
  #Introductory Econometrics: A Modern Approach by Jeffrey M. Wooldridge
  rtco <<- rtco
  rtco.orig <<- rtco
  
  # limited to 3 highest peaks
  x.co <<- colnames(rtco)
  x.co <<- x.co[grep("^rann",x.co)]
  if (length(x.co)<3) cat("Number of variables < 3")
  #x.co <<- x.co[grep("^P0.[123]",x.co)]
}

loadpath <- function(p, scale=TRUE) {
  print("Lasso.R: Load Path Coef")
  pathcoef=paste0(p,"/zMat.RData")
  load(pathcoef)
  zMat <- as.data.frame(zMat)
  print(roiNames)
  pathNames=sprintf("%s_%s",rep(roiNames,length(roiNames)), rep(roiNames,each=length(roiNames)))
  colnames(zMat) <- pathNames
  zMat<-zMat[,rep(roiNames,length(roiNames)) != rep(roiNames,each=length(roiNames))]
  x.pall <<- colnames(zMat)
  n <- length(unique(rtco.orig$Did))
  rtcopa<<-cbind(rtco.orig,zMat)
  rtcopa.orig<<-rtcopa
  print("Caution: Check if subject and task sequence match!!!")
  if (scale) rtcopa <<- cbind(rtcopa[,-grep("T_",colnames(rtcopa))],scale(rtcopa[,grep("T_",colnames(rtcopa))]))
  cat("Dimension:", dim(rtcopa)[1], "x", dim(rtcopa)[2],"\n")
  if (TRUE && exists("sigP")) {
    bf<-length(sigP)-nrow(sigP)
    sigP<-sigP*(bf) # Bonfferoni
    diag(sigP)<-NA # self-loop
    for (f in c(0.001,0.01,0.05,0.05*bf,1)) {
      tmpP<-sigP
      tmpP[tmpP>f]<-NA #threshold
      if (sum(!is.na(tmpP))>5) break; # make sure there is some variables
    }
    cat(sprintf("Filter Bonferoni corrected p <=%s\n",f))
    sigP<-tmpP #threshold
    x.pa <- pathNames[!is.na(c(sigP))]
  } else {
    cat("No filter for P\n")
    x.pa<-pathNames
  }
  x.pa<<-x.pa[order(x.pa)]
}

#### Helper Functions ####
sinkini <- function(...){
  if (!interactive()) {
    sink(...)
  }
}
stopQuietly <- function() {
  if (! interactive()) q(save="no")
}
selpath <- function(p,r) {
  return(p[ sub("T_(.*)_T_(.*)","\\1",p) %in% r & sub("T_(.*)_T_(.*)","\\2",p) %in% r ])
}
simVar <- function(x) {
  x<-x[grep("_",x)]
  x[grep("^T",x)] <- sub("^T_P.*_(.*)_T_P.*_(.*)","\\1-\\2",x[grep("^T",x)])
  x[grep("^P",x)] <- sub("^P.*_(.*)","\\1", x[grep("^P",x)])
  return(x)
}
#### Stat Function ####
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('<%s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}
variable.names.merMod <- function(fit) rownames(summary(fit)$coefficients)

r.squared <- function(...) UseMethod("r.squared")
r.squared.merMod <- function(fit) MuMIn::r.squaredGLMM(fit)[,"R2m"]
parens <- function(x) paste0("(",x,")")


posthoc.lmerTest <- function(fit,fmlbase){
  print("posthoc: lmerTest")
  fit <- as(fit,"merModLmerTest")
  kx<-variable.names(fit)
  kx<-kx[grep("_",kx)]
  if (length(kx)<1) {
    print("No variables left!")
    return(NULL)
  }
  print(sprintf("Fitted: %s", paste(simVar(kx),collapse = ", ")))
  s<-summary(fit)
  co<-data.frame(s$coefficients)
  co<-co[grep("_",rownames(co)),]
  #ci<-data.frame(confint(fit,method="boot",nsim=500)); ci<-ci[grep("_",rownames(ci)),]
  ci<-data.frame(confint(fit)); ci<-ci[grep("_",rownames(ci)),]
  sign1 <- function(n) sub("(.).*","\\1",sprintf("%+d",sign(n)))
  n<-paste(sign1(co[,1]),sub("-","_",simVar(rownames(co))), sep = "")
  # Cohen's f (Selya 2012 doi:10.3389/fpsyg.2012.00111)
  # save(fit,file="fit.RData")
  R2Full=r.squared(fit)
  f <- double(length = length(kx))
  for (i in 1:length(kx)) {
    fml2 <- as.formula(sprintf(fmlbase, paste(kx[-i], collapse = " + ")))
    fit2<-update(fit,fml2)
    R2=r.squared(fit2)
    if ( is.na(R2Full > R2) || R2Full > R2) f[i] = unname(sqrt((R2Full-R2)/(1-R2Full))) 
    else f[i] = 0
  }
  tmp <- data.frame(Name=kx,Sign=n,X=sub('.', '', n), beta=co[,1], Estimates=sprintf("%0.5f",co[,1]), 
                    lwr=sprintf("%0.5f",ci[,1]), upr=sprintf("%0.5f",ci[,2]), se=sprintf("%0.5f",co[,2]),
                    t=sprintf("%0.5f",co[,4]), p=pvalr(co[,5]), s=gtools::stars.pval(co[,5]), rawp=co[,5],
                    f=sprintf("%0.5f",f))
  tmp <- tmp[order(tmp$rawp),]
  return(tmp)
}
#### Main Function ####

runLasso <- function(y,x,data,output) {
  cat("runLasso:\n")

  # print(colnames(data))
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  { # Run Lasso
    lambda <- seq(100,0,by=-1)
    BIC_vec<-rep(Inf,length(lambda))
    family = gaussian(link = "identity")
    # specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
    # cat("Running ",j,collapse=" ")
    fml<-as.formula(sprintf("%s ~ 1+%s",y,paste(x,collapse = "+"))) #glmmlasso
    cat("Number of predictors:",length(x),"\n")
    if (ngroup==1) {
      fmlbase<-paste(y, " ~ %s + (1|RTL)", sep = "")
      re=list(RTL=~1)
      Delta.start<-as.matrix(t(rep(0,1+length(x)+length(levels(data$RTL))))) # variables in formula, levels in random effect
    } else {
      fmlbase<-paste(y, " ~ %s + (1|RTL)", sep = "")
      data$RTL<-interaction(data$RTL,data$Group)
      re=list(RTL=~1)
      Delta.start<-as.matrix(t(rep(0,1+length(x)+ length(levels(data$RTL))))) # variables in formula, levels in random effect
    }
    Q.start<-0.1
    pb <- progress_bar$new(format=paste(output," :percent :elapsed"),total=length(lambda), clear=FALSE)
    for(j in 1:length(lambda)){
      pb$tick()
      sinkini(file="/dev/null")
      glm <- glmmLasso(fix=fml, rnd=re,
                       family = family, data = data, 
                       lambda=lambda[j], switch.NR=F,final.re=TRUE,
                       control = list(start=Delta.start[j,],q_start=Q.start[j]))
      sinkini()
      iv<-length(colnames(glm$Deltamatrix)[seq(2,length(x))][glm$Deltamatrix[glm$conv.step,seq(2,length(x))]!=0])
      if (j==1 && iv > 0) {print("Increase lambda!!!")}
      BIC_vec[j]<-glm$aic
      Delta.start<-rbind(Delta.start,glm$Deltamatrix[glm$conv.step,])
      Q.start<-c(Q.start,glm$Q_long[[glm$conv.step+1]])
    }
    opt<-which.min(BIC_vec)
    glm_final <- glmmLasso(fix=fml, rnd=list(RTL=~1), family = family, data = data, lambda=lambda[opt],
                           switch.NR=F,final.re=TRUE, control = list(start=Delta.start[opt,],q_start=Q.start[opt]))
    
    #sinkini(file=paste(output,".lasso.txt",sep = ""))
    #print(summary(glm_final))
    #sinkini()
  } # Run Lasso
  
  { # Format result
    s<-summary(glm_final)
    p<-data.frame(X=rownames(s$coefficients),s$coefficients)
    p<-p[complete.cases(p),] # ### remove NA ### #
    p<-p[grep("_",p$X),]
    p$X <- as.character(p$X)
    p<-data.frame(Name=sub("-","_",simVar(p$X)),p)
    colnames(p)<-c("Name","X","Estimate","se","t","p")
    sign1 <- function(n) sub("(.).*","\\1",sprintf("%+d",sign(n)))
    p$Sign<-paste(sign1(p$Estimate),p$Name, sep = "")
    sig  <- p[p$p<0.05,]$X
    kept <- p$X
    fvar <- list("sig"=sig,"kept"=kept)
  } # Format result
  
  { # Print lasso output
    sinkini(file=paste(output,".lasso.txt",sep = ""))
    p$s <- gtools::stars.pval(p$p)
    p<-p[order(p$p),c("Name","Estimate","se","t","p","s","Sign")]
    print(sprintf("Kept: %s",paste(p$Sign,collapse = ",")))
    print(sprintf("Sig: %s",paste(p[p$p<0.05,]$Sign,collapse = ",")))
    print(p[,c("Name","Estimate","se","t","p","s")])
    sinkini()
  } # Print output

  { # Calculate anova result
    fml0 <-as.formula(sprintf(fmlbase,1))
    fml1 <-as.formula(sprintf(fmlbase,paste(kept,collapse = "+")))
    fit0 <- lme4::lmer(fml0,data)
    fit1 <- lme4::lmer(fml1,data)
    sinkini(sprintf("%s.anova.txt",output))
    print(anova(fit1))
    print(anova(fit0,fit1))
    sinkini()
  }
  
  if (FALSE) { # Calculate result with GLHT
    s <- posthoc.glht(fit1,fmlbase)
    sinkini(sprintf("%s.glht.txt",output))
    print(sprintf("Final Variables: %s, [ %s ]",paste(s$Sign[s$rawp<0.05], collapse = ","),paste(s$Sign[s$rawp>=0.05], collapse = ",")))
    print(s[,c("X","f","Estimates","lwr","upr","p","s")], row.names = FALSE)
    print(sprintf("R-square= %0.4f", r.squared(fit1)))
    sinkini()
  }
  
  { # Calculate result with lmerTest
    s <- posthoc.lmerTest(fit1,fmlbase)
    sinkini(sprintf("%s.lmer.txt",output))
    if (is.null(s)) {
      print("Nothing left!")
    } else {
      print(sprintf("Final Variables: %s, [ %s ]",paste(s$Sign[s$rawp<0.05], collapse = ","),paste(s$Sign[s$rawp>=0.05], collapse = ",")))
      print(s[,c("X","f","Estimates","se","lwr","upr","p","s")], row.names = FALSE)
      print(sprintf("R-square= %0.4f", r.squared(fit1)))
    }
    sinkini()
  } # Calculate result with lmerTest
  return(fvar)
}

modelDV <- function() { # Use global variables dv1 dv2 dv3 optXXX
  obase<-paste(optMod,sep="")
  
  cvar=cvar1=cvar2=cvar3=x.co
  
  loadpath(obase, scale=TRUE)
  print(x.pa)
  if (optMatch=="A") {
    x1 <- x2 <- x.pa
  } else if (optMatch=="M") {
    x1 <- selpath(x.pa,cvar1)
    x2 <- selpath(x.pa,cvar2)
  }
  
  # PATH
  output <- paste(obase,"d",dv1,"P",sep = "_")
  pvar1 <- runLasso(dv1,x1,rtcopa,output)

  output <- paste(obase,"e",dv2,"P",sep = "_")
  pvar2 <- runLasso(dv2,x2,rtcopa,output)
  print("Done")
  stopQuietly()
}

#### Main Script ####
options(na.action = "na.fail")
options(width=1000)
dhms()
optMod <- args[1]

#### M1: RT, CT, RCT ####
# Cope -> Path -> Cope + Path
# Load RT and COPE data
initrtco(optMod,optGrp)
dhms("Read")
print(sprintf("Begin with: %s", paste(x.co, collapse = " ")))

dv1="RT"
dv2="RCT"
modelDV()
dhms("Done")

#https://stats.stackexchange.com/questions/77891/checking-assumptions-lmer-lme-mixed-models-in-r
#qqnorm(myModel.lme, ~ranef(., level=2))
