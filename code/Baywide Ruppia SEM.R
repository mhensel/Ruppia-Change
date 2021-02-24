#Baywide Ruppia SEMs
library(tidyverse); library(readxl); library(patchwork);library(beyonce)
library(lme4); library(MuMIn); library(DHARMa); library(piecewiseSEM); library(nlme); library(semPlot)

#data to use is here. maybe the 8515 data if needed? also could load decrease/increase data. Check Baywide_RuChange_models for that
Env_Var_nosumy1sp.CBP_WQ <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_nosumy1sp.CBP_WQ.csv") #more informative and not as unweildy. Has y1 means, summer y1, spring this year. (eliminates y1 spring, summer this year). 
Env_Var.tss <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var.tss.csv")

Env_Var_nosumy1sp.CBP_WQ
Env_Var.tss

####MASTER ENV AND BAYWIDE DATASET RmZone_Env #### make sure to run the NA code below this
RmZone_Env <- Env_Var_nosumy1sp.CBP_WQ %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) 

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
#is.na(RmZone_Env) <- RmZone_Env == "NaN"
RmZone_Env[is.nan(RmZone_Env)] <- 0
is.na(RmZone_Env) <- RmZone_Env == "Inf"
is.na(RmZone_Env) <- RmZone_Env == "-Inf"
RmZone_Env <- as.data.frame(RmZone_Env)


##With TSS data now!
RmZone_Env.tss <- Env_Var.tss %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) #%>%
#filter(denscomp.max > 5) #this gets rid of any zones that are really tiny. about 200 data points less. originally i thought this would clean things up but it actually weakens the SEM? idk about for these graphs

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
#is.na(RmZone_Env) <- RmZone_Env == "NaN"
RmZone_Env.tss[is.nan(RmZone_Env.tss)] <- 0
is.na(RmZone_Env.tss) <- RmZone_Env.tss == "Inf"
is.na(RmZone_Env.tss) <- RmZone_Env.tss == "-Inf"
RmZone_Env.tss <- as.data.frame(RmZone_Env.tss)


RmZone_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env.csv")
#this is the lmer that works best for reference
rmmeanenvtest.lmer <- lmer(dens.percomp.change ~ SAVArea.y1 + ChlA.me + Sal.me + (dens.percomp.y1:Sal.spme)+ (dens.percomp.y1:Sal.y1Dpos) + (dens.percomp.y1:ChlA.spme) + (TP.sumy1me:dens.percomp.y1) +(Secc.y1Dme:dens.percomp.y1) + (1|STATION),  data = RmZone_Env) 

g<-qplot(y = log10(TP.spme), x = log10(TSS.spme), color = year, data = RmZone_Env.TSS)+ theme(legend.position="none")
h<-qplot(y = dens.percomp.change, x = Sal.Dneg, color = year, data = RmZone_Env8515)
g+h

ggplot(data = RmZone_Env.TSS) + 
  #geom_smooth(method = "lm", aes(x = (log10(TSS.spme)*dens.percomp.y1), y = dens.percomp.change)) +
  geom_point(aes(x = (log10(TSS.spme)*dens.percomp.y1), y = dens.percomp.change, color = dens.percomp.y1)) + 
  ylab("Ruppia density change") + xlab("TSS.spring * Ruppia density y-1") +
  theme_bw(base_size=20) +
  scale_color_gradient(low="blue", high="red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

head(Env_Var_year.CBP_WQ)
#chris tss modeling
f1 <- lmer(TSS.me ~ TP.me + TN.me + Secc.me + Sal.me + (1|STATION) + (1|year), Env_Var_year.CBP_WQ[Env_Var_year.CBP_WQ$year > 2003,])

r.squaredGLMM(f1)
coT <- lm(Sal.spme ~ TSS.spme, data = RmZone_Env.TSS)
cor(RmZone_Env.TSS:log10(Sal.spme), RmZone_Env.TSSlog10(TSS.spme))

##Generate the Rm_SEM dataset here####
#drop NAs, filter out nonsense values

Rm_SEM <- RmZone_Env %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0) #%>% filter(year > 1999)
Rm_SEM.tss <- RmZone_Env.tss %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0)

Rm_SEM  #1455 obs
Rm_SEM.tss #8934 obs

###sems below, listed
RmBay.sem  #p = .21, spring means no logging, TN included makes it worse but we are going to keep it in
RmBay.sem1 #p = .34, same as above but TN out. 
RmBayLOG.sem # p = .53, everything logged in RmBay.sem
RmBayLMLOG.sem # p = .38, same as RmBayLOG.sem but linear model instead of ME
RmBayLM.sem #p = .8, same as RmBay.sem but linear

RmBayLOG_03.sem
RmBayLOG.tss.sem #p = .3, RmBayLOG with TSS (uses Rm_SEM.tss)

AIC(RmBayLOG.tss.sem, RmBayLOG_03.sem)

#this first one is the skeleton w p = .2.
#TN in model, spring means only
RmBaysp.sem <- psem(
  ChlAsp <- lme(ChlA.spme ~
                  TP.spme +
                  TN.spme,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccsp <- lme(Secc.spme ~
                  ChlA.spme +
                  Sal.spme +
                  TN.spme +
                  TP.spme,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 Sal.spme + 
                 ChlA.spme + 
                 TP.spme +
                 TN.spme +
                 Secc.spme +
                 (dens.percomp.y1:Sal.spme) + 
                 (dens.percomp.y1:ChlA.spme) + 
                 (TP.spme:dens.percomp.y1) +
                 (TN.spme:dens.percomp.y1) +
                 (Secc.spme:dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  TN.spme %~~% TP.spme,
  Secc.spme %~~% Sal.spme,
  ChlA.spme %~~% Sal.spme,
  TP.spme %~~% Sal.spme,
  data = Rm_SEM)

summary(RmBaysp.sem)
coefs(RmBay.sem)
dSep(RmBay.sem)
fisherC(RmBay.sem)

#annual means (not spring) fits p = .168
RmBay.sem <- psem(
  ChlAme <- lme(ChlA.me ~
                  TP.me +
                  TN.me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccme <- lme(Secc.me ~
                  ChlA.me +
                  Sal.me +
                  TN.me +
                  TP.me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 Sal.me + 
                 ChlA.me + 
                 TP.me +
                 TN.me +
                 Secc.me +
                 (dens.percomp.y1:Sal.me) + 
                 (dens.percomp.y1:ChlA.me) + 
                 (TP.me:dens.percomp.y1) +
                 (TN.me:dens.percomp.y1) +
                 (Secc.me:dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  TN.me %~~% TP.me,
  Secc.me %~~% Sal.me,
  ChlA.me %~~% Sal.me,
  TP.me %~~% Sal.me,
  data = Rm_SEM)

summary(RmBay.sem)

#y1 means (not spring) DOES NOT FIT
RmBayy1.sem <- psem(
  ChlAy1me <- lme(ChlA.y1me ~
                  TP.y1me +
                  TN.y1me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccy1me <- lme(Secc.y1me ~
                  ChlA.y1me +
                  Sal.y1me +
                  TN.y1me +
                  TP.y1me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 Sal.y1me + 
                 ChlA.y1me + 
                 TP.y1me +
                 TN.y1me +
                 Secc.y1me +
                 (dens.percomp.y1:Sal.y1me) + 
                 (dens.percomp.y1:ChlA.y1me) + 
                 (TP.y1me:dens.percomp.y1) +
                 (TN.y1me:dens.percomp.y1) +
                 (Secc.y1me:dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  TN.y1me %~~% TP.y1me,
  Secc.y1me %~~% Sal.y1me,
  ChlA.y1me %~~% Sal.y1me,
  TP.y1me %~~% Sal.y1me,
  data = Rm_SEM)

summary(RmBayy1.sem)

plot(RmBayLOG.sem, return = FALSE, node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white", width = .7),
     edge_attrs = data.frame(style = "solid", color = "black", arrowsize = 1, penwidth = 1),
     ns_dashed = T, alpha = 0.05, show = "std", digits = 2,
     add_edge_label_spaces = T)

##TN removed p = .31 no logs
RmBay.sem1 <- psem(
  ChlAsp <- lme(ChlA.spme ~
                  TP.spme ,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccsp <- lme(Secc.spme ~
                  ChlA.spme +
                  Sal.spme +
                  TP.spme,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 Sal.spme + 
                 ChlA.spme + 
                 TP.spme +
                 Secc.spme +
                 (dens.percomp.y1:Sal.spme) + 
                 (dens.percomp.y1:ChlA.spme) + 
                 (TP.spme:dens.percomp.y1) +
                 (Secc.spme:dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  Secc.spme %~~% Sal.spme,
  ChlA.spme %~~% Sal.spme,
  TP.spme %~~% Sal.spme,
  data = Rm_SEM)

summary(RmBay.sem1)
AIC(RmBay.sem, RmBay.sem1)

####BEST SEM LOG####
#Log everytnig SEM. Highest P value here
RmBayLOG.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
    #            correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
     #           correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
     #          correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  data = Rm_SEM)

summary(RmBayLOG.sem)
coefs(RmBayLOG.sem)
dSep(RmBayLOG.sem)
fisherC(RmBayLOG.sem)

#test .tss data
#ok i just wanted to test if the best SEM we have still fits the data from the .tss dataset and it does p = .23
RmBayLOGcompare.tss.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  data = Rm_SEM.tss)

summary(RmBayLOGcompare.tss.sem)

AIC(RmBayLOGsp.tss_03.sem, RmBayLOG.tss.sem, RmBayLOG.sem)

#linear logged using lms fits fine at .363
RmBayLMLOG.sem <- psem(
  ChlAsp <- lm(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                data = Rm_SEM),
  Seccsp <- lm(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                data = Rm_SEM),
  RuInt <- lm(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               data = Rm_SEM),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  data = Rm_SEM)

summary(RmBayLMLOG.sem)
coefs(RmBayLMLOG.sem)
dSep(RmBayLMLOG.sem)
fisherC(RmBayLMLOG.sem)

#even this basic ass bitch lm no logs fits .25
RmBayLM.sem <- psem(
  ChlAsp <- lm(ChlA.spme ~
                 TP.spme +
                 TN.spme,
               data = Rm_SEM),
  Seccsp <- lm(Secc.spme ~
                 ChlA.spme +
                 Sal.spme +
                 TN.spme +
                 TP.spme,
               data = Rm_SEM),
  RuInt <- lm(dens.percomp.change ~
                dens.percomp.y1 +
                Sal.spme + 
                ChlA.spme+ 
                TP.spme +
                TN.spme + 
                Secc.spme +
                (dens.percomp.y1:Secc.spme) +
                (dens.percomp.y1:Sal.spme) + 
                (dens.percomp.y1:ChlA.spme) + 
                (TP.spme:dens.percomp.y1) +
                (TN.spme:dens.percomp.y1) +
                (Secc.spme:dens.percomp.y1),
              data = Rm_SEM),
  TN.spme %~~% TP.spme,
  Secc.spme %~~% Sal.spme,
  ChlA.spme %~~% Sal.spme,
  TP.spme %~~% Sal.spme,
  data = Rm_SEM)

summary(RmBayLM.sem)

####TSS models#####

##TSS predicted in this one. Great p = .276, 
RmBayLOGsp.tss.sem <- psem(
  TSSsp <- lme(log10(TSSr.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
        #        correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TSSr.spme) +
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
           #     correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(TSSr.spme) +
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
          #      correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(TSSr.spme) +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                 (dens.percomp.y1:log10(TSSr.spme)) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
          #     correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss),
  log10(TSSr.spme) %~~% log10(Sal.spme),
  log10(TN.spme) %~~% log10(TP.spme),
  #log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  data = Rm_SEM.tss)

summary(RmBayLOGsp.tss.sem)
coefs(RmBayLOGsp.tss.sem)
dSep(RmBayLOGsp.tss.sem)
fisherC(RmBayLOGsp.tss.sem)

AIC(RmBayLOG_03.sem, RmBayLOG.tss.sem)

### BAY-WIDE STRUCTURAL EQUATION MODEL larfchunk #########

# Create SEM
station.sem <- psem(
  # Predictors of chl-a
  chla <- lme(log10.chla ~
                tp +
                tn +
                log10.WTEMP,
              random = ~ 1 | STATION/Year,
              correlation = corARMA(form = ~ 1 | STATION/Year, q = 2),
              control = lmeControl(opt = "optim"),
              data = sav_env),
  
  # Predictors of secchi
  secchi <- lme(log10.SECCHI ~
                  log10.WTEMP +
                  log10.chla +
                  tn +
                  tp,
                random = ~ 1 | STATION/Year,
                correlation = corARMA(form = ~ 1 | STATION/Year, q = 2),
                control = lmeControl(opt = "optim"),
                data = sav_env),
  
  Richness <- lme(Richness ~
                    log10.SAVCover1 +
                    log10.chla +
                    log10.SALINITY +
                    tp +
                    log10.SECCHI +
                    log10.Hectares1MPlusSAVComp,
                  random = ~ 1 | STATION/Year,
                  correlation = corARMA(form = ~ 1 | STATION/Year, q = 2),
                  control = lmeControl(opt = "optim"),
                  data = sav_env),
  
  # Predictors of SAV
  SAV <- lme(log10.SAVCover ~
               log10.SAVCover1 +
               log10.chla +
               log10.SECCHI +
               log10.WTEMP +
               log10.SALINITY +
               tn +
               tp +
               Richness +
               log10.Hectares1MPlusSAVComp,
             random = ~ 1 | STATION/Year,
             correlation = corARMA(form = ~ 1 | STATION/Year, q = 2),
             control = lmeControl(opt = "optim"),
             data = sav_env),
  
  tn %~~% tp,
  
  log10.SALINITY %~~% log10.SECCHI,
  
  data = sav_env
  
)

# Get summary
summary(station.sem)


#
##
###
####SUBESTUARIES#####
###
##
#
#
library(patchwork)
i<-qplot(x = log10(nptp.adj)*dens.percomp.y1, y = dens.percomp.change, color = Year, data = RmSubE988_SEM)+ theme(legend.position="none")
j<-qplot(x = log10(nptn.adj)*dens.percomp.y1, y = dens.percomp.change, color = Year, data = RmSubE988_SEM)+ theme(legend.position="none")
i+j

rmSubEenv.lmer <- lmer(dens.percomp.change ~ 
                         SAVArea.y1 + 
                         log10(Developed.y1) +  log10(Agro.y1) + 
                         log10(tssx.y1) + #log10(tssx.adj)
                         log10(fertilizerp_kg.y1) + log10(fertilizern_kg.y1) +
                         log10(nptp.y1) + log10(nptn.y1) +
                         log10(Developed.y1):dens.percomp.y1 + log10(Agro.y1):dens.percomp.y1 + 
                         log10(tssx.y1):dens.percomp.y1 + log10(tssx.y1):dens.percomp.y1 +
                         log10(fertilizerp_kg.y1):dens.percomp.y1 + 
                         log10(fertilizern_kg.y1):dens.percomp.y1 +
                         log10(nptn.y1):dens.percomp.y1 + log10(nptp.y1):dens.percomp.y1 + 
                         (1|SUBEST_ID),  data = RmSubE_Env) 
Anova(rmSubEenv.lmer)
dens.percomp.change ~ SAVArea.y1 +
  Developed.y1 + 
  log10(tssx.y1) + 
  log10(nptp.y1) + log10(nptn.y1) +
  ptp.y1 + ptn.y1 +
  Developed.y1:dens.percomp.y1 + 
  log10(tssx.y1):dens.percomp.y1 + 
  log10(nptn.y1):dens.percomp.y1 + log10(nptp.y1):dens.percomp.y1 + 
  ptn.y1:dens.percomp.y1 + ptp.y1:dens.percomp.y1 + 
  (1|SUBEST_ID)

#SEM Notes: 
#dont worry about the NAs here (in that, they arent removing any important data). some Subestuaries that i have Ruppia for just dont have any data
#if ptp or ptn are in the SEM, may need to filter out the 0s bc they mess up the log10(). not that many points tho
RmSubE_SEM <- RmSubE_Env %>% drop_na() 
%>% filter(tssx.adj < 10000000)
RmSubE988_SEM <- RmSubE_Env982 %>% drop_na() ##USE THIS ONE

%>% filter(ptn.adj > 0) %>% filter(tssx.adj > 0) %>% filter(tssx.adj < 10000000,) %>% filter(ptp.adj > 0) %>% filter(ptn.adj > 0) %>% filter(ptn.y1 > 0) #%>% filter(year > 1999)
#filtering the ptn takes out some data that fucked the SEM so idk whats up here. 


#Rm_SEM.tss <- RmZone_Env.TSS %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0) %>% filter(year > 1999)

RmSubE_SEM  #868
s

###sems below, listed
RmSubENoLogtrim.sem  #p = .3 flow in there, no predicting man/fert, p/n together np

RmBay.sem1 #p = .34, same as above but TN out. 
RmBayLOG.sem # p = .53, everything logged in RmBay.sem
RmBayLMLOG.sem # p = .38, same as RmBayLOG.sem but linear model instead of ME
RmBayLM.sem #p = .8, same as RmBay.sem but linear

RmBayLOG_03.sem
RmBayLOG.tss.sem #p = .3, RmBayLOG with TSS (uses Rm_SEM.tss)

AIC(RmBayLOG.tss.sem, RmBayLOG_03.sem)



##start testing SEMS here 

#this unlogged doesnt fit anymore with the trimmed dataset .016
#FYI untrimmed and y1 does not work
#when ptn > 0 isnt filtered out, this doesnt fit anymoe
RmSubENoLogtrim.sem <- psem(
  tss <- lme(tssx.adj ~ flow + Agro ,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE988_SEM),
  
  nptp <- lme(nptp.adj ~ manurep_kg.adj + fertilizerp_kg.adj + flow + Agro,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  
  nptn <- lme(nptn.adj ~ manuren_kg.adj + fertilizern_kg.adj + flow + Agro,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 +
                    flow +
                    nptn.adj +
                    nptp.adj +
                    tssx.adj +
                    dens.percomp.y1:flow +
                    dens.percomp.y1:nptp.adj+
                    dens.percomp.y1:nptn.adj+
                    dens.percomp.y1:tssx.adj,
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  #manurep_kg.adj %~~% manuren_kg.adj,
  fertilizerp_kg.adj %~~% fertilizern_kg.adj,
  nptn.adj %~~% nptp.adj,
  nptn.adj %~~% tssx.adj, 
  nptp.adj %~~% tssx.adj,
  #log10(nptn.y1) %~~% log10(manuren_kg.y1),
  data = RmSubE988_SEM)
summary(RmSubENoLogtrim.sem)



#logged version of fits model below
#fits! p = .09, untrimmed doesnt fit
#worse w Developed. doesnt fit if flow is logged. doesnt fit w manure in it

RmSubELOGtrim.sem <-  psem(
    tss <- lme(log10(tssx.adj) ~ flow + Agro,
               random = ~ 1 | SUBEST_ID,
               correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
               data = RmSubE988_SEM),
    fertilizern <- lme(log10(fertilizern_kg.adj)  ~ Agro,
                       random = ~ 1 | SUBEST_ID,
                       correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                       data = RmSubE988_SEM),
    fertilizerp <- lme(log10(fertilizerp_kg.adj)  ~ Agro,
                       random = ~ 1 | SUBEST_ID,
                       correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                       data = RmSubE988_SEM),
    
    nptp <- lme(log10(nptp.adj) ~ log10(fertilizerp_kg.adj) + flow,
                random = ~ 1 | SUBEST_ID,
                correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                data = RmSubE988_SEM),
    
    nptn <- lme(log10(nptn.adj) ~ log10(fertilizern_kg.adj)  + flow ,
                random = ~ 1 | SUBEST_ID,
                correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                data = RmSubE988_SEM),
    
    Ruchange <- lme(dens.percomp.change ~
                      dens.percomp.y1 +
                      flow +
                      log10(nptn.adj) +
                      log10(nptp.adj) +
                      log10(tssx.adj) +
                      dens.percomp.y1:flow +
                      dens.percomp.y1:log10(nptp.adj)+
                      dens.percomp.y1:log10(nptn.adj)+
                      dens.percomp.y1:log10(tssx.adj),
                    random = ~ 1 | SUBEST_ID,
                    correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                    data = RmSubE988_SEM), 
    log10(fertilizern_kg.adj) %~~% log10(fertilizerp_kg.adj),
    #log10(manuren_kg.adj) %~~% log10(manurep_kg.adj),
    log10(nptn.adj) %~~% log10(nptp.adj),
    log10(nptn.adj) %~~% log10(tssx.adj), 
    log10(nptp.adj) %~~% log10(tssx.adj),
    data = RmSubE988_SEM)
summary(RmSubELOGtrim.sem)

AIC(RmSubELOGtrim.sem, RmSubELOGtrimflowx.sem)

#just N fits p = .84
#note- manure and fertilizer breaks all these models
RmSubELOGtrim.N.sem <-  psem(
  tss <- lme(log10(tssx.adj) ~ flow + Agro,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE988_SEM),
  
  nptn <- lme(log10(nptn.adj) ~  flow + Agro ,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 +
                    flow +
                    log10(nptn.adj) +
                    log10(tssx.adj) +
                    dens.percomp.y1:flow +
                    dens.percomp.y1:log10(nptn.adj)+
                    dens.percomp.y1:log10(tssx.adj),
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  log10(nptn.adj) %~~% log10(tssx.adj), 
  data = RmSubE988_SEM)
summary(RmSubELOGtrim.N.sem)

#this model fits well but nptp isnt sig p = .9
RmSubELOGtrim.P.sem <- psem(
  tss <- lme(log10(tssx.adj) ~ flow + Agro,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE988_SEM),
  
  nptp <- lme(log10(nptp.adj) ~ Agro + flow,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 +
                    flow +
                    log10(nptp.adj) +
                    log10(tssx.adj) +
                    dens.percomp.y1:flow +
                    dens.percomp.y1:log10(nptp.adj)+
                    dens.percomp.y1:log10(tssx.adj),
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  log10(nptp.adj) %~~% log10(tssx.adj),
  data = RmSubE988_SEM)
summary(RmSubELOGtrim.P.sem)

#flow*nptp/n in RuChange increases P to .7 but those terms arent sig. do that interaction term but take flow out of that model and its .619 and the NPTP and NPTN are stronger predictors
RmSubELOGtrimflowx.sem <-  psem(
  tss <- lme(log10(tssx.adj) ~ flow + Agro,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE988_SEM),
  
  nptp <- lme(log10(nptp.adj) ~ log10(fertilizerp_kg.adj) +log10(manurep_kg.adj)+ flow,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  
  nptn <- lme(log10(nptn.adj) ~ log10(fertilizern_kg.adj) + log10(manuren_kg.adj) + flow ,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 +
                    #flow +
                    flow*log10(nptn.adj) +
                    flow*log10(nptp.adj) +
                    log10(tssx.adj) +
                    #dens.percomp.y1:flow +
                    dens.percomp.y1:log10(nptp.adj)+
                    dens.percomp.y1:log10(nptn.adj)+
                    dens.percomp.y1:log10(tssx.adj),
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  log10(fertilizern_kg.adj) %~~% log10(fertilizerp_kg.adj),
  log10(manuren_kg.adj) %~~% log10(manurep_kg.adj),
  #log10(fertilizern_kg.adj) %~~% Agro,
  #log10(manuren_kg.adj) %~~% Agro,
  #Agro %~~% log10(fertilizerp_kg.adj),
  #Agro %~~% log10(manurep_kg.adj),
  log10(nptn.adj) %~~% log10(nptp.adj),
  log10(nptn.adj) %~~% log10(tssx.adj), 
  log10(nptp.adj) %~~% log10(tssx.adj),
  data = RmSubE988_SEM)
summary(RmSubELOGtrimflowx.sem)


#this unlogged y1 almost fits
RmSubEY1NoLog.sem <- psem(
  tss <- lme(tssx.y1 ~ flow.y1 ,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE988_SEM),
  
  nptp <- lme(nptp.y1 ~ manurep_kg.y1 + fertilizerp_kg.y1 + flow.y1,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE988_SEM),
  
  nptn <- lme(nptn.y1 ~ manuren_kg.y1 + fertilizern_kg.y1 + flow.y1,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             data = RmSubE988_SEM),
  
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 +
                    flow.y1 +
                    nptn.y1 +
                    nptp.y1 +
                    tssx.y1 +
                    dens.percomp.y1:flow.y1 +
                    dens.percomp.y1:nptp.y1+
                    dens.percomp.y1:nptn.y1+
                    dens.percomp.y1:tssx.y1,
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  manurep_kg.y1 %~~% fertilizerp_kg.y1,
  nptn.y1 %~~% nptp.y1,
  nptn.y1 %~~% tssx.y1, 
  nptp.y1 %~~% tssx.y1,
  #log10(nptn.y1) %~~% log10(manuren_kg.y1),
  data = RmSubE988_SEM)
summary(RmSubEY1NoLog.sem)




RmSubENoLogLM.sem <- psem(
  tss <- lm(tssx.adj ~ logq0(flow) + Year +nptn.adj + nptp.adj,
             data = RmSubE988_SEM),
  nptn <- lm(nptn.adj ~ manuren_kg.adj + fertilizern_kg.adj + flow,
              data = RmSubE988_SEM),
  nptp <- lm(nptp.adj ~ manurep_kg.adj + fertilizerp_kg.adj + flow,
             data = RmSubE988_SEM),

  #manureN <- lm(manuren_kg.adj ~ Agro + WatershedHa,
  #               data =RmSubE988_SEM),
  
  #fertilizerN <- lm(fertilizern_kg.adj ~ Agro + WatershedHa,
  #                   data =RmSubE988_SEM),
  
  Ruchange <- lm(dens.percomp.change ~
                    dens.percomp.y1 + log10(flow) +
                    #flow.adj*nptn.adj +
                    nptn.adj +
                    nptp.adj +
                    tssx.adj +
                    dens.percomp.y1:log10(flow) +
                    dens.percomp.y1:nptn.adj+
                    dens.percomp.y1:nptp.adj+
                    dens.percomp.y1:tssx.adj,
                  #log10.Hectares1MPlusSAVComp,
                  data = RmSubE988_SEM), 
  fertilizerp_kg.adj %~~% fertilizern_kg.adj,
  manurep_kg.adj %~~% manuren_kg.adj,
  nptn.adj %~~% nptp.adj,
  nptn.adj %~~% tssx.adj,
  nptp.adj %~~% tssx.adj,
  data = RmSubE988_SEM)
summary(RmSubENoLogLM.sem)





###Testing everything ####
RmSubELOG.sem <- psem(
  nptp <- lme(log10(nptp.adj) ~ Developed + log10(manurep_kg.adj) + log10(fertilizerp_kg.adj) + log10(WatershedHa),
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE_SEM),
  nptn <- lme(log10(nptn.adj) ~ Developed + log10(manuren_kg.adj) + log10(fertilizern_kg.adj) + log10(WatershedHa),
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE_SEM),
  
  ptp <- lme(log10(ptp.adj) ~ Developed + log10(WatershedHa),
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             data = RmSubE_SEM),
  
  ptn <- lme(log10(ptn.adj) ~ Developed + log10(WatershedHa),
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             data = RmSubE_SEM),
  
  manurep <- lme(log10(manurep_kg.adj) ~ Agro + log10(WatershedHa),
                 random = ~ 1 | SUBEST_ID,
                 correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                 data =RmSubE_SEM),
  manuren <- lme(log10(manuren_kg.adj) ~ Agro + log10(WatershedHa),
                 random = ~ 1 | SUBEST_ID,
                 correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                 data =RmSubE_SEM),
  
  fertilizerp <- lme(log10(fertilizerp_kg.adj) ~ Agro + log10(WatershedHa),
                     random = ~ 1 | SUBEST_ID,
                     correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                     data =RmSubE_SEM),
  fertilizern <- lme(log10(fertilizern_kg.adj) ~ Agro + log10(WatershedHa),
                     random = ~ 1 | SUBEST_ID,
                     correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                     data =RmSubE_SEM),
  
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 + Developed + Agro +
                    log10(flow)*log10(nptn.adj) +
                    log10(flow)*log10(nptp.adj) +
                    log10(fertilizerp_kg.adj) + log10(fertilizern_kg.adj) +
                    log10(ptp.adj) + log10(ptn.adj) +
                    log10(nptp.adj) + log10(nptn.adj) +
                    log10(tssx.adj) + log10(flow) +
                    dens.percomp.y1:flow +
                    #dens.percomp.y1:Agro + dens.percomp.y1:Developed +
                    dens.percomp.y1:log10(fertilizerp_kg.adj)+ dens.percomp.y1:log10(fertilizern_kg.adj)+ 
                    dens.percomp.y1:log10(ptp.adj)+ dens.percomp.y1:log10(ptn.adj)+
                    dens.percomp.y1:log10(nptp.adj) + dens.percomp.y1:log10(nptp.adj) +
                    dens.percomp.y1:log10(tssx.adj),
                  #log10.Hectares1MPlusSAVComp,
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE_SEM), 
  log10(ptn.adj) %~~% log10(fertilizern_kg.adj),
  #log10(ptp.adj) %~~% log10(ptn.adj),
  #log10(nptp.adj) %~~% log10(nptn.adj),
  log10(nptn.y1) %~~% log10(manuren_kg.y1),
  log10(ptp.adj) %~~% log10(fertilizerp_kg.adj),
  log10(fertilizerp_kg.adj) %~~% log10(fertilizern_kg.adj),
  data = RmSubE_SEM)

summary(RmSubELOG.sem)

coefs(RmSubELOG.sem)
dSep(RmSubELOG.sem)
fisherC(RmSubELOG.sem)

#same as above but Y1
RmSubELOGy1.sem <- psem(
  nptp <- lme(log10(nptp.y1) ~ Developed + log10(manurep_kg.y1) + log10(fertilizerp_kg.y1) + log10(WatershedHa),
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE_SEM),
  nptn <- lme(log10(nptn.y1) ~ Developed + log10(manuren_kg.y1) + log10(fertilizern_kg.y1) + log10(WatershedHa),
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE_SEM),
  
  ptp <- lme(log10(ptp.y1) ~ Developed + log10(WatershedHa),
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             data = RmSubE_SEM),
  
  ptn <- lme(log10(ptn.y1) ~ Developed + log10(WatershedHa),
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             data = RmSubE_SEM),
  
  manurep <- lme(log10(manurep_kg.y1) ~ Agro + log10(WatershedHa),
                 random = ~ 1 | SUBEST_ID,
                 correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                 data =RmSubE_SEM),
  manuren <- lme(log10(manuren_kg.y1) ~ Agro + log10(WatershedHa),
                 random = ~ 1 | SUBEST_ID,
                 correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                 data =RmSubE_SEM),
  
  fertilizerp <- lme(log10(fertilizerp_kg.y1) ~ Agro + log10(WatershedHa),
                     random = ~ 1 | SUBEST_ID,
                     correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                     data =RmSubE_SEM),
  fertilizern <- lme(log10(fertilizern_kg.y1) ~ Agro + log10(WatershedHa),
                     random = ~ 1 | SUBEST_ID,
                     correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                     data =RmSubE_SEM),
  
  Ruchange <- lme(dens.percomp.change ~
                    dens.percomp.y1 + Developed + Agro +
                    #log10(flow):log10(nptn.y1) +
                    #log10(flow):log10(nptp.y1) +
                    log10(fertilizerp_kg.y1) + log10(fertilizern_kg.y1) +
                    log10(ptp.y1) + log10(ptn.y1) +
                    log10(nptp.y1) + log10(nptn.y1) +
                    log10(tssx.y1) + log10(flow) +
                    dens.percomp.y1:flow +
                    dens.percomp.y1:Agro + dens.percomp.y1:Developed +
                    dens.percomp.y1:log10(fertilizerp_kg.y1)+ dens.percomp.y1:log10(fertilizern_kg.y1)+ 
                    dens.percomp.y1:log10(ptp.y1)+ dens.percomp.y1:log10(ptn.y1)+
                    dens.percomp.y1:log10(nptp.y1) + dens.percomp.y1:log10(nptp.y1) +
                    dens.percomp.y1:log10(tssx.y1),
                  #log10.Hectares1MPlusSAVComp,
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE_SEM), 
  log10(ptn.y1) %~~% log10(fertilizern_kg.y1),
  #log10(ptp.y1) %~~% log10(ptn.y1),
  #log10(nptp.y1) %~~% log10(nptn.y1),
  log10(nptn.y1) %~~% log10(manuren_kg.y1),
  log10(ptp.y1) %~~% log10(fertilizerp_kg.y1),
  log10(fertilizerp_kg.y1) %~~% log10(fertilizern_kg.y1),
  data = RmSubE_SEM)

summary(RmSubELOGy1.sem)



#####leffy mods####
#lefcharnk Subes SEMs 
# Create SEMs for each salinity zone with TN
TFOH.TN.sem <- psem(
  
  # Predictors of manure and fertilizer
  manure <- lme(log10.manureN ~ Agro + log10.WatershedHa,
                random = ~ 1 | SUBEST_ID,
                correlation = corARMA(form = ~ 1 | SUBEST_ID/Year, q = 2),
                sav_wsmodelTFOH),
  
  
  fertilizer <- lme(log10.fertilizerN ~ Agro + log10.WatershedHa,
                    random = ~ 1 | SUBEST_ID,
                    correlation = corARMA(form = ~ 1 | SUBEST_ID/Year, q = 2),
                    sav_wsmodelTFOH),
  
  # Predictors of nutrients
  nptp <- lme(log10.nptn ~ Developed + log10.manureN + log10.fertilizerN + log10.WatershedHa,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID/Year, q = 2),
              sav_wsmodelTFOH),
  
  ptp <- lme(log10.ptn ~ Developed + log10.WatershedHa,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID/Year, q = 2),
             sav_wsmodelTFOH),
  
  # Predictors of SAV
  SAV <- lme(log10.Area_HA ~
               log10.Area_HA1 +
               log10.flow *
               log10.nptn +
               log10.ptn +
               log10.tssx +
               log10.Richness +
               log10.Hectares1MPlusSAVComp,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID/Year, q = 2),
             sav_wsmodelTFOH),
  
  sav_wsmodelTFOH
  
)

##
MH.TN.sem <- psem(
  
  # Predictors of manure and fertilizer
  manure <- lme(log10(manuren_kg.adj) ~ Agro + log10(WatershedHa),
                random = ~ 1 | SUBEST_ID,
                correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                RmSubE988_SEM),
  
  
  fertilizer <- lme(log10(fertilizern_kg.adj) ~ Agro + log10(WatershedHa),
                    random = ~ 1 | SUBEST_ID,
                    correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                    RmSubE988_SEM),
  
  # Predictors of nutrients
  nptp <- lme(log10(nptn.adj) ~ Developed + log10(manuren_kg.adj)+ log10(fertilizern_kg.adj) + log10(WatershedHa),
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              RmSubE988_SEM),
  
  ptn <- lme(log10(ptn.adj) ~ Developed + log10(WatershedHa),
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             RmSubE988_SEM),
  
  # Predictors of SAV
  SAV <- lme(dens.percomp.change ~
               dens.percomp.y1 +
               log10(flow) *
               log10(nptn.adj) +
               log10(ptn.adj) +
               log10(tssx.adj), #+
               #log10.Richness +
               #log10.Hectares1MPlusSAVComp,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
             RmSubE988_SEM),
  
  RmSubE988_SEM
  
)
summary(MH.TN.sem)

