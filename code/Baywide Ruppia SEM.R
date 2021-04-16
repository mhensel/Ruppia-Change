
#Baywide Ruppia SEMs
library(tidyverse); library(readxl); library(patchwork);library(beyonce)
library(lme4); library(MuMIn); library(DHARMa); library(piecewiseSEM); library(nlme); library(semPlot)



RmZone_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env.csv")
RmZone_Env <- read.csv("../data/RmZone_Env.csv")
RmZone_Env.tss <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env.tss.csv")
RmZone_spmeEnv.tss <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_spmeEnv.tss.csv")



##Generate the Rm_SEM dataset here####
#drop NAs, filter out nonsense values

Rm_SEM <- RmZone_Env %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0) #%>% filter(year > 1999)
Rm_SEM.tss <- RmZone_Env.tss  %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0)
#USE ME!! below
Rm_SEMspme.tss <- RmZone_spmeEnv.tss  %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0) %>% filter(TSSr.spme < 200)

colSums(is.na(Rm_SEMspme.tss))

Rm_SEM85 <- RmZone_Env8515 %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0) #%>% filter(year > 1999)
Rm_SEM.tss85 <- RmZone_Env8515.tss  %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0)
Rm_SEMspme8515.tss <- RmZone_spmeEnv8515.tss  %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0)

Rm_SEMspme9010.tss <- RmZone_spmeEnv9010.tss  %>% drop_na() %>% filter(Sal.spme > 0) %>% filter(Secc.spme > 0) %>% filter(ChlA.spme > 0)

write_csv(Rm_SEM, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/Rm_SEM.csv")
write_csv(Rm_SEM, "./data/Rm_SEM.csv")
write_csv(Rm_SEM.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/Rm_SEM.tss.csv")
write_csv(Rm_SEM.tss, "./data/Rm_SEM.tss.csv")

write_csv(Rm_SEMspme.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/Rm_SEMspme.tss.csv")
write_csv(Rm_SEMspme.tss, "./data/Rm_SEMspme.tss.csv")

Rm_SEM  #1455 obs
Rm_SEM.tss #1368 obs 

###sems below, listed
RmBay.sem  #p = .21, spring means no logging, TN included makes it worse but we are going to keep it in
RmBay.sem1 #p = .34, same as above but TN out. 
RmBayLOG.sem # p = .53, everything logged in RmBay.sem
RmBayLMLOG.sem # p = .38, same as RmBayLOG.sem but linear model instead of ME
RmBayLM.sem #p = .8, same as RmBay.sem but linear

RmBayLOG_03.sem
RmBayLOG.tss.sem #p = .3, RmBayLOG with TSS (uses Rm_SEM.tss)


qplot(y = log10(TSSr.spme), x = log10(ChlA.spme), color = dens.percomp.change, data = Rm_SEMspme.tss) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "right")


AIC(RmBayLOG.sem, RmBayLOGsp.sem)


Rm_SEMspme.tss <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/Rm_SEMspme.tss.csv")
Rm_SEM <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/Rm_SEM.csv")

##p = .157, C = 11.8
#N, Sal have direct effect on RuChange
#include chla~sal and its .4, c = 4.0
RmBayLOGsp.tss.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TSSr.spme) +
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEMspme.tss),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(TSSr.spme) +
                  log10(ChlA.spme) +
                 # log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEMspme.tss),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
             #    log10(TSSr.spme) +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
           #      (dens.percomp.y1:log10(TSSr.spme)) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION , q = 1),
               control = lmeControl(opt = "optim"),
               data = Rm_SEMspme.tss),
  log10(TSSr.spme) %~~% log10(Sal.spme),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(TN.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  log10(TN.spme) %~~% log10(TSSr.spme),
  log10(TP.spme) %~~% log10(TSSr.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
#  log10(ChlA.spme) %~~% log10(Sal.spme),
  data = Rm_SEMspme.tss)

summary(RmBayLOGsp.tss.sem)
coefs(RmBayLOGsp.tss.sem)
dSep(RmBayLOGsp.tss.sem)
fisherC(RmBayLOGsp.tss.sem)

#TSS removed, slightly different DF used (RM_SEM)

#which DF to use acutally. better fit with spme.tss and TSS in but not as a predictor of change. also Rm_SEM fits p = .05, no tss but using tss dataset p = .07

#ChlA ~ Sal appears to be a big mover of the fit (include and its very good)

RmBayLOGsp.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
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
               correlation = corARMA(form = ~ 1 | STATION , q = 1),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(TN.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  data = Rm_SEM)

summary(RmBayLOGsp.sem)

##BEST SEM LOG no TSS##
#Log everytnig SEM. Highest P value here
RmBayLOG.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  #log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 1),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  #log10(TP.spme) %~~% log10(Sal.spme),
  #log10(TN.spme) %~~% log10(Sal.spme),
  data = Rm_SEM)

summary(RmBayLOG.sem)
coefs(RmBayLOG.sem)
dSep(RmBayLOG.sem)
fisherC(RmBayLOG.sem)


#no secci fits fine p = .2 c = 5.9
RmBayLOG.chla.sem <- psem(
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TSSr.spme)
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 1),
                control = lmeControl(opt = "optim"),
                data = Rm_SEMspme.tss),
  # Seccsp <- lme(log10(Secc.spme) ~
  #                log10(ChlA.spme) +
  #log10(Sal.spme) +
  #                log10(TN.spme) +
  #                log10(TP.spme),
  #              random = ~ 1 | STATION,
  #              correlation = corARMA(form = ~ 1 | STATION, q = 1),
  #             control = lmeControl(opt = "optim"),
  #              data = Rm_SEM),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TSSr.spme) +
                 #          log10(TP.spme) +
                 log10(TN.spme) + 
                          log10(TSSr.spme) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 #        (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) ,
               #    (log10(Secc.spme):dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 1),
               control = lmeControl(opt = "optim"),
               data = Rm_SEMspme.tss),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  #log10(TP.spme) %~~% log10(Sal.spme),
  #log10(TN.spme) %~~% log10(Sal.spme),
  data = Rm_SEMspme.tss)

#year as fixed effect not station, fits a little better actually
RmBayLOGsp.tss.semYR <- psem(
  TSSsp <- lme(log10(TSSr.spme) ~
                 log10(TP.spme) +
                 log10(TN.spme),
               random = ~ 1 | year,
               correlation = corARMA(form = ~ 1 | year, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEMspme.tss),
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TSSr.spme) +
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | year,
                correlation = corARMA(form = ~ 1 | year, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEMspme.tss),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(TSSr.spme) +
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | year,
                correlation = corARMA(form = ~ 1 | year, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEMspme.tss),
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
               random = ~ 1 | year,
               correlation = corARMA(form = ~ 1 | year, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEMspme.tss),
  log10(TSSr.spme) %~~% log10(Sal.spme),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Secc.spme),
  log10(Secc.spme) %~~% log10(Sal.spme),
  data = Rm_SEMspme.tss)

summary(RmBayLOGsp.tss.semYR)

AIC(RmBayLOGsp.tss.semYR, RmBayLOGsp.tss.sem)



#annual means (not spring) doesnt fit
#Chla.me*y1 signif, no sal
RmBayAnn.sem <- psem(
  TSS <- lme(log10(TSSr.me) ~
                 log10(TP.me) +
                 log10(TN.me),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss),
  ChlA <- lme(log10(ChlA.me) ~
                  log10(TSSr.me) +
                  log10(TP.me) +
                  log10(TN.me),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  Secc <- lme(log10(Secc.me) ~
                  log10(TSSr.me) +
                  # log10(ChlA.me) +
                  log10(Sal.me) +
                  log10(TN.me) +
                  log10(TP.me),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  RuInt <- lme(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(TSSr.me) +
                 log10(Sal.me) + 
                 log10(ChlA.me) + 
                 log10(TP.me) +
                 log10(TN.me) + 
                 log10(Secc.me) +
                 (dens.percomp.y1:log10(TSSr.me)) +
                 (dens.percomp.y1:log10(Secc.me)) +
                 (dens.percomp.y1:log10(Sal.me)) + 
                 (dens.percomp.y1:log10(ChlA.me)) + 
                 (log10(TP.me):dens.percomp.y1) +
                 (log10(TN.me):dens.percomp.y1) +
                 (log10(Secc.me):dens.percomp.y1),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss),
  log10(TSSr.me) %~~% log10(Sal.me),
  log10(TN.me) %~~% log10(TP.me),
  log10(ChlA.me) %~~% log10(Sal.me),
  log10(ChlA.me) %~~% log10(Secc.me),
  #log10(Secc.spme) %~~% log10(TSSr.spme),
  data = Rm_SEM.tss)

summary(RmBayAnn.sem)

#y1 means (not spring) DOES NOT FIT
RmBayy1.sem <- psem(
  ChlAy1me <- lme(ChlA.y1me ~
                  TP.y1me +
                  TN.y1me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
  Seccy1me <- lme(Secc.y1me ~
                  TSSr.y1me +
                  ChlA.y1me +
                  Sal.y1me +
                  TN.y1me +
                  TP.y1me,
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss),
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

#plot(RmBayLOG.sem, return = FALSE, node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white", width = .7),
#     edge_attrs = data.frame(style = "solid", color = "black", arrowsize = 1, penwidth = 1),
#     ns_dashed = T, alpha = 0.05, show = "std", #digits = 2,
#     add_edge_label_spaces = T)

##NoTssTN removed p = .31 no logs
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

#linear logged using lms fits fine at .363 but not for TSS data
RmBayLMLOG.sem <- psem(
  TSS <- lm(log10(TSSr.me) ~
               log10(TP.me) +
               log10(TN.me), 
            data = Rm_SEM.tss),
  ChlAsp <- lm(log10(ChlA.spme) ~
                  log10(TP.spme) +
                  log10(TN.spme) + 
                  log10(TSSr.spme),
                data = Rm_SEM.tss),
  Seccsp <- lm(log10(Secc.spme) ~
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme) +
                  log10(TSSr.spme),
                data = Rm_SEM.tss),
  RuInt <- lm(dens.percomp.change ~
                 dens.percomp.y1 +
                 log10(Sal.spme) + 
                 log10(ChlA.spme) + 
                 log10(TP.spme) +
                 log10(TN.spme) + 
                 log10(Secc.spme) +
                log10(TSSr.spme) +
                 (dens.percomp.y1:log10(Secc.spme)) +
                 (dens.percomp.y1:log10(Sal.spme)) + 
                 (dens.percomp.y1:log10(ChlA.spme)) + 
                 (log10(TP.spme):dens.percomp.y1) +
                 (log10(TN.spme):dens.percomp.y1) +
                 (log10(Secc.spme):dens.percomp.y1) +
                (log10(TSSr.spme):dens.percomp.y1),
               data = Rm_SEM.tss),
  log10(TN.spme) %~~% log10(TP.spme),
  log10(Secc.spme) %~~% log10(ChlA.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  data = Rm_SEM.tss)

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

###tests
#this one dont fit as well 
RmBayLOGsp.tss8515.sem <- psem(
  TSSsp <- lme(log10(TSSr.spme) ~
                 log10(TP.spme) +
                 log10(TN.spme),
               random = ~ 1 | STATION,
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss85),
  ChlAsp <- lme(log10(ChlA.spme) ~
                  log10(TSSr.spme) +
                  log10(TP.spme) +
                  log10(TN.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss85),
  Seccsp <- lme(log10(Secc.spme) ~
                  log10(TSSr.spme) +
                  log10(ChlA.spme) +
                  log10(Sal.spme) +
                  log10(TN.spme) +
                  log10(TP.spme),
                random = ~ 1 | STATION,
                correlation = corARMA(form = ~ 1 | STATION, q = 2),
                control = lmeControl(opt = "optim"),
                data = Rm_SEM.tss85),
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
               correlation = corARMA(form = ~ 1 | STATION, q = 2),
               control = lmeControl(opt = "optim"),
               data = Rm_SEM.tss85),
  log10(TSSr.spme) %~~% log10(Sal.spme),
  log10(TN.spme) %~~% log10(TP.spme),
  #log10(Secc.spme) %~~% log10(Sal.spme),
  log10(ChlA.spme) %~~% log10(Sal.spme),
  log10(TP.spme) %~~% log10(Sal.spme),
  data = Rm_SEM.tss85)

summary(RmBayLOGsp.tss8515.sem)



### BAY-WIDE STRUCTURAL EQUATION MODEL larfchunk #########

# Create SEM
station.sem <- psem(
  # Predictors of chl-a
  chla <- lme(log10.chla ~
                tp +
                tn +
                l  og10.WTEMP,
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
RmSubE_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE_Env.csv")
RmSubE982_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE_Env982.csv")


RmSubE_SEM <- RmSubE_Env %>% drop_na() %>% filter(tssx.adj < 10000000)
RmSubE988_SEM <- RmSubE_Env982 %>% drop_na() ##USE THIS ONE

%>% filter(ptn.adj > 0) %>% filter(tssx.adj > 0) %>% filter(tssx.adj < 10000000,) %>% filter(ptp.adj > 0) %>% filter(ptn.adj > 0) %>% filter(ptn.y1 > 0) #%>% filter(year > 1999)
#filtering the ptn takes out some data that fucked the SEM so idk whats up here. 

write_csv(RmSubE988_SEM, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE988_SEM.csv")
write_csv(RmSubE988_SEM, "./data/RmSubE988_SEM.csv")



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
RmSubE988_SEM <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE988_SEM.csv")

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

#just N fits p = .84, .97 on untrimmed. 
#note- manure and fertilizer breaks all these models
RmSubELOGtrim.N.sem <-  psem(
  tss <- lme(log10(tssx.adj) ~ flow + Agro,
             random = ~ 1 | SUBEST_ID,
             correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
             data = RmSubE_SEM),
  
  nptn <- lme(log10(nptn.adj) ~  flow + Agro ,
              random = ~ 1 | SUBEST_ID,
              correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
              data = RmSubE_SEM),
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
                  data = RmSubE_SEM), 
  log10(nptn.adj) %~~% log10(tssx.adj), 
  data = RmSubE_SEM)
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
#do
RmSubELOGtrimflowx.sem <-  psem(
 # tss <- lme(log10(tssx.adj) ~ flow + Agro,
 #            random = ~ 1 | SUBEST_ID,
  #           correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2), 
  #           data = RmSubE988_SEM),
  #
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
     #               log10(tssx.adj) +
                    #dens.percomp.y1:flow +
                    dens.percomp.y1:log10(nptp.adj)+
                    dens.percomp.y1:log10(nptn.adj)+
     #               dens.percomp.y1:log10(tssx.adj),
                  random = ~ 1 | SUBEST_ID,
                  correlation = corARMA(form = ~ 1 | SUBEST_ID, q = 2),
                  data = RmSubE988_SEM), 
  log10(fertilizern_kg.adj) %~~% log10(fertilizerp_kg.adj),
  log10(manuren_kg.adj) %~~% log10(manurep_kg.adj),
  #log10(fertilizern_kg.adj) %~~% Agro,
  #log10(manuren_kg.adj) %~~% Agro,
  #Agro %~~% log10(fertilizerp_kg.adj),
  #Agro %~~% log10(manurep_kg.adj),
  #log10(nptn.adj) %~~% log10(nptp.adj),
  #log10(nptn.adj) %~~% log10(tssx.adj), 
  #log10(nptp.adj) %~~% log10(tssx.adj),
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

