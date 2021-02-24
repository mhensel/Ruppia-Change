#Random Forest Models for Baywide Ruppia Stations. Corresponds to Baywide_RuChange_models.R

library(tidyverse); library(randomForest); library(pdp)
##Below segments: run them if you havent run Baywide_RuChange_models. if not, skip to the Random Forest Sections.####
#master RuppiaStation Zone data over time, created in Baywide_RuChange_models file
RmZone_Env<- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env.csv")

#trimmed the tails, so that stations close to their max or min dont skew the data
RmZone_Env8515 <- RmZone_Env %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
write_csv(RmZone_Env8515, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env8515.csv")

RmZone_Env9010 <- RmZone_Env %>%
  filter(!dens.percomp.y1 > 0.90) %>% filter(!dens.percomp.y1 < 0.10)
#only stations in years that have increases. .tail is better to use
RmZone_Env_Inc <- RmZone_Env %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmZone_Env_Inc.tail <- RmZone_Env %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)

#only stations in years that have decreases .tail is better to use
RmZone_Env_Dec <- RmZone_Env %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmZone_Env_Dec.tail <- RmZone_Env %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)

##here are all of the babies
RmZone_Env_ALL #all RM Zone and all Env data, no filters
RmZone_Env_ALL8515 #all Rm Zone and all Env data, tails filtered
RmZone_Env #all RM Zone and filtered USE Env data, no tails filtered
RmZone_Env8515 #all RM Zone and filtered USE Env data, tails filtered
RmZone_Env9010 #all RM Zone and filtered USE Env data, tails filtered 90 10
RmZone_Env_Inc #increases in RM Zone and filtered USE Env data, tails filtered
RmZone_Env_Inc.tail #inc in RM Zone and filtered USE Env data, top tail filtered 
RmZone_Env_Dec #decreases in RM Zone and filtered USE Env data, tails filtered
RmZone_Env_Dec.tail #dec in RM Zone and filtered USE Env data, bottom tail filtered

#Subest
Subestuary_WQ<- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/subestEnv15.csv")

RmSubE_Env <- Subestuary_WQ %>%
  filter(SUBEST_ID %in% RuppiaSubestuaries$SUBEST_ID) %>%
  full_join(RmZoneSubestuaries)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
#is.na(RmSubE_Env) <- RmSubE_Env == "NaN"
RmSubE_Env[is.nan(RmSubE_Env)] <- 0
is.na(RmSubE_Env) <- RmSubE_Env == "Inf"
is.na(RmSubE_Env) <- RmSubE_Env == "-Inf"
RmSubE_Env <- as.data.frame(RmSubE_Env)

#trimmed the tails, so that Subes close to their max or min dont skew the data
RmSubE_Env8515 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)

RmSubE_Env9010 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.90) %>% filter(!dens.percomp.y1 < 0.10)
#only stations in years that have increases. .tail is better to use
RmSubE_Env_Inc <- RmSubE_Env %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmSubE_Env_Inc.tail <- RmSubE_Env %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)
#only stations in years that have decreases .tail is better to use
RmSubE_Env_Dec <- RmSubE_Env %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmSubE_Env_Dec.tail <- RmSubE_Env %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)

##here are all of the SubE babies
RmSubE_Env #all RM Zone and filtered USE Env data, no tails filtered
RmSubE_Env8515 #all RM Zone and filtered USE Env data, tails filtered
RmSubE_Env9010 #all RM Zone and filtered USE Env data, tails filtered 90 10
RmSubE_Env_Inc #increases in RM Zone and filtered USE Env data, tails filtered
RmSubE_Env_Inc.tail #inc in RM Zone and filtered USE Env data, top tail filtered 
RmSubE_Env_Dec #decreases in RM Zone and filtered USE Env data, tails filtered
RmSubE_Env_Dec.tail #dec in RM Zone and filtered USE Env data, bottom tail filtered


#####Random Forest Baywide Ruppia and VIPs#####
#make sure you take out correlated variables

##Ru STATION Random Forests####

#per.max change unfiltered#
max.compdenschangeRF <- randomForest(dens.percomp.change ~ ., RmZone_Env %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max))  #41%
max.compareachangeRF <- randomForest(SAVArea.percomp.change ~ ., RmZone_Env %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max)) 
varImpPlot(max.compdenschangeRF, type = 2) #solid, Sal.y1Dpos and Sal.me, ChlA.me too
varImpPlot(max.compareachangeRF, type = 2)

#per.composite tails filtered, LIKELY MASTER!
max.compdenschangeRF72 <- randomForest(dens.percomp.change ~ ., RmZone_Env8515 %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #27%
max.compareachangeRF72 <- randomForest(SAVArea.percomp.change ~ ., RmZone_Env8515 %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max)) 
varImpPlot(max.compdenschangeRF72, type = 2) #solid, Sal.me, ChlA.me too. Sal.spmin, Sal.min, Salspme
varImpPlot(max.compareachangeRF72, type = 2)

#per.composite tails filtered, LIKELY MASTER! WITH TSS
max.compdenschangeRF72.tss <- randomForest(dens.percomp.change ~ ., RmZone_Env8515.tss %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #27%
max.compareachangeRF72.tss <- randomForest(SAVArea.percomp.change ~ ., RmZone_Env8515.tss %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max)) 
varImpPlot(max.compdenschangeRF72.tss, type = 2) #solid, Sal.me, ChlA.me too. Sal.spmin, Sal.min, Salspme
varImpPlot(max.compareachangeRF72.tss, type = 2)

###lets examine what factors influences INCREASES in RM
max.dens_InctailRF <- randomForest(dens.percomp.change ~ ., RmZone_Env_Inc.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #18%
max.area_InctailRF <- randomForest(SAVArea.percomp.change ~ ., RmZone_Env_Inc.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max))

varImpPlot(max.dens_InctailRF, type = 2) 
varImpPlot(max.area_InctailRF, type = 2)

###lets examine what factors influences DECREASES in RM
max.dens_DectailRF <- randomForest(dens.percomp.change ~ ., RmZone_Env_Dec.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max))
max.area_DectailRF <- randomForest(SAVArea.percomp.change ~ ., RmZone_Env_Dec.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max))

varImpPlot(max.dens_DectailRF, type = 2) #39% but driven by dens/area
varImpPlot(max.area_DectailRF, type = 2)

####PDPs for important RFs
#par(mar=c(1,1,1,1))
#print out top 15 partial dependence plots to R drive/Ruppia/MH_Ruppia
##likely master ones here
importanceOrder=order(-max.compdenschangeRF72$importance)
names=rownames(max.compdenschangeRF72$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.compdenschangeRF72, pred.data = (RmZone_Env8515 %>% drop_na()), eval(name), main=name, xlab=name)

importanceOrder=order(-max.compareachangeRF72$importance)
names=rownames(max.compareachangeRF72$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.compareachangeRF72, pred.data = as.data.frame(RmZone_Env8515 %>% drop_na()), eval(name), main=name, xlab=name)

#increases...not very inf in the RF
importanceOrder=order(-max.dens_InctailRF$importance)
names=rownames(max.dens_InctailRF$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.dens_InctailRF, pred.data = as.data.frame(RmZone_Env_Inc.tail %>% drop_na()), eval(name), main=name, xlab=name)
#decreases, ok on fit but highly dependent on area/density
importanceOrder=order(-max.dens_DectailRF$importance)
names=rownames(max.dens_DectailRF$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.dens_DecRF, pred.data = as.data.frame(RmZone_Env_Dec.tail %>% drop_na()), eval(name), main=name, xlab=name)

##biplots, another way to visualize
library(pdp)
#par(mar=c(1,2,1,2))

SalmeChlameBi <- partial(max.compdenschangeRF72, pred.var = c("Sal.me", "ChlA.me"), chull = TRUE)
plotPartial(SalmeChlameBi)

SalSalDmeBi <- partial(max.compdenschangeRF72, pred.var = c("Sal.me", "Sal.spme"), chull = TRUE)
plotPartial(SalSalDmeBi)

NPBi <- partial(max.compdenschangeRF72, pred.var = c("TP.me", "TN.me"), chull = TRUE)
plotPartial(NPBi)
#increase/decrease data
SprPCHL <- partial(max.dens_DectailRF, pred.var = c("TP.spme", "ChlA.spme"), chull = TRUE)
plotPartial(SprPCHL)

#
##Ru Subestuaries Random Forests####
#
library(leaps)
RmSub_leap <- RmSubE_Env %>% drop_na()


#per.max change unfiltered#
max.compdenschangeRF_Sub <- randomForest(dens.percomp.change ~ ., RmSubE_Env %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max))  #30%
max.compareachangeRF_Sub <- randomForest(SAVArea.percomp.change ~ ., RmSubE_Env %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max)) 
varImpPlot(max.compdenschangeRF_Sub, type = 2) 
varImpPlot(max.compareachangeRF_Sub, type = 2)

#per.composite tails filtered
max.compdenschangeRF72_Sub <- randomForest(dens.percomp.change ~ ., RmSubE_Env8515 %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #12%
max.compareachangeRF72_Sub <- randomForest(SAVArea.percomp.change ~ ., RmSubE_Env8515 %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max)) 
varImpPlot(max.compdenschangeRF72_Sub, type = 2)
varImpPlot(max.compareachangeRF72_Sub, type = 2)

###lets examine what factors influences INCREASES in RM
max.dens_InctailRF_Sub <- randomForest(dens.percomp.change ~ ., RmSubE_Env_Inc.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #20
max.area_InctailRF_Sub <- randomForest(SAVArea.percomp.change ~ ., RmSubE_Env_Inc.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max))

varImpPlot(max.dens_InctailRF_Sub, type = 2) 
varImpPlot(max.area_InctailRF_Sub, type = 2)

###lets examine what factors influences DECREASES in RM
max.dens_DectailRF_Sub <- randomForest(dens.percomp.change ~ ., RmSubE_Env_Dec.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp,  -SAVArea.percomp.change, -dens.percomp, -SAVArea.percomp.y1, -denscomp.max)) #56
max.area_DectailRF_Sub <- randomForest(SAVArea.percomp.change ~ ., RmSubE_Env_Dec.tail %>% drop_na() %>% select(-dens.change, -dens.weight.mean, -SAVArea.change, -SAVArea.prop.change, -SAVArea, -dens.prop.change, -SAVArea.percomp, -dens.percomp.change, -dens.percomp, -dens.percomp.y1, -denscomp.max))

varImpPlot(max.dens_DectailRF_Sub, type = 2) #39% but driven by dens/area
varImpPlot(max.area_DectailRF_Sub, type = 2)

####PDPs for important RFs
#par(mar=c(1,1,1,1))
#print out top 15 partial dependence plots to R drive/Ruppia/MH_Ruppia
##likely master ones here
importanceOrder=order(-max.compdenschangeRF72_Sub$importance)
names=rownames(max.compdenschangeRF72_Sub$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.compdenschangeRF72_Sub, pred.data = (RmSubE_Env8515 %>% drop_na()), eval(name), main=name, xlab=name)

importanceOrder=order(-max.compareachangeRF72_Sub$importance)
names=rownames(max.compareachangeRF72_Sub$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.compareachangeRF72_Sub, pred.data = as.data.frame(RmSubE_Env8515 %>% drop_na()), eval(name), main=name, xlab=name)

#increases
importanceOrder=order(-max.dens_InctailRF_Sub$importance)
names=rownames(max.dens_InctailRF_Sub$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.dens_InctailRF_Sub, pred.data = as.data.frame(RmSubE_Env_Inc.tail %>% drop_na()), eval(name), main=name, xlab=name)
#decreases, 
importanceOrder=order(-max.dens_DectailRF_Sub$importance)
names=rownames(max.dens_DectailRF_Sub$importance)[importanceOrder][1:15]
par(mfrow=c(5, 3), xpd=NA)
for (name in names)
  partialPlot(max.dens_DecRF_Sub, pred.data = as.data.frame(RmSubE_Env_Dec.tail %>% drop_na()), eval(name), main=name, xlab=name)

##biplots, another way to visualize
library(pdp)
#par(mar=c(1,2,1,2))

SalmeChlameBi <- partial(max.compdenschangeRF72, pred.var = c("Sal.me", "ChlA.me"), chull = TRUE)
plotPartial(SalmeChlameBi)

SalSalDmeBi <- partial(max.compdenschangeRF72, pred.var = c("Sal.me", "Sal.spme"), chull = TRUE)
plotPartial(SalSalDmeBi)

NPBi <- partial(max.compdenschangeRF72, pred.var = c("TP.me", "TN.me"), chull = TRUE)
plotPartial(NPBi)
#increase/decrease data
SprPCHL <- partial(max.dens_DectailRF, pred.var = c("TP.spme", "ChlA.spme"), chull = TRUE)
plotPartial(SprPCHL)

