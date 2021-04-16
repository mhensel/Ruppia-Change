##Calculate Baywide Ruppia change variables (part 1), merge with Env Data from CBP (part 2), graph up trends (part 3), do some mixed models (part 4)
#First half of file is Station Zones, second half is Subestuaries data
#this file creates datasets used in the SEM and in the Random Forests. 
#RmZone_Env, RmSubE_Env are the merged Env/Ruppia cover files. 

library(tidyverse); library(readxl); library(patchwork);library(beyonce)
library(randomForest); library(leaps); library(lme4)
library(MuMIn); library(DHARMa); library(piecewiseSEM); library(nlme)

#CAN SKIP THIS first 200 lines to step 3, IF YOU RUN merge Ruppia Env.R after tidyCBPWQ data 2019.R

#
##
###
####STATION ZONES#####
###
##
#

############Part 1: Filtering RM_Zone Stations/SubEs, calculating Ru cover and change variables#######

#this is all the stations, with % of it in the Rm_Zone and a composite max HA variable which I will use to create a maximum density weighted area composite (denscomp.max). as of now 11/11/20, its just area * .85 but dave is going to get a more accurate number
RuppiaOverlap_StationZone <- read_excel("/Volumes/savshare2/Current Projects/Ruppia/Ruppia areas in Chesapeake Bay/Ruppia SAV Zones Overlap with Station Zones.xlsx")

#Use this to filter out the station zones that arent in the RM_Zone and calculate the maximum composite area in any given year (denscomp.max)
RuppiaStations <- RuppiaOverlap_StationZone %>% 
  dplyr::filter(RMZoneSAV_HA > 0) %>%
  mutate(denscomp.max = (RMZoneD4_Ha*.85) + (RMZoneD3_Ha*.55) + (RMZoneD2_Ha*.25) + (RMZoneD1_Ha*.05)) %>%
  select(STATION, denscomp.max, RMZoneSAV_HA)

#Ruppia coverage in each Station zone per year is "SAVAreaHa". with density class too
SAVyear_StationZone <- read_excel("/Volumes/savshare2/Current Projects/Ruppia/Ruppia areas in Chesapeake Bay/Ruppia SAV Area by Year by Station Zone.xlsx")

#calculate Ruppia coverage in each station zone for each year and year - 1
#Ruppia response variables created here: 
#dens.weight.mean = total density-weighted-mean area in zone, SAVArea = total area of SAV in zone, dens.change = change in dens weight mean area from y1, SAVArea.change,  dens.prop.change = scaled density weighted area change, SAVArea.prop.change = scaled area change
#dens.permax.change density weighted means but scaled from the maximum extent of the zone. UPDATE 11/10: now dens.percomp.change and related dens.percomp are dens weighted means, scaled from the maxiumum composite extent in each Station zone. 

####Station Ruppia change over time####
RmZoneStations <- SAVyear_StationZone %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  #filter(between(Year, 1990, 2019)) %>% use this if you want to create a certain time dataset
  mutate(per.cov = case_when(Density == 1 ~ .05,
                             Density == 2 ~ .25, 
                             Density == 3 ~ .55, 
                             Density == 4 ~ .85)) %>% #convert to density weighted means using these Braun-Bloc metrics
  mutate(dens_cov = SAVAreaHa * per.cov) %>%
  group_by(STATION, Year) %>%
  summarize(dens.weight.mean = sum(dens_cov), SAVArea = sum(SAVAreaHa)) %>%
  mutate(dens.weight.mean.y1 = lag(dens.weight.mean, order_by = Year, k = 1), SAVArea.y1 = lag(SAVArea, order_by = Year, k = 1)) %>%
  mutate(dens.change = dens.weight.mean - dens.weight.mean.y1, SAVArea.change = SAVArea - SAVArea.y1) %>%
  mutate(dens.prop.change = (dens.weight.mean - dens.weight.mean.y1)/dens.weight.mean.y1, SAVArea.prop.change = (SAVArea - SAVArea.y1)/SAVArea.y1) %>% ungroup() %>% 
  full_join(RuppiaStations) %>% group_by(STATION) %>% 
  mutate(#dens.permax = dens.weight.mean/max(dens.weight.mean), dens.permax.y1 = dens.weight.mean.y1/max(dens.weight.mean), dens.permax.change = (dens.weight.mean-dens.weight.mean.y1)/max(dens.weight.mean), 
         dens.percomp = dens.weight.mean/denscomp.max, dens.percomp.y1 = dens.weight.mean.y1/denscomp.max, dens.percomp.change = (dens.weight.mean-dens.weight.mean.y1)/denscomp.max,
         #SAVArea.permax = SAVArea/max(SAVArea), SAVArea.permax.y1 = SAVArea.y1/max(SAVArea), SAVArea.permax.change = (SAVArea-SAVArea.y1)/max(SAVArea), 
         SAVArea.percomp = SAVArea/RMZoneSAV_HA, SAVArea.percomp.y1 = SAVArea.y1/RMZoneSAV_HA, SAVArea.percomp.change = (SAVArea-SAVArea.y1)/RMZoneSAV_HA) %>% 
  rename("year" = "Year") %>% select(-RMZoneSAV_HA)

#if you want to rewrite, go for it.
write_csv(RmZoneStations, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZoneStations.csv")
write_csv(RmZoneStations, "../data/RmZoneStations.csv")


############Part 2: Merging Station Zone with Environmental Data#######
#load in the CBP data, updated to 2019 as of October 2020. this file comes from Ruppia/Data/CBP WQ Station data/tidyCBPWQ data 2019.R





#NOTE: this is a little TOO much data, so the one below is going to be used instead of this. But can load it anyways
Env_Var_ALL.CBP_WQ<- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_ALL.CBP_WQ.csv")
#filter out the bad stations in the ENV data 
##DF of baywide Ruppia and All ENV#
RmZone_Env_ALL <- Env_Var_ALL.CBP_WQ %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) %>%
  #filter(dens.weight.mean > 10) %>% filter(dens.weight.mean.y1 > 10) %>% use this to filter if needed
  select(-CBSEG_2003)
#is.na(RmZone_Env_ALL) <- RmZone_Env_ALL == "NaN"
RmZone_Env_ALL[is.nan(RmZone_Env_ALL)] <- 0
is.na(RmZone_Env_ALL) <- RmZone_Env_ALL == "Inf"
is.na(RmZone_Env_ALL) <- RmZone_Env_ALL == "-Inf"
RmZone_Env_ALL <- as.data.frame(RmZone_Env_ALL)
write_csv(RmZone_Env_ALL, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env_ALL.csv")

###Merge RuppiaStation coverage data with CBP WQ Station data####
Env_Var_nosumy1sp.CBP_WQ <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_nosumy1sp.CBP_WQ.csv") #more informative and not as unweildy. Has y1 means, summer y1, spring this year. (eliminates y1 spring, summer this year). 

#here is the TSS data. 
Env_Var.tss <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var.tss.csv")

####MASTER ENV AND BAYWIDE DATASET RmZone_Env ####
#NO TSS
RmZone_Env <- Env_Var_nosumy1sp.CBP_WQ %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) #%>%
  #filter(denscomp.max > 5) #this gets rid of any zones that are really tiny. about 200 data points less. originally i thought this would clean things up but it actually weakens the SEM? idk about for these graphs

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
#is.na(RmZone_Env) <- RmZone_Env == "NaN"
RmZone_Env[is.nan(RmZone_Env)] <- 0
is.na(RmZone_Env) <- RmZone_Env == "Inf"
is.na(RmZone_Env) <- RmZone_Env == "-Inf"
RmZone_Env <- as.data.frame(RmZone_Env)

write_csv(RmZone_Env, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env.csv")

#if u wanna look at just the data in the segments. this has y and y1 for both summer and spring variables 
RmZone_Env_yearly <- Env_Var_yearly %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) 
RmZone_Env_spring <- Env_Var_spring %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations)
RmZone_Env_grow <- Env_Var_grow %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) 

####Our USE dataframes are these babies below. 11/3/20#####
RmZone_Env_ALL8515 <- RmZone_Env_ALL %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#this is here in case needed
#RmZone_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env.csv")

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

###Repeat Step 2 with TSSr data####
#here is the TSS only data
Env_Var.tss <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var.tss.csv")
#TSSr is RAW. TSS was FINAL_VALUE which is just nothing before 99 so its not in this dataset
#
Env_Var_nosumy1sp.CBP_WQ.tss <- full_join(Env_Var_nosumy1sp.CBP_WQ, Env_Var.tss)


####MASTER ENV AND BAYWIDE DATASET with TSS RmZone_Env.tss ####
RmZone_EnvUNFILT <- Env_Var_nosumy1sp.CBP_WQ.tss %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) %>%
filter(denscomp.max > 5) #this gets rid of any zones that are really tiny. about 200 data points less. originally i thought this would clean things up but it actually weakens the SEM? idk about for these graphs

#is.nan.data.frame <- function(x)
#  do.call(cbind, lapply(x, is.nan))
#is.na(RmZone_Env) <- RmZone_Env == "NaN"
RmZone_EnvUNFILT[is.nan(RmZone_EnvUNFILT)] <- 0
is.na(RmZone_EnvUNFILT) <- RmZone_EnvUNFILT == "Inf"
is.na(RmZone_EnvUNFILT) <- RmZone_EnvUNFILT == "-Inf"
RmZone_EnvUNFILT <- as.data.frame(RmZone_EnvUNFILT)

##going to go ahead and do the station filters now: LE3.6, LE3.7,CB5.4W,CB7.1,CB7.1N,CB7.1S,EE3.4, EE3.5 are the stations you cant use any data TSSr before 1999. 
#CB4.1E, CB5.1, CB5.2,CB5.3,LE2.3,EE1.1,EE2.1,EE2.2,EE3.0,EE3.1,EE3.2,EE3.3,ET4.2,ET5.2,ET8.1,ET9.1,LE2.2,RET2.4,LE3.2,LE3.3, LE3.4, CB7.3E are the stations where it is ok to use the TSSr data.
#need to: filter everything before 1999 out on LE3.6, LE3.7,CB5.4W,CB7.1,CB7.1N, CB7.1S,EE3.4, EE3.5

RmZone_Env1999 <- RmZone_EnvUNFILT %>%
  filter(STATION %in% c("LE3.6", "LE3.7","CB5.4W","CB7.1","CB7.1N", "CB7.1S","EE3.4", "EE3.5")) %>% filter(year > 1999) 

RmZone_EnvFILT <- RmZone_EnvUNFILT %>% #group_by(STATION, year) %>%
  filter(!STATION %in% c("LE3.6", "LE3.7","CB5.4W","CB7.1","CB7.1N", "CB7.1S","EE3.4", "EE3.5"))

RmZone_Env.tss <- bind_rows(RmZone_Env1999, RmZone_EnvFILT)

write_csv(RmZone_Env.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env.tss.csv")

####Our USE dataframes are these babies below. 1/15/21#####
#trimmed the tails, so that stations close to their max or min dont skew the data
RmZone_Env8515.tss <- RmZone_Env.tss %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
write_csv(RmZone_Env8515.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZone_Env8515.tss.csv")

RmZone_Env9010.tss <- RmZone_Env.tss %>%
  filter(!dens.percomp.y1 > 0.90) %>% filter(!dens.percomp.y1 < 0.10)
#only stations in years that have increases. .tail is better to use
RmZone_Env_Inc.tss <- RmZone_Env.tss %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmZone_Env_Inc.tail.tss <- RmZone_Env.tss %>%
  filter(dens.percomp.change > 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)

#only stations in years that have decreases .tail is better to use
RmZone_Env_Dec.tss <- RmZone_Env.tss %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
RmZone_Env_Dec.tail.tss <- RmZone_Env.tss %>%
  filter(dens.percomp.change < 0) %>%
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.percomp.y1 < 0.15)

##here are all of the TSS babies
RmZone_Env.tss #all RM Zone and filtered USE Env data, no tails filtered
RmZone_Env8515.tss #all RM Zone and filtered USE Env data, tails filtered
RmZone_Env9010.tss #all RM Zone and filtered USE Env data, tails filtered 90 10
RmZone_Env_Inc.tss #increases in RM Zone and filtered USE Env data, tails filtered
RmZone_Env_Inc.tail.tss #inc in RM Zone and filtered USE Env data, top tail filtered 
RmZone_Env_Dec.tss #decreases in RM Zone and filtered USE Env data, tails filtered
RmZone_Env_Dec.tail.tss #dec in RM Zone and filtered USE Env data, bottom tail filtered





##
###
####
####################START HERE NOW! 3/1/21#########
####
###
##

RmZone_spmeEnv.tss <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_spmeEnv.tss.csv")

#Original Workflow: take these babies and go over to Baywide Ruppia RFs.R, and run some random forests. Then come back to the below Part 3

######Part 3: Correlations and Exploratory Station Zone Figures#######
Rudenstime <- qplot(x = year, y = dens.weight.mean, data = RmZone_Env8515.tss) + 
  #stat_summary(fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.5) +
  ylab("Density weighted area/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")

RudensSTATIONtime <-qplot(x = year, y = dens.weight.mean, data = RmZone_Env8515.tss) + 
  #stat_summary(fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(aes(color = STATION), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.5) +
  ylab("Density weighted area/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right")


ggplot(data = RmZone_Env8515.tss) + 
  stat_summary(aes(x = dens.percomp, y = dens.percomp.change, color = dens.percomp.y1), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6) +
  ylab("Density weighted area mean change/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right")

ggplot(data = RmSubE_Env982.tss) + 
  stat_summary(aes(x = dens.percomp, y = dens.percomp.change, color = dens.percomp.y1), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6) +
  ylab("Density weighted area mean change/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right")

#density change over time
meanchange <- qplot(x = year, y = dens.change, color = STATION, data = RmZone_Env8515.tss) + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  ylab("Density weighted area mean change/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")

propchange <- qplot(x = year, y = dens.prop.change, color = STATION, data = RmZone_Env8515.tss) + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  ylab("Prop dens weight area change/year") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")

####mean maxcomp change over time####
meancompchange <- 
  ggplot(data = Rm_SEM) + 
  stat_summary(aes(x = year, y = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "black") +
  stat_summary(aes(x = year, y = dens.percomp.change), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = 1.5, color = "black") +
  #geom_smooth(aes(x = year, y = dens.percomp.change, group = STATION, color = STATION), method = "lm", alpha = 0.5, size = .5) +
  geom_hline(yintercept = 0, color = "red") +
  ylim(-.3, .3) +
  theme_bw(base_size=20) + 
  ylab("Ruppia change") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "")


densweightmeTime <- 
ggplot(data = RmZone_Env9010.tss) + 
  #geom_point(aes(x = year, y = dens.percomp.change), color = "black", alpha = 0.5, size = 1.5) +
  stat_summary(aes(x = year, y = dens.weight.mean), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "black") +
  stat_summary(aes(x = year, y = dens.weight.mean), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = 1.5, color = "black") +
  # geom_hline(yintercept = 0, color = "red") +
  theme_bw(base_size=20) + 
  ylab("Ruppia DWM (HA/Zone)") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")



propchange + meanchange + meancompchange
densweightmeTime
#connowingo dam discharge data, comes from tidyCBPWQ data 2019.R file at the end
meancompchange / conno_monthly
meancompchange / conno_Spmonthly
meancompchange / conno_Wimonthly


#####correlatoons and stuff##
library(corrplot)
cor(RmZone_Env8515, method = "pearson", use = "complete.obs")
cor <- rcorr(as.matrix(RmZone_Env8515 %>% drop_na))
cor(RmZone_Env8515$Sal.me, RmZone_Env8515$Sal.Dneg)

##correlation plots ####
qplot(y = Secc.sumy1me, x = ChlA.y1me, color = year, data = RmZone_Env8515.tss%>% drop_na() )

qplot(y = Temp.y1me, x = log(ChlA.me), color = STATION, data = RmZone_Env8515%>% drop_na() ) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "")

qplot(y = Sal.y1ran, x = Sal.sumy1ran, color = year, data = RmZone_Env8515%>% drop_na() ) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "")

qplot(y = log10(TSSr.spme), x = log10(ChlA.spme), color = year, data = Rm_SEMspme.tss%>% drop_na() ) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "")

qplot(y = TP.spme, x = TP.spmax, color = STATION, data = RmZone_Env8515%>% drop_na() ) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "")

qplot(y = log10(TN.spme), x = log10(TP.spme), color = dens.percomp.change, data = Rm_SEM%>% drop_na() ) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "")

ggplot(data = RmZone_Env8515 %>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = Sal.me, y = ChlA.me)) +
  geom_point(aes(x = Sal.me, y = ChlA.me, color = STATION)) +
  theme_bw(base_size=14) + theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")

#######Scatterplots of possible predictors and dens.percomp.change ######
qplot(y = Temp.me, x = ChlA.me, data = RmZone_Env_ALL%>% drop_na() )
a<-qplot(y = dens.percomp.change, x = TN.me, color = year, data = RmZone_Env8515)+ theme(legend.position="none")
b<-qplot(y = dens.percomp.change, x = TP.sumy1me, color = year, data = RmZone_Env8515)
a+b

qplot(y = dens.permax.change, x = Sal.y1Dpos, color = STATION, data = RmZone_Env_ALL%>% drop_na() )
##CHL A plots####
#ChlA.me (log or not log) ChlAme is pretty good. ChlA sp ran, spDme also solid predictors
#permax looks better than percomp. but can switch if needed 
logchladens <-
ggplot(data = RmZone_Env8515.tss ) + 
  geom_smooth(method = "lm", aes(x = log10(ChlA.me), y = dens.percomp.change)) +
  geom_point(aes(x = log10(ChlA.me), y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

logchlaspmedens <- 
  ggplot(data = Rm_SEMspme8515.tss ) + 
  geom_smooth(method = "lm", aes(x = log10(ChlA.spme), y = dens.percomp.change)) +
  geom_point(aes(x = log10(ChlA.spme), y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
chladens <-
  ggplot(data = RmZone_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = ChlA.me, y = dens.percomp.change)) +
  geom_point(aes(x = ChlA.me, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
chlaspmedens <- 
  ggplot(data = RmZone_spmeEnv8515.tss ) + 
  geom_smooth(method = "lm", aes(x = log10(ChlA.spme), y = dens.percomp.change), color = "black") +
  geom_point(aes(x = log10(ChlA.spme), y = dens.percomp.change), color = "green") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
chlaDnegdens <- 
ggplot(data = RmZone_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = ChlA.Dneg, y = dens.percomp.change)) +
  geom_point(aes(x = ChlA.Dneg, y = dens.percomp.change, color = year)) +
  xlim(-100, 10) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

logchladens + logchlaspmedens
chladens + chlaspmedens
chlaDnegdens





###SALINITY plots####
Salmedens <- 
ggplot(data = RmZone_Env8515) +# %>% filter(year > 1990)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = Sal.me, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.me, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "right")
#Sal.sp me
Salspmedens <- 
  ggplot(data = RmZone_spmeEnv8515.tss) + # %>% filter(year > 1991)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = Sal.spme, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.spme, y = dens.percomp.change), color = "blue") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "")
Salmindens <- 
  ggplot(data = RmZone_Env8515 %>% filter(year > 1991)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = Sal.min, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.min, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "right")

salspDnegdens<-
  ggplot(data =RmZone_Env8515.tss) + 
  #geom_smooth(method = "lm", aes(x = Sal.spDneg, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.spDneg, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
salspDposdens<-
ggplot(data =RmZone_Env8515.tss) + 
  #geom_smooth(method = "lm", aes(x = Sal.spDpos, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.spDpos, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

salspDposdens + salspDnegdens
meancompchange / salDtime
salDtime <- 
ggplot(data =RmZone_Env8515.tss) + 
  #geom_smooth(method = "lm", aes(x = Sal.spDpos, y = dens.percomp.change)) +
  stat_summary(fun.data = mean_se, geom = "line", aes(y = Sal.spDpos, x = year)) +
  stat_summary(fun.data = mean_se, geom = "line", aes(y = Sal.spDneg, x = year), color = "red") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

#Sal.Dneg predicts drops in Ru. Dme not as good
salDnegdensDec <-
ggplot(data = RmZone_Env_Dec) + 
  geom_smooth(method = "lm", aes(x = Sal.Dneg, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.Dneg, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

salmindensDec <-
  ggplot(data = RmZone_Env_Dec) + 
  geom_smooth(method = "lm", aes(x = Sal.min, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.min, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

salspDmedensDec <-
  ggplot(data = RmZone_Env_Dec) + 
  geom_smooth(method = "lm", aes(x = Sal.spDneg, y = dens.percomp.change)) +
  geom_point(aes(x = Sal.spDneg, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

##TP plots####
TPspmedens2K <- 
  ggplot(data = RmZone_Env8515 %>% filter(year > 2000)) + 
  geom_smooth(method = "lm", aes(x = log10(TP.spme), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.spme), y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")
TPspmedens <- 
  ggplot(data = RmZone_Env8515) + 
  geom_smooth(method = "lm", aes(x = log10(TP.spme), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.spme), y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

ggplot(data = RmZone_Env8515) + 
  geom_smooth(method = "lm", aes(x = log10(TP.spme), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.spme), y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TPy1summedens <- 
  ggplot(data = RmZone_Env8515) + 
  geom_smooth(method = "lm", aes(x = TP.sumy1me, y = dens.percomp.change)) +
  geom_point(aes(x = TP.sumy1me, y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")
TPmedens <- 
  ggplot(data = RmZone_Env8515) + 
  geom_smooth(method = "lm", aes(x = log10(TP.me), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.me), y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")
TPmaxdens <- 
  ggplot(data = RmZone_Env8515) + 
  geom_smooth(method = "lm", aes(x = log10(TP.max), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.max), y = dens.percomp.change, color = year)) + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TP.dposdensInc <- 
ggplot(data = RmZone_Env_Inc %>% filter(year > 2000)) + 
  geom_smooth(method = "lm", aes(x = log10(TP.Dpos), y = dens.percomp.change)) +
  geom_point(aes(x = log10(TP.Dpos), y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

#TN plots####
TNmedens2K <- 
ggplot(data = RmZone_Env9010.tss %>% filter(year > 2000)) + 
  #geom_smooth(method = "lm", aes(x = TN.me, y = dens.permax.change)) +
  geom_point(aes(x = TN.me, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TNspmedens <- 
  ggplot(data = RmZone_Env8515.tss) + 
  #geom_smooth(method = "lm", aes(x = TN.me, y = dens.permax.change)) +
  geom_point(aes(x = TN.spme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TNspDmedens <- #notrends
  ggplot(data = RmZone_Env8515.tss) + 
  #geom_smooth(method = "lm", aes(x = TN.me, y = dens.permax.change)) +
  geom_point(aes(x = TN.spDme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

#TEMP plots#### 
#not much here w Temp at all. 
Tempmedens <- 
ggplot(data = RmZone_Env8515) + 
  #geom_smooth(method = "lm", aes(x = Temp.spDme, y = dens.percomp.change)) +
  geom_point(aes(x = Temp.me, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) + xlim(10, 25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
Tempy1DposInc <-  ##this one was supposedly corr w increases but looks like nah
  ggplot(data = RmZone_Env_Inc) + 
  #geom_smooth(method = "lm", aes(x = Temp.spDme, y = dens.percomp.change)) +
  geom_point(aes(x = Temp.y1Dpos, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

#Secc plots####
seccspmedens <- 
  ggplot(data = RmZone_spmeEnv9010.tss) + 
  #geom_smooth(method = "lm", aes(x = Secc.y1Dme, y = dens.percomp.change)) +
  geom_point(aes(x = Secc.spme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")
#nothing really here
seccy1Dmedens2k <- 
ggplot(data = RmZone_Env8515 %>% filter(year >1999)) + 
  #geom_smooth(method = "lm", aes(x = Secc.y1Dme, y = dens.percomp.change)) +
  geom_point(aes(x = Secc.y1Dme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")
seccDposdens <-  #.Dpos, sumy1Dme supposedly corr w decreases
  ggplot(data = RmZone_Env8515) + 
  # geom_smooth(method = "lm", aes(x = Secc.spme, y = dens.percomp.change)) +
  geom_point(aes(x = Secc.sumy1Dme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "")

####TSS plots####
RmZone_Env.tss
TSSmedens <- 
  ggplot(data = RmZone_Env8515.tss) +# %>% filter(year > 1990)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = TSSr.me, y = dens.percomp.change)) +
  geom_point(aes(x = TSSr.me, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "right")
#TSS.sp me
TSSspmedens <- 
  ggplot(data = RmZone_spmeEnv8515.tss) +# %>% filter(year > 1990)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = log10(TSSr.spme), y = dens.percomp.change), color = "brown") +
  geom_point(aes(x = log10(TSSr.spme), y = dens.percomp.change), color = "brown") +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "right")

TSSspmedens99 <- 
  ggplot(data = RmZone_Env8515.tss %>% filter(year > 1999)) + #%>% drop_na() ) + 
  geom_smooth(method = "lm", aes(x = TSSr.spme, y = dens.percomp.change)) +
  geom_point(aes(x = TSSr.spme, y = dens.percomp.change, color = year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = "right")



#########ENV VAR PLOTS ###############

(Salmedens + Salspmedens) 
(logchladens + logchlaspmedens) 
chladens + chlaspmedens 
TPspmedens2K 
TPspmedens + TPmedens
TSSspmedens
Tempmedens + TNmedens2K + seccy1Dmedens2k


##ENV vars over time plots####
meancompchange
densweightmeTime

RuSalChla_time <- meancompchange / 
salmetime /chlametime 
RuTPTSS_time <- meancompchange / 
TPmetime / TSSsptime

EnvRu_time <- RuSalChla_time | RuTPTSS_time


salmetimeVSBay <- 
ggplot(data = Env_Var_nosumy1sp.CBP_WQ) + 
  geom_point(aes(x = year, y = Sal.spme, color = dens.percomp.change)) +
  stat_summary(aes(x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(aes(x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(data = RmZone_Env8515,aes(x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "blue") +
  stat_summary(data = RmZone_Env8515,aes( x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .6, color = "blue") +
  scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = 14, color = "red") +
  ylab("Salinity (mean PPT/year)") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

salmetime <- 
  ggplot(data = RmZone_Env8515) + 
  #geom_point(aes(x = year, y = Sal.spme), alpha = .5, size = 1, color = "blue") +
  stat_summary(data = RmZone_Env8515,aes(x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "blue") +
  stat_summary(data = RmZone_Env8515,aes( x = year, y = Sal.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .9, color = "blue") +
  scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 14, color = "red") +
  ylab("Salinity") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

chlametime <- 
  ggplot(data = RmZone_Env8515) + 
  #geom_point(aes(x = year, y = ChlA.spme), alpha = .5, size = 1, color = "green") +
  stat_summary(data = RmZone_Env8515,aes(x = year, y = ChlA.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "green") +
  stat_summary(data = RmZone_Env8515,aes( x = year, y = ChlA.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .9, color = "green") +
  scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 14, color = "red") +
  ylab("ChlA") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TPmetime <- 
  ggplot(data = RmZone_Env8515) + 
  #geom_point(aes(x = year, y = ChlA.spme), alpha = .5, size = 1, color = "green") +
  stat_summary(data = RmZone_Env8515,aes(x = year, y = TP.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "purple") +
  stat_summary(data = RmZone_Env8515,aes( x = year, y = TP.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .9, color = "purple") +
  scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 14, color = "red") +
  ylab("TP ug/L") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

TNtime <-
  ggplot(data = RmZone_Env8515) + 
  #geom_point(aes(x = year, y = TN.me, color = dens.percomp.change)) +
  stat_summary(aes(x = year, y = TN.me, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(aes(x = year, y = TN.me, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .6, color = "black") +
  scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = .7, color = "red") +
  ylab("TN (mean ug/L/year)") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

Temptime <-
  ggplot(data = RmZone_Env8515) + 
 # geom_point(aes(x = year, y = Temp.me, color = dens.percomp.change)) +
  stat_summary(aes(x = year, y = Temp.spme), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(aes(x = year, y = Temp.spme), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .6, color = "black") +
  #scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 11, color = "red") +
  #ylim(10,25) +
  ylab("Temperature (mean C/year") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

Secctime <-
  ggplot(data = RmZone_Env8515) + 
  #geom_point(aes(x = year, y = Secc.me, color = dens.percomp.change)) +
  stat_summary(aes(x = year, y = Secc.me, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .6, color = "black") +
  stat_summary(aes(x = year, y = Secc.me, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .6, color = "black") +
  scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 11, color = "red") +
  ylab("Secci depth (mean m/year)") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

TSSsptime <-
  ggplot(data = RmZone_Env8515.tss) + 
  #geom_point(aes(x = year, y = ChlA.spme), alpha = .5, size = 1, color = "green") +
  stat_summary(data = RmZone_Env8515.tss,aes(x = year, y = TSSr.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "brown") +
  stat_summary(data = RmZone_Env8515.tss,aes( x = year, y = TSSr.spme, color = dens.percomp.change), fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = .9, color = "brown") +
  scale_color_gradient(low="blue", high="red") +
  #geom_hline(yintercept = 14, color = "red") +
  ylab("TSS") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme_bw(base_size=20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

######ENV change over time graphs####
meancompchange4box + TPmetime + chlametime + salmetime 
meancompchange4box + Temptime +  Secctime + chlametime + TNtime + salmetime + TPmetime

multiWQtimeLOG <- 
  ggplot(data = RmZone_Env8515.tss) + 
  stat_summary(aes(x = year, y = log10(ChlA.spme)), color = "green", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = year, y = log10(Sal.spme)), color = "grey", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = year, y = log10(TSSr.spme)), color = "brown", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  #stat_summary(aes(x = year, y = dens.percomp), color = "black", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 2) +
  #geom_hline(yintercept = 14, color = "red") +
  xlab("") + 
  scale_y_continuous(name = "log Spring ChlA/TSS/Salinity")+ #,sec.axis = sec_axis(~., name="Prop. Ruppia Cover")) +
  scale_x_continuous(n.breaks = 20) +
  theme_bw(base_size=24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

coeff <- 1
multiWQtime <- 
  ggplot(data = RmZone_Env8515.tss) + 
  #stat_summary(aes(x = year, y = dens.weight.mean), color = "black", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 2) +
  stat_summary(aes(x = year, y = ChlA.spme), color = "green", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = year, y = Sal.spme), color = "brown", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  #stat_summary(aes(x = year, y = TSSr.spme), color = "brown", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  scale_y_continuous(
    # Features of the first axis
    name = "Spring ChlA",
    # Add a second axis and specify its features *coeff
    sec.axis = sec_axis(~., name="Spring Salinity")) + 
  theme(
    axis.title.y = element_text(color = "green"),
    axis.title.y.right = element_text(color = "brown")) + 
  xlab("") + 
  scale_x_continuous(n.breaks = 10) +
  theme_bw(base_size=24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")


#####ENV Vars and Density change scatterplots to keep HERE####
meanmaxchange
logchladens
Salmedens
salDnegdensDec
(Salmedens + Salspmedens) 
(logchladens + logchlaspmedens) 
chladens + chlaspmedens 
chlametime + salmetime 
TPspmedens2K + TPy1summedens
TNmedens2K + TNspmedens
Tempmedens
seccy1Dmedens2k

logchlaspmedens + Salspmedens
TPy1summedens + TNmedens2K
#patch2b <- (bur2r + bur2a) / (muss2)/ (lit2) / (mel2r + mel2a) 
#patch2b + plot_annotation(tag_levels = 'A')

######Interaction plots#####
chlaXdmw <- (RmZone_spmeEnv.tss %>% mutate(ChlAxDMW = log10(ChlA.spme)*dens.percomp.y1) %>%
  qplot(data = ., x = ChlAxDMW, y = dens.percomp.change, color = dens.percomp.y1) +
    geom_point(size = 2.5) +scale_color_gradient(low="blue", high="red") +
    xlab("Chla spring * Ru y-1") + ylab("Ruppia change") +
    theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right"))

a <- (RmZone_spmeEnv.tss %>% mutate(SalxDMW = log10(Sal.spme)*dens.percomp.y1) %>%
    qplot(data = ., x = SalxDMW, y = dens.percomp.change, color = Sal.spme))

a + salXdmw

salXdmw <- (RmZone_spmeEnv.tss %>% mutate(SalxDMW = log10(Sal.spme)*dens.percomp.y1) %>%
  qplot(data = ., x = SalxDMW, y = dens.percomp.change, color = dens.percomp.y1)+scale_color_gradient(low="blue", high="red") + 
    geom_point(size = 2.5) +scale_color_gradient(low="blue", high="red") +
    xlab("Salinity spring * Ru y-1") + ylab("Ruppia change") +theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))

hist(Rm_SEMspme.tss$TSSr.spme)

b + TSSXdmw

b <- (Rm_SEMspme.tss %>% mutate(TSSxDMW = log10(TSSr.spme)*dens.percomp.y1) %>%
        qplot(data = ., x = TSSxDMW, y = dens.percomp.change, color = log10(TSSr.spme)))

TSSXdmw <- (Rm_SEMspme.tss %>% mutate(TSSxDMW = log10(TSSr.spme)*dens.percomp.y1) %>%
              qplot(data = ., x = TSSxDMW, y = dens.percomp.change, color = dens.percomp.y1)+scale_color_gradient(low="blue", high="red") + 
              xlab("TSS spring * Ru y-1") + ylab("Ruppia change") + theme_bw(base_size=20)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right"))

tpXdmw <- (RmZone_spmeEnv.tss %>% mutate(TPxDMW = log10(TP.spme)*dens.percomp.y1) %>%
              qplot(data = ., x = TPxDMW, y = dens.percomp.change, color = dens.percomp.y1)+ scale_color_gradient(low="blue", high="red") + 
             xlab("TP spring * Ru y-1") + ylab("Ruppia change") + theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))

tnXdmw <- (RmZone_spmeEnv.tss %>% mutate(TNxDMW = log(TN.sumy1me)*dens.percomp.y1) %>%
             qplot(data = ., x = TNxDMW, y = dens.percomp.change, color = dens.percomp.y1)+ scale_color_gradient(low="blue", high="red") + 
             xlab("TN spring * Ru y-1") + ylab("Ruppia change") + theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))
secXdmw <- (RmZone_Env8515 %>% mutate(SecxDMW = Secc.y1Dme*dens.percomp.y1) %>%
             qplot(data = ., x = SecxDMW, y = dens.percomp.change, color = dens.percomp.y1)+ scale_color_gradient(low="blue", high="red") + 
              xlab("Secc spring * Ru y-1") + ylab("Ruppia change") + theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))

SalDmeXdmw <- (RmZone_Env8515 %>% mutate(SalDmexDMW = Sal.y1Dpos*dens.percomp.y1) %>%
              qplot(data = ., x = SalDmexDMW, y = dens.percomp.change, color = dens.percomp.y1)+ scale_color_gradient(low="blue", high="red") + 
              xlab("Saly1Dpos spring * Ru y-1") + ylab("Ruppia change") + theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))

chlaXdmw + salXdmw
tpXdmw + tnXdmw + secXdmw


####Part 4: Modelling Ruppia Station Zone Change Over time#######
#
c<-qplot(y = dens.percomp.change, x = log10(TSSr.spme), color = year, data = RmZone_Env8515.tss )
d<-qplot(y = dens.percomp.change, x = TSSr.spme, color = year, data = RmZone_Env8515.tss)
c+d
#note- dens.percomp.change is predicted better than permax. also RmZone_Env is the dataset to use
#no environmental R2 = .29
rmmeanNOenv.lmer <- lmer(dens.percomp.change ~ dens.weight.mean.y1 +SAVArea.percomp.y1 +SAVArea.y1 + dens.percomp.y1 + (1|STATION),  data = RmZone_Env)
##testing lefcheck baywide model
rmmean.jon.lmer1 <- lmer(dens.percomp.change ~ dens.weight.mean.y1 + log10(ChlA.me) + log10(Secc.me) + log10(Temp.me) + log10(Sal.me) + TN.me + TP.me + (1|STATION),  data = RmZone_Env)

#### BEST MODEL W INTERACTION TERMS (AIC)
#mix means and spring r2 .39
rmmeanenv.lmer <- lmer(dens.percomp.change ~ 
                         SAVArea.y1 + ChlA.me + Sal.me + Secc.me + 
                         (dens.percomp.y1:Sal.spme)+ 
                         (dens.percomp.y1:Sal.y1Dpos) + 
                         (dens.percomp.y1:ChlA.spme) + 
                         (TP.spme:dens.percomp.y1)+
                         (TP.sumy1me:dens.percomp.y1) +
                         (Secc.y1Dme:dens.percomp.y1) + 
                         (1|STATION) + (1|year),  data = RmZone_Env) 
#log rmmean r2 .40
rmmeanenv.loglmer <- lmer(dens.percomp.change ~ 
                            SAVArea.y1 + log10(ChlA.me) + Sal.me + log10(Secc.me) + 
                            (dens.percomp.y1:Sal.spme)+ 
                            (dens.percomp.y1:Sal.y1Dpos) + 
                            (dens.percomp.y1:log10(ChlA.spme)) + 
                            (log10(TP.spme):dens.percomp.y1)+
                            (log10(TP.sumy1me):dens.percomp.y1) +
                            (Secc.y1Dme:dens.percomp.y1) + 
                            (1|STATION),  data = RmZone_Env)
summary(rmmeanenv.loglmer)


ranef(rmmeanenv.lmer)
as.data.frame(VarCorr(rmmeanenv.lmer))
mgcv::extract.lme.cov(rmmeanenv.lmer,start.level=1)

#spring only r2 40
rmmeanenv.splmer <- lmer(dens.percomp.change ~ SAVArea.y1 + ChlA.spme + Sal.spme + Secc.spme + (dens.percomp.y1:Sal.spme)+ (dens.percomp.y1:Sal.y1Dpos) + (dens.percomp.y1:ChlA.spme) + (log10(TP.spme):dens.percomp.y1)+(TP.sumy1me:dens.percomp.y1) +(Secc.y1Dme:dens.percomp.y1) + (1|STATION),  data = RmZone_Env)


AICc(rmmeanenv.lmer, rmmeanenvtest.lmer1, rmmeanenv.splmer, rmmeanenv.loglmer)
r.squaredGLMM(rmmeanenv.lmer)
r.squaredGLMM(rmmeanenv.loglmer)
r.squaredGLMM(rmmeanenvtest.lmer1)
r.squaredGLMM(rmmeanenv.splmer)
###chris leaps model
chrismodel <- lmer(dens.percomp.change ~ ChlA.me:dens.percomp.y1 + Sal.spmin:dens.percomp.y1 + Sal.y1Dpos:dens.percomp.y1 + SAVArea.y1 + Secc.y1Dme:dens.percomp.y1 + TP.y1max:dens.percomp.y1 + (1|STATION), RmZone_Env)

AICc(rmmeanNOenv.lmer,rmmeanenv.lmer, rmmeanenvtest.lmer1, chrismodel, rmmean.jon.lmer1, predicttest.lmer)

Anova(rmmeanenv.lmer, test.statistic = "F")
Anova(rmmeanenvtest.lmer1, test.statistic = "F")
anova(rmmeanenv.lmer, rmmeanenvtest.lmer1)

vif(rmmeanenv.lmer)
drop1(rmmeanenv.lmer)

plot(simulateResiduals(rmmeanenv.lmer))
qqnorm(residuals(rmmeanenv.lmer))
testZeroInflation(simulateResiduals(rmmeanenv.lmer))




##predictive check on lmer
y1<-predict(rmmeanenv.lmer)
predy <- y1 * daw$denscomp.max

daw <- RmZone_Env %>% select(year, STATION, denscomp.max, dens.weight.mean, dens.weight.mean.y1, dens.percomp.y1, SAVArea.percomp.y1, SAVArea.y1, ChlA.me, Sal.me, Sal.spme, Sal.y1Dpos, ChlA.spme, TP.y1max, Secc.y1Dme) %>% na.omit()

prd.lm <- lm(dens.weight.mean ~ dens.weight.mean.y1 + predy, data = daw)
summary(prd.lm)

prd.lmer <- lmer(dens.weight.mean ~ dens.weight.mean.y1 + predy + (1|STATION), data = daw)
summary(prd.lmer)
r.squaredGLMM(prd.lmer)

predicttest.lmer <- lmer(dens.percomp.change ~ dens.percomp.y1+ log(ChlA.y1me) + Secc.y1me + Sal.y1me + TN.y1me + Temp.y1me + TP.y1me + log(ChlA.spme) + Secc.spme + Sal.spme + TN.spme + Temp.spme + TP.spme + log(ChlA.sumy1me) + Secc.sumy1me + Sal.sumy1me + TN.sumy1me + Temp.sumy1me + log(TP.sumy1me) + (1|STATION),  data = RmZone_Env)
summary(predicttest.lmer)
r.squaredGLMM(predicttest.lmer)

##
y2<-predict(rmmeanenvTSS.lmer)
predy <- y2 * daw$denscomp.max

daw2 <- RmZone_Env.tss %>% select(year, STATION, denscomp.max, dens.weight.mean, dens.weight.mean.y1, dens.percomp.y1, SAVArea.percomp.y1, SAVArea.y1, ChlA.me, Sal.me, Sal.spme, Sal.y1Dpos, ChlA.spme, TP.y1max, TP.spme, TSSr.spme, Secc.y1Dme) %>% na.omit()

SAVArea.y1, ChlA.me, Sal.me, Secc.me, TSSr.spme +
  (dens.percomp.y1:Sal.spme)+ 
  (dens.percomp.y1:Sal.y1Dpos) + 
  (dens.percomp.y1:ChlA.spme) + 
  (TP.spme:dens.percomp.y1)+
  (TP.sumy1me:dens.percomp.y1) +
  (Secc.y1Dme:dens.percomp.y1) + 
  (TSSr.spme:dens.percomp.y1) 

prd.lm <- lm(dens.weight.mean ~ dens.weight.mean.y1 + predy, data = daw)
summary(prd.lm)

prd.lmer <- lmer(dens.weight.mean ~ dens.weight.mean.y1 + predy + (1|STATION), data = daw)
summary(prd.lmer)
r.squaredGLMM(prd.lmer)

predicttest.lmer <- lmer(dens.percomp.change ~ dens.percomp.y1+ log(ChlA.y1me) + Secc.y1me + Sal.y1me + TN.y1me + Temp.y1me + TP.y1me + log(ChlA.spme) + Secc.spme + Sal.spme + TN.spme + Temp.spme + TP.spme + log(ChlA.sumy1me) + Secc.sumy1me + Sal.sumy1me + TN.sumy1me + Temp.sumy1me + log(TP.sumy1me) + (1|STATION),  data = RmZone_Env)
summary(predicttest.lmer)
r.squaredGLMM(predicttest.lmer)


##TSS Station ME models####
rmmeanenvTSS.lmer <- lmer(dens.percomp.change ~ 
                         SAVArea.y1 + ChlA.me + Sal.me + Secc.me + TSSr.spme +
                         (dens.percomp.y1:Sal.spme)+ 
                         (dens.percomp.y1:Sal.y1Dpos) + 
                         (dens.percomp.y1:ChlA.spme) + 
                         (TP.spme:dens.percomp.y1)+
                         (TP.sumy1me:dens.percomp.y1) +
                         (Secc.y1Dme:dens.percomp.y1) + 
                           (TSSr.spme:dens.percomp.y1) +
                         (1|STATION) + (1|year),  data = RmZone_Env.tss) 

#spring only r2 40
rmmeanenvTSS.splmer <- lmer(dens.percomp.change ~ SAVArea.y1 + ChlA.spme + Sal.spme + Secc.spme + TSSr.spme + (dens.percomp.y1:Sal.spme)+ (dens.percomp.y1:Sal.y1Dpos) + (dens.percomp.y1:ChlA.spme) + (log10(TP.spme):dens.percomp.y1)+(TP.sumy1me:dens.percomp.y1) +(Secc.y1Dme:dens.percomp.y1) +(dens.percomp.y1:Sal.spme) + (1|STATION),  data = RmZone_Env.tss)


AICc(rmmeanenvTSS.lmer,rmmeanenvTSS.splmer, rmmeanenv.lmer)
r.squaredGLMM(rmmeanenvTSS.lmer)
r.squaredGLMM(rmmeanenvTSS.splmer)

#extra mixed effects models
#decreases is modeled very well 
rmmeanenvtest.lmerdec1 <- lmer(dens.percomp.change ~ dens.percomp.y1 + SAVArea.y1+ Sal.Dneg + Sal.spme + log(ChlA.me) + log(ChlA.spme) + TP.sumy1me + (1|STATION) + (1|year),  data = RmZone_Env_Dec.tail)
r.squaredGLMM(rmmeanenvtest.lmerdec1)
#
rmmeanenvtest.lm <- lm(dens.percomp.change ~ dens.weight.mean.y1+ dens.percomp.y1 + (dens.percomp.y1*Sal.spme) + log(ChlA.me) + (dens.percomp.y1*log(ChlA.spme)) + (TP.sumy1me*dens.percomp.y1), data = RmZone_Env)
Anova(rmmeanenvtest.lm, test.statistic = "F")
summary(rmmeanenvtest.lm)

rmmeanenvtest.glm <- glm(dens.percomp.change ~ 
                                 SAVArea.y1 + ChlA.me + Sal.me + Secc.me + 
                                 (dens.percomp.y1:Sal.spme)+ 
                                 (dens.percomp.y1:Sal.y1Dpos) + 
                                 (dens.percomp.y1:ChlA.spme) + 
                                 (log10(TP.spme):dens.percomp.y1)+
                                 (TP.sumy1me:dens.percomp.y1) +
                                 (Secc.y1Dme:dens.percomp.y1),  
                         data = RmZone_Env, family = negative.binomial(theta = 5, link = "log"))
Anova(rmmeanenvtest.glmer)

rmmeanenvtest.glmer <- glmer(dens.percomp.change ~ 
                         SAVArea.y1 + ChlA.me + Sal.me + Secc.me + 
                         (dens.percomp.y1:Sal.spme)+ 
                         (dens.percomp.y1:Sal.y1Dpos) + 
                         (dens.percomp.y1:ChlA.spme) + 
                         (log10(TP.spme):dens.percomp.y1)+
                         (TP.sumy1me:dens.percomp.y1) +
                         (Secc.y1Dme:dens.percomp.y1) + 
                         (1|STATION),  data = RmZone_Env, family = quasipoisson) 
Anova(rmmeanenvtest.glmer)



#
##
###
####SUBESTUARIES#####
###
##
#
#
###PART 1: Filter, Merge Ruppia Subestuary coverage data with Subestuary Env data####
RuppiaOverlap_SubestuaryZone <- read_excel("/Volumes/savshare2/Current Projects/Ruppia/Ruppia areas in Chesapeake Bay/Ruppia SAV Zones Overlap with Subestuaries.xlsx")

#Use this to filter out the station zones that arent in the RM_Zone and calculate the maximum composite area in any given year (denscomp.max)
RuppiaSubestuaries <- RuppiaOverlap_SubestuaryZone %>% 
  dplyr::filter(RMZoneSAV_HA > 0) %>%
  mutate(denscomp.max = (RMZoneD4_Ha*.85) + (RMZoneD3_Ha*.55) + (RMZoneD2_Ha*.25) + (RMZoneD1_Ha*.05)) %>%
  select(SUBEST_ID, denscomp.max, RMZoneSAV_HA)

#Ruppia coverage in each Station zone per year is "SAVAreaHa". with density class too
SAVyear_Subestuary <- read_excel("/Volumes/savshare2/Current Projects/Ruppia/Ruppia areas in Chesapeake Bay/Ruppia SAV Area by Year by Subestuary.xlsx")

#####PArt2 = Subestuary Ruppia Change over time#####
RmZoneSubestuaries <- SAVyear_Subestuary %>%
  filter(SUBEST_ID %in% RuppiaSubestuaries$SUBEST_ID) %>%
  #filter(between(Year, 1990, 2019)) %>% use this if you want to create a certain time dataset
  mutate(per.cov = case_when(Density == 1 ~ .05,
                             Density == 2 ~ .25, 
                             Density == 3 ~ .55, 
                             Density == 4 ~ .85)) %>% #convert to density weighted means 
  mutate(dens_cov = SAVAreaHa * per.cov) %>%
  group_by(SUBEST_ID, Year) %>%
  summarize(dens.weight.mean = sum(dens_cov), SAVArea = sum(SAVAreaHa)) %>%
  mutate(dens.weight.mean.y1 = lag(dens.weight.mean, order_by = Year, k = 1), SAVArea.y1 = lag(SAVArea, order_by = Year, k = 1)) %>%
  mutate(dens.change = dens.weight.mean - dens.weight.mean.y1, SAVArea.change = SAVArea - SAVArea.y1) %>%
  mutate(dens.prop.change = (dens.weight.mean - dens.weight.mean.y1)/dens.weight.mean.y1, SAVArea.prop.change = (SAVArea - SAVArea.y1)/SAVArea.y1) %>% ungroup() %>% 
  full_join(RuppiaSubestuaries) %>% group_by(SUBEST_ID) %>% 
  mutate(#dens.permax = dens.weight.mean/max(dens.weight.mean), dens.permax.y1 = dens.weight.mean.y1/max(dens.weight.mean), dens.permax.change = (dens.weight.mean-dens.weight.mean.y1)/max(dens.weight.mean), 
    dens.percomp = dens.weight.mean/denscomp.max, dens.percomp.y1 = dens.weight.mean.y1/denscomp.max, dens.percomp.change = (dens.weight.mean-dens.weight.mean.y1)/denscomp.max,
    #SAVArea.permax = SAVArea/max(SAVArea), SAVArea.permax.y1 = SAVArea.y1/max(SAVArea), SAVArea.permax.change = (SAVArea-SAVArea.y1)/max(SAVArea), 
    SAVArea.percomp = SAVArea/RMZoneSAV_HA, SAVArea.percomp.y1 = SAVArea.y1/RMZoneSAV_HA, SAVArea.percomp.change = (SAVArea-SAVArea.y1)/RMZoneSAV_HA) %>% 
  select(-RMZoneSAV_HA)

#if you want to rewrite, go for it.
#write_csv(RmZoneSubestuaries, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZoneSubestuaries.csv")
write_csv(RmZoneSubestuaries, "../Ruppia/RmZoneSubestuaries.csv")

Subestuary_WQ<- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/Subestuary_WQ.csv")

####MASTER ENV AND Subestuary DATASET RmSubE_Env ####
RmSubE_Env <- Subestuary_WQ %>%
  filter(SUBEST_ID %in% RuppiaSubestuaries$SUBEST_ID) %>%
  full_join(RmZoneSubestuaries) #%>%
  # filter(denscomp.max > 1)  #use this line to filter out tiny zones. takes out about 300 points
#NOTE: filtering out the tiny zones makes the SEM not fit quite as well. idk

RmSubE_Env[is.nan(RmSubE_Env)] <- 0
is.na(RmSubE_Env) <- RmSubE_Env == "Inf"
is.na(RmSubE_Env) <- RmSubE_Env == "-Inf"
RmSubE_Env <- as.data.frame(RmSubE_Env)

write_csv(RmSubE_Env, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSubE_Env.csv")

#trimmed the tails, so that Subes close to their max or min dont skew the data
RmSubE_Env8515 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
#write_csv(RmSubE_Env8515, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSubE_Env8515.csv")

RmSubE_Env982 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.98) %>% filter(!dens.percomp.y1 < 0.02)
write_csv(RmSubE_Env982, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSubE_Env982.csv")
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
  filter(!dens.percomp.y1 > 0.85) #%>% filter(!dens.permax.y1 < 0.15)

##here are all of the SubE babies
RmSubE_Env #all RM SubEZone and filtered USE Env data, no tails filtered
RmSubE_Env8515 #all RM SubEZone and filtered USE Env data, tails filtered
RmSubE_Env9010 #all RM SubEZone and filtered USE Env data, tails filtered 90 10
RmSubE_Env_Inc #increases in RM SubEZone and filtered USE Env data, tails filtered
RmSubE_Env_Inc.tail #inc in RM SubEZone and filtered USE Env data, top tail filtered 
RmSubE_Env_Dec #decreases in RM SubEZone and filtered USE Env data, tails filtered
RmSubE_Env_Dec.tail #dec in RM SubEZone and filtered USE Env data, bottom tail filtered

#Original Workflow: take these babies and go over to Baywide Ruppia RFs.R, and run some random forests. Then come back to the below Part 3


######Part 3: Correlations and Exploratory Subestuaries #####

RmSubE_Env982 <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE_Env982.csv")

#######Scatterplots of possible predictors and dens.percomp.change ######
e<-qplot(y = log10(tssx.adj), x = log10(nptp.adj), color = Year, data = RmSubE_Env982 %>% drop_na()) + theme(legend.position="none")
f<-qplot(y = dens.percomp.change, x = log10(nptn.adj), color = Year, data = RmSubE_Env8515%>% drop_na()) + theme(legend.position="none")
e+f

##Subes Ruppia plots####
SubERu <- 
  RmZoneSubestuaries %>% group_by(Year) %>% summarise(meanRu = sum(SAVArea)) %>%
  ggplot(data = .) + 
  geom_point(aes(x = Year, y = meanRu), alpha = 0.5, size = 3) +
  stat_summary(aes(x = Year, y = meanRu), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "black") +
  ylab("Subestuary Ruppia Area") +
  theme_bw(base_size=12)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right")

meancompchangeSubE <- 
  ggplot(data = RmSubE988_SEM) + #%>% filter(Year < 2012)) + 
  #geom_point(aes(x = Year, y = dens.percomp.change), color = "black", alpha = 0.5, size = 1.5) +
  stat_summary(aes(x = Year, y = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "purple") +
  stat_summary(aes(x = Year, y = dens.percomp.change), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .9, color = "purple") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw(base_size=20) + 
  ylab("Ruppia change") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

(meancompchangeSubE / meancompchange)
####Statoion amd Subest togethe####
StationSubest_RuChangetime <- 
  ggplot(data = RmZone_Env8515) + 
    stat_summary(data = RmSubE_Env982, aes(x = Year, y = dens.percomp.change), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .9, color = "purple") +
    stat_summary(data = RmSubE_Env982, aes(x = Year, y = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .5, color = "purple") +
  stat_summary(data = RmZone_Env8515, aes(x = year, y = dens.percomp.change), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .9, color = "black") +
    stat_summary(data = RmZone_Env8515, aes(x = year, y = dens.percomp.change), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .5, color = "black") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw(base_size=20) + #scale_color_gradient(low="light green", high="brown") +
  ylab("Ruppia change") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

%>% filter(Year <2015)
StationSubest_Rutime <- 
  ggplot(data = RmZone_Env8515) + 
  stat_summary(data = RmSubE_Env8515 , aes(x = Year, y = dens.percomp), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .9, color = "purple") +
  geom_smooth(method = "lm", se = F, data = RmSubE_Env8515, aes(x = Year, y = dens.percomp), color = "red") +
  #stat_summary(data = RmSubE_Env8515, aes(x = Year, y = dens.percomp), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .5, color = "purple") +
  stat_summary(data = RmZone_Env8515, aes(x = year, y = dens.percomp), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .9, color = "black") +
  #stat_summary(data = RmZone_Env8515, aes(x = year, y = dens.percomp), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .5, color = "black") +
  geom_smooth(method = "lm", se= F, data = RmZone_Env8515, aes(x = year, y = dens.percomp), color = "green") +
  theme_bw(base_size=20) + #scale_color_gradient(low="light green", high="brown") +
  ylab("Ruppia ") + 
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

StationSubest_RutimeALL

###Subest Env over time####
meancompchangeSubE / 
multiWQtimeLOGSubes <- 
  ggplot(data = RmSubE_Env982) + 
  stat_summary(aes(x = Year, y = log10(nptn.adj)), color = "orange", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = Year, y = log10(flow)), color = "blue", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = Year, y = log10(tssx.adj)), color = "brown", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = Year, y = log10(nptp.adj)), color = "red", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  #stat_summary(aes(x = Year, y = dens.percomp), color = "black", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 2) +
  #geom_hline(yintercept = 14, color = "red") +
  xlab("") + 
  scale_y_continuous(name = "log NPTN/Flow/Fert")+ #,sec.axis = sec_axis(~., name="Prop. Ruppia Cover")) +
  scale_x_continuous(n.breaks = 20) +
  theme_bw(base_size=24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

coeff <- 1
multiWQtimeSubes <- 
  ggplot(data = RmSubE_Env982) + 
  stat_summary(aes(x = Year, y = nptn.adj), color = "orange", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = Year, y = flow), color = "blue", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  stat_summary(aes(x = Year, y = fertilizern_kg.adj), color = "brown", fun.data = mean_cl_normal, geom = "line", fun.args = list(mult = 1), size = 1.2) +
  scale_y_continuous(
    # Features of the first axis
    name = "Spring ChlA",
    # Add a second axis and specify its features *coeff
    sec.axis = sec_axis(~., name="Spring Salinity")) + 
  theme(
    axis.title.y = element_text(color = "green"),
    axis.title.y.right = element_text(color = "brown")) + 
  xlab("") + 
  scale_x_continuous(n.breaks = 20) +
  theme_bw(base_size=24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")
  
##load plots####
tssxy1dens <-
  ggplot(data = RmSubE_Env982) + 
  geom_smooth(method = "lm", aes(x = log10(tssx.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(tssx.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

tssxdens <- 
  ggplot(data = RmSubE_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = log10(tssx.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(tssx.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
nptndens <-
  ggplot(data = RmSubE_Env982) + 
  geom_smooth(method = "lm", aes(x = log10(nptn.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(nptn.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
nptny1dens <-
  ggplot(data = RmSubE_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = log10(nptn.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(nptn.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
nptpdens <-
  ggplot(data = RmSubE_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = log10(nptp.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(nptp.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")

tssxy1dens + tssxdens


###Agro/Devo####
devodens <-
  ggplot(data = RmSubE_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = Developed, y = dens.percomp.change)) +
  geom_point(aes(x = Developed, y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")
#Agro
agrodens <-
  ggplot(data = RmSubE_Env8515 ) + 
  geom_smooth(method = "lm", aes(x = Agro, y = dens.percomp.change)) +
  geom_point(aes(x = Agro, y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")


#Manure plots####
manureN <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(manuren_kg.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(manuren_kg.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

manureNy1 <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(manuren_kg.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(manuren_kg.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")


manureP <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(manurep_kg.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(manurep_kg.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

manurePy1 <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(manurep_kg.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(manurep_kg.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

#Fertilizer plots####
fertN <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(fertilizern_kg.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(fertilizern_kg.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

fertNy1 <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(fertilizern_kg.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(fertilizern_kg.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")


fertP <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(fertilizerp_kg.adj), y = dens.percomp.change)) +
  geom_point(aes(x = log10(fertilizerp_kg.adj), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

fertPy1  <- 
  ggplot(data = RmSubE_Env9010) +
  geom_smooth(method = "lm", aes(x = log10(fertilizerp_kg.y1), y = dens.percomp.change)) +
  geom_point(aes(x = log10(fertilizerp_kg.y1), y = dens.percomp.change, color = Year)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right")

####interaction plots SubE####
tssXdmw <- (RmSubE_Env982 %>% mutate(TSSxDMW = log10(tssx.adj)*dens.percomp.y1) %>%
               qplot(data = ., x = TSSxDMW, y = dens.percomp.change, color = dens.percomp.y1) +
               geom_point(size = 2.5) +scale_color_gradient(low="blue", high="red") +
               xlab("TSS * Ru y-1") + ylab("Ruppia change") +
               theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right"))

nptnXdmw <- (RmSubE_Env982 %>% mutate(nptnxDMW = log10(nptn.adj)*dens.percomp.y1) %>%
              qplot(data = ., x = nptnxDMW, y = dens.percomp.change, color = dens.percomp.y1) +
              geom_point(size = 2.5) +scale_color_gradient(low="blue", high="red") +
              xlab("N conc (nonpoint) * Ru y-1") + ylab("Ruppia change") +
              theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))

flowXdmw <- (RmSubE_Env982 %>% mutate(flowxDMW = log10(flow)*dens.percomp.y1) %>%
               qplot(data = ., x = flowxDMW, y = dens.percomp.change, color = dens.percomp.y1) +
               geom_point(size = 2.5) +scale_color_gradient(low="blue", high="red") +
               xlab("Flow* Ru y-1") + ylab("Ruppia change") +
               theme_bw(base_size=25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = ""))


##Part 4: Subestuaries models####
#dens.percomp.y1 explains 25% of the variation
rmSubEenv.lmer <- lmer(dens.percomp.change ~ 
                          SAVArea.y1 + 
                          log10(Developed) +  log10(Agro) + 
                          log10(tssx.y1) + #log10(tssx.adj)
                          log10(fertilizerp_kg.adj) + log10(fertilizern_kg.adj) +
                          log10(nptp.y1) + log10(nptn.y1) +
                          log10(Developed):dens.percomp.y1 + log10(Agro):dens.percomp.y1 + 
                          log10(tssx.y1):dens.percomp.y1 + log10(tssx.adj):dens.percomp.y1 +
                          log10(fertilizerp_kg.adj):dens.percomp.y1 + 
                          log10(fertilizern_kg.adj):dens.percomp.y1 +
                        log10(nptn.y1):dens.percomp.y1 + log10(nptp.y1):dens.percomp.y1 + 
                             (1|SUBEST_ID),  data = RmSubE_Env) 
#best AIC
rmSubEenvY1.lmer <- lmer(dens.percomp.change ~ 
                             Developed.y1 +  Agro.y1 + 
                             log10(tssx.y1) + #log10(tssx.y1) +
                             log10(nptp.y1)+ log10(nptn.y1) +
                             Developed.y1:dens.percomp.y1 + Agro.y1:dens.percomp.y1 + 
                             log10(tssx.y1):dens.percomp.y1 + 
                             log10(nptn.y1):dens.percomp.y1 + 
                           log10(nptp.y1):dens.percomp.y1 + 
                             (1|SUBEST_ID),  data = RmSubE_Env) 

rmSubEenvadj.lmer <- lmer(dens.percomp.change ~ 
                           Developed +  Agro + 
                           log10(tssx.adj) + #log10(tssx.y1) +
                           log10(nptp.adj)+ log10(nptn.adj) +
                           Developed:dens.percomp.y1 + Agro:dens.percomp.y1 + 
                           log10(tssx.adj):dens.percomp.y1 + 
                           log10(nptn.adj):dens.percomp.y1 + 
                           log10(nptp.adj):dens.percomp.y1 + 
                           (1|SUBEST_ID),  data = RmSubE_Env) 
#bestr2
rmSubEenvLeaps.lmer <- lmer(dens.percomp.change ~ SAVArea.y1 +
                              Developed.y1 + #log10(Agro) +
                              log10(tssx.y1) + #log10(tssx.adj)
                              #log10(fertilizerp_kg.adj) + log10(fertilizern_kg.adj) +
                              log10(nptp.y1) + log10(nptn.y1) +
                              ptp.y1 + ptn.y1 +
                              Developed.y1:dens.percomp.y1 + #log10(Agro):dens.percomp.y1 + 
                              log10(tssx.y1):dens.percomp.y1 + #log10(tssx.adj):dens.percomp.y1 +
                              #log10(fertilizerp_kg.adj):dens.percomp.y1 + 
                              #log10(fertilizern_kg.adj):dens.percomp.y1 +
                              log10(nptn.y1):dens.percomp.y1 + log10(nptp.y1):dens.percomp.y1 + 
                              ptn.y1:dens.percomp.y1 + ptp.y1:dens.percomp.y1 + 
                              (1|SUBEST_ID),  data = RmSubE_Env)
rmSubEenvLeapsadj.lmer <- lmer(dens.percomp.change ~ SAVArea.y1 +
                              Developed + #log10(Agro) +
                              tssx.adj + #log10(tssx.adj)
                              #log10(fertilizerp_kg.adj) + log10(fertilizern_kg.adj) +
                              nptp.adj + nptn.adj +
                              ptp.adj + ptn.adj +
                              Developed:dens.percomp.y1 + #log10(Agro):dens.percomp.y1 + 
                              tssx.adj:dens.percomp.y1 + #log10(tssx.adj):dens.percomp.y1 +
                              #log10(fertilizerp_kg.adj):dens.percomp.y1 + 
                              #log10(fertilizern_kg.adj):dens.percomp.y1 +
                              nptn.adj:dens.percomp.y1 + nptp.adj:dens.percomp.y1 + 
                              ptn.adj:dens.percomp.y1 + ptp.adj:dens.percomp.y1 + 
                              (1|SUBEST_ID) + (1|Year),  data = RmSubE_Env)

rmSubE.mini <- lmer(dens.percomp.change ~ SAVArea.y1 +
                      Developed + tssx.adj + nptp.adj + nptn.adj +
                      Developed:dens.percomp.y1 +tssx.adj:dens.percomp.y1 + 
                      nptn.adj:dens.percomp.y1 + nptp.adj:dens.percomp.y1 + 
                      (1|SUBEST_ID) + (1|Year),  data = RmSubE_Env)
nptn.adj + nptp.adj + nptn.y1 + nptp.y1 + ptp.y1

leapersub <- lmer(dens.percomp.change ~ Developed*dens.percomp.y1 + nptn.adj*dens.percomp.y1 + nptp.adj*dens.percomp.y1 + nptn.y1*dens.percomp.y1 + nptp.y1*dens.percomp.y1 + ptp.y1*dens.percomp.y1 + Developed.y1*dens.percomp.y1 +(1|Year) + (1|SUBEST_ID), data = RmSubE_Env)
r.squaredGLMM(leapersub)

r.squaredGLMM(rmSubEenv.lmer)
r.squaredGLMM(rmSubEenvY1.lmer)
r.squaredGLMM(rmSubEenvadj.lmer)
r.squaredGLMM(rmSubEenvLeaps.lmer)
r.squaredGLMM(rmSubE.mini)
car::Anova(rmSubEenv.lmer, test.statistic = "F")
Anova(rmSubEenvLeapsadj.lmer)
AICc(rmSubEenv.lmer, rmSubEenvY1.lmer, rmSubEenvadj.lmer, rmSubEenvLeaps.lmer, rmSubEenvLeapsadj.lmer, rmSubE.mini, leapersub)

vif(rmSubEenv.lmer)
drop1(rmSubEenv.lmer)

plot(simulateResiduals(rmSubEenv.lmer))
qqnorm(residuals(rmSubEenv.lmer))
testZeroInflation(simulateResiduals(rmSubEenv.lmer))

##predictive check
y1<-predict(rmSubEenv.lmer)
predy <- y1 * daw$denscomp.max

dawE <- RmSubE_Env %>% select(year, STATION, denscomp.max, dens.weight.mean, dens.weight.mean.y1, dens.percomp.y1, SAVArea.percomp.y1, SAVArea.y1, ChlA.me, Sal.me, Sal.spme, Sal.y1Dpos, ChlA.spme, TP.y1max, Secc.y1Dme) %>% na.omit()

prd.lm <- lm(dens.weight.mean ~ dens.weight.mean.y1 + predy, data = daw)
summary(prd.lm)

prd.lmer <- lmer(dens.weight.mean ~ dens.weight.mean.y1 + predy + (1|STATION), data = daw)
summary(prd.lmer)
r.squaredGLMM(prd.lmer)

