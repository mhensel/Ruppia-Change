##code to build the datasets needed for Ruppia change analyses, merging the Env data and the Ruppia zone (Station, Subestuary) data over time (formerly Part 1 and 2 of Baywide_RuChange_models.R)

library(tidyverse); library(readxl); library(patchwork);library(beyonce)
library(randomForest); library(leaps); library(lme4)
library(MuMIn); library(DHARMa); library(piecewiseSEM); library(nlme)

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
#denspercomp only... use this for master code file
RmZoneStations_denspercomp <- SAVyear_StationZone %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  #filter(between(Year, 1990, 2019)) %>% use this if you want to create a certain time dataset
  mutate(per.cov = case_when(Density == 1 ~ .05,
                             Density == 2 ~ .25, 
                             Density == 3 ~ .55, 
                             Density == 4 ~ .85)) %>% #convert to density weighted means 
  mutate(dens_cov = SAVAreaHa * per.cov) %>%
  group_by(STATION, Year) %>%
  summarize(dens.weight.mean = sum(dens_cov), SAVArea = sum(SAVAreaHa)) %>%
  mutate(dens.weight.mean.y1 = lag(dens.weight.mean, order_by = Year, k = 1), SAVArea.y1 = lag(SAVArea, order_by = Year, k = 1)) %>%
  full_join(RuppiaStations) %>% group_by(STATION) %>% 
  mutate(dens.percomp = dens.weight.mean/denscomp.max, dens.percomp.y1 = dens.weight.mean.y1/denscomp.max, SAVArea.percomp = SAVArea/RMZoneSAV_HA, SAVArea.percomp.y1 = SAVArea.y1/RMZoneSAV_HA) %>% 
  mutate(dens.percomp.change = (dens.percomp-dens.percomp.y1), SAVArea.percomp.change = (SAVArea.percomp-SAVArea.percomp.y1)) %>% 
  rename("year" = "Year") %>% select(STATION, year, dens.weight.mean, dens.weight.mean.y1, dens.percomp, dens.percomp.y1, dens.percomp.change, SAVArea.percomp.change)

#if you want to rewrite, go for it.
write_csv(RmZoneStations_denspercomp, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZoneStations_denspercomp.csv")
write_csv(RmZoneStations_denspercomp, "./data/RmZoneStations_denspercomp.csv")


#Original RmZoneStations that has a lot more metrics. Dont discard this! 
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
#write_csv(RmZoneStations, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmZoneStations.csv")
#write_csv(RmZoneStations, "../Ruppia/RmZoneStations.csv")


############Part 2: Merging Station Zone with Environmental Data#######
#load in the CBP data, updated to 2019 as of October 2020. this file comes from Ruppia/Data/CBP WQ Station data/tidyCBPWQ data 2019.R

###Merge RuppiaStation coverage data with CBP WQ Station data####
Env_Var_nosumy1sp.CBP_WQ <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_nosumy1sp.CBP_WQ.csv") #more informative and not as unweildy. Has y1 means, summer y1, spring this year. (eliminates y1 spring, summer this year). 

####MASTER ENV AND BAYWIDE DATASET RmZone_Env ####

###Repeat Step 2 with TSSr data####
#here is the TSS only data
Env_Var.tss <-read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var.tss.csv")
#TSSr is RAW. TSS was FINAL_VALUE which is just nothing before 99 so its not in this dataset
#
Env_Var_nosumy1sp.CBP_WQ.tss <- full_join(Env_Var_nosumy1sp.CBP_WQ, Env_Var.tss)

#RmZoneStationZ <-read.csv("./data/RmZoneStations.csv")

####MASTER ENV AND BAYWIDE DATASET with TSS RmZone_Env.tss ####
RmZone_EnvUNFILT <- Env_Var_nosumy1sp.CBP_WQ.tss %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations) #%>%
#filter(denscomp.max > 5) #this gets rid of any zones that are really tiny. about 200 data points less. originally i thought this would clean things up but it actually weakens the SEM? idk about for these graphs

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

write_csv(RmZone_Env.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env.tss.csv")
write_csv(RmZone_Env.tss, "./data/RmZone_Env.tss.csv")

####Our USE dataframes are these babies below. 1/15/21#####
#trimmed the tails, so that stations close to their max or min dont skew the data
RmZone_Env8515.tss <- RmZone_Env.tss %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
write_csv(RmZone_Env8515.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env8515.tss.csv")
write_csv(RmZone_Env8515.tss, "./data/RmZone_Env8515.tss.csv")

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

#Original Workflow: take these babies and go over to Baywide Ruppia RFs.R, and run some random forests. Then come back to the below Part 3






########Spring Means only dataset#######
spme_CBPWQ <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/spme_CBPWQ.csv")
#use the edited RmZoneStations
RmZone_spmeEnvUNFILT <- spme_CBPWQ %>%
  filter(STATION %in% RuppiaStations$STATION) %>%
  full_join(RmZoneStations_denspercomp) #%>%
#filter(denscomp.max > 5) #this gets rid of any zones that are really tiny. about 200 data points less. originally i thought this would clean things up but it actually weakens the SEM? idk about for these graphs

#is.nan.data.frame <- function(x)
#  do.call(cbind, lapply(x, is.nan))
#is.na(RmZone_Env) <- RmZone_Env == "NaN"
RmZone_spmeEnvUNFILT[is.nan(RmZone_spmeEnvUNFILT)] <- 0
is.na(RmZone_spmeEnvUNFILT) <- RmZone_spmeEnvUNFILT == "Inf"
is.na(RmZone_spmeEnvUNFILT) <- RmZone_spmeEnvUNFILT == "-Inf"
RmZone_spmeEnvUNFILT <- as.data.frame(RmZone_spmeEnvUNFILT)

##going to go ahead and do the station filters now: LE3.6, LE3.7,CB5.4W,CB7.1,CB7.1N,CB7.1S,EE3.4, EE3.5 are the stations you cant use any data TSSr before 1999. 
#CB4.1E, CB5.1, CB5.2,CB5.3,LE2.3,EE1.1,EE2.1,EE2.2,EE3.0,EE3.1,EE3.2,EE3.3,ET4.2,ET5.2,ET8.1,ET9.1,LE2.2,RET2.4,LE3.2,LE3.3, LE3.4, CB7.3E are the stations where it is ok to use the TSSr data.
#need to: filter everything before 1999 out on LE3.6, LE3.7,CB5.4W,CB7.1,CB7.1N, CB7.1S,EE3.4, EE3.5

RmZone_spmeEnv1999 <- RmZone_spmeEnvUNFILT %>%
  filter(STATION %in% c("LE3.6", "LE3.7","CB5.4W","CB7.1","CB7.1N", "CB7.1S","EE3.4", "EE3.5")) %>% filter(year > 1999) 

RmZone_spmeEnvFILT <- RmZone_spmeEnvUNFILT %>% #group_by(STATION, year) %>%
  filter(!STATION %in% c("LE3.6", "LE3.7","CB5.4W","CB7.1","CB7.1N", "CB7.1S","EE3.4", "EE3.5"))

RmZone_spmeEnv.tss <- bind_rows(RmZone_spmeEnv1999, RmZone_spmeEnvFILT)

write_csv(RmZone_spmeEnv.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_spmeEnv.tss.csv")
write_csv(RmZone_spmeEnv.tss, "./data/RmZone_spmeEnv.tss.csv")

RmZone_spmeEnv8515.tss<- RmZone_spmeEnv.tss %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
write_csv(RmZone_spmeEnv8515.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmZone_Env8515.tss.csv")
write_csv(RmZone_spmeEnv8515.tss, "./data/RmZone_Env8515.tss.csv")

RmZone_spmeEnv9010.tss <- RmZone_spmeEnv.tss %>%
  filter(!dens.percomp.y1 > 0.90) %>% filter(!dens.percomp.y1 < 0.10)

colSums(is.na(RmZone_spmeEnv.tss))


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

write_csv(RmZoneSubestuaries, "../data/RmZoneSubestuaries.csv")
write_csv(RmZoneSubestuaries, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets//RmZoneSubestuaries.csv")

Subestuary_WQ <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/Subestuary WQ data/Subestuary_WQ.csv") #this comes from tidyCBPWQ data 2019.R
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

write_csv(RmSubE_Env, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE_Env.csv")
write_csv(RmSubE_Env, "./data/RmSubE_Env.csv")

#trimmed the tails, so that Subes close to their max or min dont skew the data
RmSubE_Env8515 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.85) %>% filter(!dens.percomp.y1 < 0.15)
#writing this one bc its what we use for a bunch of the Random Forests below
#write_csv(RmSubE_Env8515, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSubE_Env8515.csv")

RmSubE_Env982 <- RmSubE_Env %>%
  filter(!dens.percomp.y1 > 0.98) %>% filter(!dens.percomp.y1 < 0.02)
write_csv(RmSubE_Env982, "/Volumes/savshare2/Current Projects/Ruppia/Data/RmSEM datasets/RmSubE_Env982.csv")
write_csv(RmSubE_Env982, "./data/RmSubE_Env982.csv")
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

#Original Workflow: take these babies and go over to Baywide Ruppia RFs.R, and run some random forests. Then come back to Part 3 in Baywide RuChange models or go to Baywide Ruppia SEM