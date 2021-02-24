#cleaning and calculating for the Chesapeake Bay Program Water Quality Station data from Rebecca Murphy 2019, in addition the Subestuary dataset is here too
#This file generates the master environmental database for all of the Ruppia analysis currently (12/2020) but also is a master Env database for all of the Chesapeake Bay. 
#Raw data by sampling date can be found CBP_WQ_2019.csv, Env_Var_nosumy1sp.CBP_WQ is what i used for Ruppia analys.
library(tidyverse); library(car)
library(piecewiseSEM); library(lme4)
library(readxl); library(lubridate); library(naniar)
############Station Env data#########
#load data! 
#from Marcs computer
load("~/Documents/R projects/Ruppia/wtemp_sal_sec_chla_to2019.rda")
#from the R drive
load("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/wtemp_sal_sec_chla_to2019.rda")
# data1  contains fairly self-explanatory columns including DATE, STATION, PARAMETER and VALUE. The only processing I did of this set from the data hub was to average duplicates and remove a couple erroneous values

#from Marcs computer
load("~/Documents/R projects/Ruppia/tn_tp_tss_to2019.rda")
#from the R drive
load("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/tn_tp_tss_to2019.rda")
#has a dataframe called data2. This has fewer columns than the other because I pulled this from already processed data sets. Duplicates are averaged also.There are some cases of values reported at <DL or a range due to the components being <DL. Those will be reported here in the middle of the range (or ½ DL). Also, you’ll see the RAW_VALUE and FINAL_VALUE columns. FINAL_VALUE has the adjustments I mentioned in my email below – some data values cut-out, some adjusted due to method changes. All units of these parameters are mg/L. I recommend you use FINAL_VALUE column

LatLongCBPWQ <- read_csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/TidalLongTermStations.csv")

###format the dates so we can get month, day, year seperate (day prob isnt important but whatever)

CBPstationNPSS <- data2 %>% 
  mutate(date = ymd(DATE))
CBPstationNPSS$year <- year(CBPstationNPSS$date)
CBPstationNPSS$month <- month(CBPstationNPSS$date)
CBPstationNPSS$day <- day(CBPstationNPSS$date)

CBPstationWTSSC <- data1 %>% 
  mutate(date = ymd(DATE))
CBPstationWTSSC$year <- year(CBPstationWTSSC$date)
CBPstationWTSSC$month <- month(CBPstationWTSSC$date)
CBPstationWTSSC$day <- day(CBPstationWTSSC$date)

#there might be a better way of doing this but idk how. taking the FINAL_VALUE and converting it to a wide format
CBPtn <- CBPstationNPSS %>%
  unique() %>%
  filter(PARAMETER == "tn") %>%
  mutate(TN= case_when(PARAMETER == "tn" ~ FINAL_VALUE)) %>%
  select(-PARAMETER, -DATE, -LAYER, -FINAL_VALUE, -RAW_VALUE, -LAYER, -DATE) %>%
  filter()

CBPtp <- CBPstationNPSS %>%
  unique() %>%
  filter(PARAMETER == "tp") %>%
  mutate(TP= case_when(PARAMETER == "tp" ~ FINAL_VALUE)) %>%
  select(-PARAMETER, -DATE, -LAYER, -FINAL_VALUE, -RAW_VALUE, -LAYER, -DATE)
CBPtp[!is.na(CBPtp$TP) & CBPtp$TP > 0.6, "TP"] <- 0.6  #deal with outliers there is a couple enormous numbers that arent real

#BAD TSS bc so many NAs from 84-99. RAW_VALUE can be used for some stations: CB4.1E CB5.1 CB5.2CB5.3LE2.3EE1.1EE2.1EE2.2EE3.0EE3.1EE3.2EE3.3ET4.2ET5.2ET8.1ET9.1LE2.2RET2.4LE3.2LE3.3 LE3.4 CB7.3E 
#check this though.
#1/15/21 update: RAW_VALUE and FINAL_VALUE are not different. Final is just an NA before 1999 in all stations. We need that TSS data from the trustworthy stations (will depend on station zone).
CBPtss <- CBPstationNPSS %>%
  unique() %>%
  filter(PARAMETER == "tss") %>%
  mutate(TSSr= case_when(PARAMETER == "tss" ~ RAW_VALUE)) %>%
  mutate(TSS= case_when(PARAMETER == "tss" ~ FINAL_VALUE)) %>%
  select(-PARAMETER, -DATE, -LAYER, -FINAL_VALUE, -RAW_VALUE, -LAYER, -DATE) #%>%
  #filter(TSS > 0, TSSr >0)  #about 20 TSS were negative. use this line and you go to 38593 points (lots of NAs too)

#using the next dataset now. these ones have a DEPTH column that needs to be averaged over
CBPwtemp <- CBPstationWTSSC %>%
  unique() %>%
  filter(PARAMETER == "WTEMP") %>%
  mutate(WTEMP= case_when(PARAMETER == "WTEMP" ~ VALUE)) %>%
  group_by(year, month, day, date, STATION) %>% summarise(WTEMP_me = mean(WTEMP)) %>% ungroup()
#NOTE: 11/17/20- 1984 and 1985 have some stations (eg WT2.1) where only winter months are measured, making any mean calculations really small - MH

CBPsal <- CBPstationWTSSC %>%
  unique() %>%
  filter(PARAMETER == "SALINITY") %>%
  mutate(SALINITY= case_when(PARAMETER == "SALINITY" ~ VALUE)) %>%
  group_by(year, month, day, date, STATION) %>% summarise(SALINITY_me = mean(SALINITY)) %>% ungroup()

CBPsecc <- CBPstationWTSSC %>%
  unique() %>%
  filter(PARAMETER == "SECCHI") %>%
  mutate(SECCHI= case_when(PARAMETER == "SECCHI" ~ VALUE)) %>%
  group_by(year, month, day, date, STATION) %>% summarise(SECCHI_me = mean(SECCHI)) %>% ungroup()

CBPchla <- CBPstationWTSSC %>%
  unique() %>%
  filter(PARAMETER == "CHLA") %>%
  mutate(CHLA= case_when(PARAMETER == "CHLA" ~ VALUE)) %>%
  group_by(year, month, day, date, STATION) %>% summarise(CHLA_me = mean(CHLA)) %>% ungroup() 
CBPchla[CBPchla$CHLA_me > 250, "CHLA_me"] <- 250  #deal with outliers there is a couple enormous numbers that arent real

#join those togehter now. gives complete dataset
#NOTE:11/6/2020 I'm excluding TSS data. it has too many NAs and isnt trustworthy
#NOTE: 1/1/2021 TSS data might be important lol. can filter out the bad ones based on info.
CBP_WQdata1 <- full_join(CBPtn, CBPtp) %>% 
  #full_join(CBPtss) %>%
  full_join(CBPwtemp) %>%
  full_join(CBPsal) %>%
  full_join(CBPsecc) %>%
  full_join(CBPchla)

CBP_WQdata2 <- full_join(CBPtn, CBPtp) %>% 
  full_join(CBPtss) %>%
  full_join(CBPwtemp) %>%
  full_join(CBPsal) %>%
  full_join(CBPsecc) %>%
  full_join(CBPchla)
  
#and then final: TN TP TSS WTEMP SALINITY SECCHI CHLA 
#with Lat/Longs in there too

LatLongCBPWQ <- read_csv("/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/TidalLongTermStations.csv") #also loaded above

#join up and add the Delta Change metric. ChlA.d = change in ChlA between now and last sampling point.

CBP_WQ_2019 <- full_join(CBP_WQdata1, LatLongCBPWQ) %>%
  rename(WTEMP = WTEMP_me, SECCHI = SECCHI_me, CHLA = CHLA_me, SALINITY = SALINITY_me) %>%
  mutate(ChlA.D = CHLA - lag(CHLA), Secc.D = SECCHI - lag(SECCHI), Sal.D = SALINITY - lag(SALINITY), Temp.D = WTEMP - lag(WTEMP), TP.D = TP - lag(TP), TN.D = TN - lag(TN))#, TSS.D = TSS - lag(TSS)) 

CBP_WQ_2019.tss <- full_join(CBP_WQdata2, LatLongCBPWQ) %>%
  rename(WTEMP = WTEMP_me, SECCHI = SECCHI_me, CHLA = CHLA_me, SALINITY = SALINITY_me) %>%
  mutate(ChlA.D = CHLA - lag(CHLA), Secc.D = SECCHI - lag(SECCHI), Sal.D = SALINITY - lag(SALINITY), Temp.D = WTEMP - lag(WTEMP), TP.D = TP - lag(TP), TN.D = TN - lag(TN), TSS.D = TSS - lag(TSS), TSSr.D = TSSr - lag(TSSr)) 

#%>% mutate(ChlA.y1 = lag(CHLA), Secc.y1 = lag(SECCHI), Sal.y1 = lag(SALINITY), Temp.y1 = lag(WTEMP), TP.y1 = lag(TP), TN.y1 = lag(TN), TSS.y1 = lag(TSS))

#examine the stations: are these right? some look weird
unique(CBP_WQ_2019$STATION) #Rebecca said 145 STATIONS so this fits

CBP_WQ_2019 %>% summarise_all(funs(sum(is.na(.))))
sum(is.na(CBP_WQ_2019$TSS))

#write it back into the R drive. this is RAW, w no filtering
write_csv(CBP_WQ_2019, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/CBP_WQ_2019.csv")
write_csv(CBP_WQ_2019.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/CBP_WQ_2019.tss.csv")

#write to Marcs computer (dont tell dave ;) )
getwd()
write_csv(CBP_WQ_2019, "../Ruppia/CBP_WQ_2019.csv")


ggplot(data = CBP_WQ_2019 %>% filter(STATION == "CB5.1")) + 
  stat_summary(aes(x = date, y = WTEMP, color = WTEMP), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .4) +
  stat_summary(aes(x = date, y = WTEMP), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = .6, color = "black") +
  scale_color_gradient(low="blue", high="red") +
  #ylim(10,25) +
  ylab("Temperature (daily C)") + 
  theme_bw(base_size=14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "")


####Calculating Env Variables of interest MH######
#What do we want to know: BIG swings in Env variables and on Ru coverage
#Monthly max and min for everything.
#any time that it deviated bigly from the monthly mean
#NOTE ON TSS: the code to create the TSS dataset is in there but hashtagged out. that file is in the folders too. 

#calculate yearly max and mins and means and range. also calculate D change max min mean
Env_Var_year.CBP_WQ <- CBP_WQ_2019 %>%
  # na.omit() %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(ChlA.max = max(CHLA, na.rm = T), ChlA.min = min(CHLA, na.rm = T), ChlA.me = mean(CHLA, na.rm = T), ChlA.ran = ChlA.max - ChlA.min, 
            ChlA.Dpos = max(ChlA.D, na.rm = T), ChlA.Dneg = min(ChlA.D, na.rm = T), ChlA.Dme = mean(ChlA.D, na.rm = T), 
            Secc.max = max(SECCHI, na.rm = T), Secc.min = min(SECCHI, na.rm = T), Secc.me = mean(SECCHI, na.rm = T), Secc.ran = Secc.max - Secc.min, Secc.Dpos = max(Secc.D, na.rm = T), Secc.Dneg = min(Secc.D, na.rm = T), Secc.Dme = mean(Secc.D, na.rm = T),
            Sal.max = max(SALINITY, na.rm = T), Sal.min = min(SALINITY, na.rm = T), Sal.me = mean(SALINITY, na.rm = T), Sal.ran = Sal.max - Sal.min, Sal.Dpos = max(Sal.D, na.rm = T), Sal.Dneg = min(Sal.D, na.rm = T), Sal.Dme = mean(Sal.D, na.rm = T),
            Temp.max = max(WTEMP, na.rm = T), Temp.min = min(WTEMP, na.rm = T), Temp.me = mean(WTEMP, na.rm = T), Temp.ran = Temp.max - Temp.min, Temp.Dpos = max(Temp.D, na.rm = T), Temp.Dneg = min(Temp.D, na.rm = T), Temp.Dme = mean(Temp.D, na.rm = T),
            TP.max = max(TP, na.rm = T), TP.min = min(TP, na.rm = T), TP.me = mean(TP, na.rm = T), TP.ran = TP.max - TP.min, 
            TP.Dpos = max(TP.D, na.rm = T), TP.Dneg = min(TP.D, na.rm = T), TP.Dme = mean(TP.D, na.rm = T),
            TN.max = max(TN, na.rm = T), TN.min = min(TN, na.rm = T), TN.me = mean(TN, na.rm = T), TN.ran = TN.max - TN.min,
            TN.Dpos = max(TN.D, na.rm = T), TN.Dneg = min(TN.D, na.rm = T), TN.Dme = mean(TN.D, na.rm = T))
            #,TSS.max = max(TSS, na.rm = T), TSS.min = min(TSS, na.rm = T), TSS.me = mean(TSS, na.rm = T), TSS.ran = TSS.max - TSS.min, 
            #TSS.Dpos = max(TSS.D, na.rm = T), TSS.Dneg = min(TSS.D, na.rm = T), TSS.Dme = mean(TSS.D, na.rm = T))
         

write_csv(Env_Var_year.CBP_WQ, "../Ruppia/Env_Var_year.CBP_WQ.csv") 

#calculate y-1 yearly
Env_Var_yearly <- Env_Var_year.CBP_WQ %>% group_by(STATION) %>%
  mutate(ChlA.y1max = lag(ChlA.max), ChlA.y1min = lag(ChlA.min), ChlA.y1me = lag(ChlA.me), ChlA.y1ran = lag(ChlA.ran), ChlA.y1Dpos = lag(ChlA.Dpos), ChlA.y1Dneg = lag(ChlA.Dneg), ChlA.y1Dme = lag(ChlA.Dme),
         Secc.y1max = lag(Secc.max), Secc.y1min = lag(Secc.min), Secc.y1me = lag(Secc.me), Secc.y1ran = lag(Secc.ran), Secc.y1Dpos = lag(Secc.Dpos), Secc.y1Dneg = lag(Secc.Dneg), Secc.y1Dme = lag(Secc.Dme),
         Sal.y1max = lag(Sal.max), Sal.y1min = lag(Sal.min), Sal.y1me = lag(Sal.me), Sal.y1ran = lag(Sal.ran), Sal.y1Dpos = lag(Sal.Dpos), Sal.y1Dneg = lag(Sal.Dneg), Sal.y1Dme = lag(Sal.Dme), 
         Temp.y1max = lag(Temp.max), Temp.y1min = lag(Temp.min), Temp.y1me = lag(Temp.me), Temp.y1ran = lag(Temp.ran), Temp.y1Dpos = lag(Temp.Dpos), Temp.y1Dneg = lag(Temp.Dneg), Temp.y1Dme = lag(Temp.Dme),
         TN.y1max = lag(TN.max), TN.y1min = lag(TN.min), TN.y1me = lag(TN.me), TN.y1ran = lag(TN.ran), TN.y1Dpos = lag(TN.Dpos), TN.y1Dneg = lag(TN.Dneg), TN.y1Dme = lag(TN.Dme), 
         TP.y1max = lag(TP.max), TP.y1min = lag(TP.min), TP.y1me = lag(TP.me), TP.y1ran = lag(TP.ran), TP.y1Dpos = lag(TP.Dpos), TP.y1Dneg = lag(TP.Dneg), TP.y1Dme = lag(TP.Dme))
         #,
        # TSS.y1max = lag(TSS.max), TSS.y1min = lag(TSS.min), TSS.y1me = lag(TSS.me), TSS.y1ran = lag(TSS.ran), TSS.y1Dpos = lag(TSS.Dpos), TSS.y1Dneg = lag(TSS.Dneg), TSS.y1Dme = lag(TSS.Dme))

#checkyaself, just chekcing that this shit works
#LE2.2day <- Env_Var_year.CBP_WQ %>% filter(STATION == "LE2.2")
#LE2.2 <- CBP_WQ_2019 %>% filter(STATION == "LE2.2")
#CBP_WQ_2019 %>% filter(STATION == "LE2.2") %>% filter(year == "2000") %>% summarize(mee = mean(ChlA.D))

colSums(is.na(Env_Var_AprAug.CBP_WQ))
#calc growing season means (called "sum" as in summer but you get the point)
Env_Var_AprAug.CBP_WQ <- CBP_WQ_2019 %>%
  filter(between(month, 4, 8)) %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(ChlA.summax = max(CHLA, na.rm = T), ChlA.summin = min(CHLA, na.rm = T), ChlA.summe = mean(CHLA, na.rm = T), ChlA.sumran = ChlA.summax - ChlA.summin, 
            Secc.summax = max(SECCHI, na.rm = T), Secc.summin = min(SECCHI, na.rm = T), Secc.summe = mean(SECCHI, na.rm = T), Secc.sumran = Secc.summax - Secc.summin, 
            Sal.summax = max(SALINITY, na.rm = T), Sal.summin = min(SALINITY, na.rm = T), Sal.summe = mean(SALINITY, na.rm = T), Sal.sumran = Sal.summax - Sal.summin,
            Temp.summax = max(WTEMP, na.rm = T), Temp.summin = min(WTEMP, na.rm = T), Temp.summe = mean(WTEMP, na.rm = T), Temp.sumran = Temp.summax - Temp.summin,
            TP.summax = max(TP, na.rm = T), TP.summin = min(TP, na.rm = T), TP.summe = mean(TP, na.rm = T), TP.sumran = TP.summax - TP.summin,
            TN.summax = max(TN, na.rm = T), TN.summin = min(TN, na.rm = T), TN.summe = mean(TN, na.rm = T), TN.sumran = TN.summax - TN.summin,
            #TSS.summax = max(TSS, na.rm = T), TSS.summin = min(TSS, na.rm = T), TSS.summe = mean(TSS, na.rm = T), TSS.sumran = TSS.summax - TSS.summin, 
            ChlA.sumDpos = max(ChlA.D, na.rm = T), ChlA.sumDneg = min(ChlA.D, na.rm = T), ChlA.sumDme = mean(ChlA.D, na.rm = T), 
            Secc.sumDpos = max(Secc.D, na.rm = T), Secc.sumDneg = min(Secc.D, na.rm = T), Secc.sumDme = mean(Secc.D, na.rm = T),
            Sal.sumDpos = max(Sal.D, na.rm = T), Sal.sumDneg = min(Sal.D, na.rm = T), Sal.sumDme = mean(Sal.D, na.rm = T),
            Temp.sumDpos = max(Temp.D, na.rm = T), Temp.sumDneg = min(Temp.D, na.rm = T), Temp.sumDme = mean(Temp.D, na.rm = T),
            TP.sumDpos = max(TP.D, na.rm = T), TP.sumDneg = min(TP.D, na.rm = T), TP.sumDme = mean(TP.D, na.rm = T),
            TN.sumDpos = max(TN.D, na.rm = T), TN.sumDneg = min(TN.D, na.rm = T), TN.sumDme = mean(TN.D, na.rm = T))
            #,
            #TSS.sumDpos = max(TSS.D, na.rm = T), TSS.sumDneg = min(TSS.D, na.rm = T), TSS.sumDme = mean(TSS.D, na.rm = T))
#growing season y-1
Env_Var_grow <- Env_Var_AprAug.CBP_WQ %>% group_by(STATION) %>%
  mutate(ChlA.sumy1max = lag(ChlA.summax), ChlA.sumy1min = lag(ChlA.summin), ChlA.sumy1me = lag(ChlA.summe), ChlA.sumy1ran = lag(ChlA.sumran), ChlA.sumy1Dpos = lag(ChlA.sumDpos), ChlA.sumy1Dneg = lag(ChlA.sumDneg), ChlA.sumy1Dme = lag(ChlA.sumDme),
         Secc.sumy1max = lag(Secc.summax), Secc.sumy1min = lag(Secc.summin), Secc.sumy1me = lag(Secc.summe), Secc.sumy1ran = lag(Secc.sumran), Secc.sumy1Dpos = lag(Secc.sumDpos), Secc.sumy1Dneg = lag(Secc.sumDneg), Secc.sumy1Dme = lag(Secc.sumDme),
         Sal.sumy1max = lag(Sal.summax), Sal.sumy1min = lag(Sal.summin), Sal.sumy1me = lag(Sal.summe), Sal.sumy1ran = lag(Sal.sumran), Sal.sumy1Dpos = lag(Sal.sumDpos), Sal.sumy1Dneg = lag(Sal.sumDneg), Sal.sumy1Dme = lag(Sal.sumDme), 
         Temp.sumy1max = lag(Temp.summax), Temp.sumy1min = lag(Temp.summin), Temp.sumy1me = lag(Temp.summe), Temp.sumy1ran = lag(Temp.sumran), Temp.sumy1Dpos = lag(Temp.sumDpos), Temp.sumy1Dneg = lag(Temp.sumDneg), Temp.sumy1Dme = lag(Temp.sumDme),
         TN.sumy1max = lag(TN.summax), TN.sumy1min = lag(TN.summin), TN.sumy1me = lag(TN.summe), TN.sumy1ran = lag(TN.sumran), TN.sumy1Dpos = lag(TN.sumDpos), TN.sumy1Dneg = lag(TN.sumDneg), TN.sumy1Dme = lag(TN.sumDme), 
         TP.sumy1max = lag(TP.summax), TP.sumy1min = lag(TP.summin), TP.sumy1me = lag(TP.summe), TP.sumy1ran = lag(TP.sumran), TP.sumy1Dpos = lag(TP.sumDpos), TP.sumy1Dneg = lag(TP.sumDneg), TP.sumy1Dme = lag(TP.sumDme)) #,
         #TSS.sumy1max = lag(TSS.summax), TSS.sumy1min = lag(TSS.summin), TSS.sumy1me = lag(TSS.summe), TSS.sumy1ran = lag(TSS.sumran), TSS.sumy1Dpos = lag(TSS.sumDpos), TSS.sumy1Dneg = lag(TSS.sumDneg), TSS.sumy1Dme = lag(TSS.sumDme))

#need a just grow season y1 because grow season of Y might be uninformative, esp compared to spring
Env_Var_sumy1 <- Env_Var_grow %>% select(-ChlA.summax, -ChlA.summin , -ChlA.summe, -ChlA.sumran, -Secc.summax , -Secc.summin , -Secc.summe, -Secc.sumran, -Sal.summax, -Sal.summin, -Sal.summe, -Sal.sumran, -Temp.summax, -Temp.summin, -Temp.summe, -Temp.sumran, -TP.summax, -TP.summin, -TP.summe, -TP.sumran, -TN.summax, -TN.summin, -TN.summe, -TN.sumran,  -ChlA.sumDpos, -ChlA.sumDneg , -ChlA.sumDme, -Secc.sumDpos, -Secc.sumDneg, -Secc.sumDme, -Sal.sumDpos, -Sal.sumDneg, -Sal.sumDme, -Temp.sumDpos, -Temp.sumDneg, -Temp.sumDme, -TP.sumDpos, -TP.sumDneg, -TP.sumDme, -TN.sumDpos, -TN.sumDneg, -TN.sumDme) #, -TSS.sumDpos, -TSS.sumDneg, -TSS.sumDme, -TSS.summax, -TSS.summin, -TSS.summe, -TSS.sumran)

##calc spring means
Env_Var_MarchMay.CBP_WQ <- CBP_WQ_2019 %>% 
  filter(between(month, 3, 5)) %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(ChlA.spmax = max(CHLA, na.rm = T), ChlA.spmin = min(CHLA, na.rm = T), ChlA.spme = mean(CHLA, na.rm = T), ChlA.spran = ChlA.spmax - ChlA.spmin, 
            Secc.spmax = max(SECCHI, na.rm = T), Secc.spmin = min(SECCHI, na.rm = T), Secc.spme = mean(SECCHI, na.rm = T), Secc.spran = Secc.spmax - Secc.spmin, 
            Sal.spmax = max(SALINITY, na.rm = T), Sal.spmin = min(SALINITY, na.rm = T), Sal.spme = mean(SALINITY, na.rm = T), Sal.spran = Sal.spmax - Sal.spmin,
            Temp.spmax = max(WTEMP, na.rm = T), Temp.spmin = min(WTEMP, na.rm = T), Temp.spme = mean(WTEMP, na.rm = T), Temp.spran = Temp.spmax - Temp.spmin,
            TP.spmax = max(TP, na.rm = T), TP.spmin = min(TP, na.rm = T), TP.spme = mean(TP, na.rm = T), TP.spran = TP.spmax - TP.spmin,
            TN.spmax = max(TN, na.rm = T), TN.spmin = min(TN, na.rm = T), TN.spme = mean(TN, na.rm = T), TN.spran = TN.spmax - TN.spmin,
            #TSS.spmax = max(TSS, na.rm = T), TSS.spmin = min(TSS, na.rm = T), TSS.spme = mean(TSS, na.rm = T), TSS.spran = TSS.spmax - TSS.spmin, 
            ChlA.spDpos = max(ChlA.D, na.rm = T), ChlA.spDneg = min(ChlA.D, na.rm = T), ChlA.spDme = mean(ChlA.D, na.rm = T), 
            Secc.spDpos = max(Secc.D, na.rm = T), Secc.spDneg = min(Secc.D, na.rm = T), Secc.spDme = mean(Secc.D, na.rm = T),
            Sal.spDpos = max(Sal.D, na.rm = T), Sal.spDneg = min(Sal.D, na.rm = T), Sal.spDme = mean(Sal.D, na.rm = T),
            Temp.spDpos = max(Temp.D, na.rm = T), Temp.spDneg = min(Temp.D, na.rm = T), Temp.spDme = mean(Temp.D, na.rm = T),
            TP.spDpos = max(TP.D, na.rm = T), TP.spDneg = min(TP.D, na.rm = T), TP.spDme = mean(TP.D, na.rm = T),
            TN.spDpos = max(TN.D, na.rm = T), TN.spDneg = min(TN.D, na.rm = T), TN.spDme = mean(TN.D, na.rm = T))
            #,
           # TSS.spDpos = max(TSS.D, na.rm = T), TSS.spDneg = min(TSS.D, na.rm = T), TSS.spDme = mean(TSS.D, na.rm = T))

##Spring y-1
Env_Var_spring <- Env_Var_MarchMay.CBP_WQ %>% group_by(STATION) %>%
  mutate(ChlA.spy1max = lag(ChlA.spmax), ChlA.spy1min = lag(ChlA.spmin), ChlA.spy1me = lag(ChlA.spme), ChlA.spy1ran = lag(ChlA.spran), ChlA.spy1Dpos = lag(ChlA.spDpos), ChlA.spy1Dneg = lag(ChlA.spDneg), ChlA.spy1Dme = lag(ChlA.spDme),
         Secc.spy1max = lag(Secc.spmax), Secc.spy1min = lag(Secc.spmin), Secc.spy1me = lag(Secc.spme), Secc.spy1ran = lag(Secc.spran), Secc.spy1Dpos = lag(Secc.spDpos), Secc.spy1Dneg = lag(Secc.spDneg), Secc.spy1Dme = lag(Secc.spDme),
         Sal.spy1max = lag(Sal.spmax), Sal.spy1min = lag(Sal.spmin), Sal.spy1me = lag(Sal.spme), Sal.spy1ran = lag(Sal.spran), Sal.spy1Dpos = lag(Sal.spDpos), Sal.spy1Dneg = lag(Sal.spDneg), Sal.spy1Dme = lag(Sal.spDme), 
         Temp.spy1max = lag(Temp.spmax), Temp.spy1min = lag(Temp.spmin), Temp.spy1me = lag(Temp.spme), Temp.spy1ran = lag(Temp.spran), Temp.spy1Dpos = lag(Temp.spDpos), Temp.spy1Dneg = lag(Temp.spDneg), Temp.spy1Dme = lag(Temp.spDme),
         TN.spy1max = lag(TN.spmax), TN.spy1min = lag(TN.spmin), TN.spy1me = lag(TN.spme), TN.spy1ran = lag(TN.spran), TN.spy1Dpos = lag(TN.spDpos), TN.spy1Dneg = lag(TN.spDneg), TN.spy1Dme = lag(TN.spDme), 
         TP.spy1max = lag(TP.spmax), TP.spy1min = lag(TP.spmin), TP.spy1me = lag(TP.spme), TP.spy1ran = lag(TP.spran), TP.spy1Dpos = lag(TP.spDpos), TP.spy1Dneg = lag(TP.spDneg), TP.spy1Dme = lag(TP.spDme))
       #  ,TSS.spy1max = lag(TSS.spmax), TSS.spy1min = lag(TSS.spmin), TSS.spy1me = lag(TSS.spme), TSS.spy1ran = lag(TSS.spran), TSS.spy1Dpos = lag(TSS.spDpos), TSS.spy1Dneg = lag(TSS.spDneg), TSS.spy1Dme = lag(TSS.spDme))

#merge all the CBP data together. this should have 7 (Vars) * 7 (max mins delt change) * 6 (time frame) variables. Merging the Yearly, growing season, and spring 

Env_Var_ALL.CBP_WQ <- full_join(Env_Var_yearly, Env_Var_grow) %>% full_join(Env_Var_spring) %>%
  mutate(STATION = replace(STATION, STATION == "LE5.5", "LE5.5-W"))

is.na(Env_Var_ALL.CBP_WQ) <- Env_Var_ALL.CBP_WQ == "NaN"
is.na(Env_Var_ALL.CBP_WQ) <- Env_Var_ALL.CBP_WQ == "Inf"
is.na(Env_Var_ALL.CBP_WQ) <- Env_Var_ALL.CBP_WQ == "-Inf"


#here is the "more informative" one with grow season and spring y1 out of it
Env_Var_nosumy1sp.CBP_WQ <- full_join(Env_Var_yearly, Env_Var_sumy1) %>% 
  full_join(Env_Var_MarchMay.CBP_WQ) %>%
  mutate(STATION = replace(STATION, STATION == "LE5.5", "LE5.5-W"))

is.na(Env_Var_nosumy1sp.CBP_WQ) <- Env_Var_nosumy1sp.CBP_WQ == "NaN"
is.na(Env_Var_nosumy1sp.CBP_WQ) <- Env_Var_nosumy1sp.CBP_WQ == "Inf"
is.na(Env_Var_nosumy1sp.CBP_WQ) <- Env_Var_nosumy1sp.CBP_WQ == "-Inf"

colSums(is.na(Env_Var_ALL.CBP_WQ))
colSums(is.na(Env_Var_nosumy1sp.CBP_WQ))

write_csv(Env_Var_ALL.CBP_WQ, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_ALL.CBP_WQ.csv")
write_csv(Env_Var_ALL.CBP_WQ, "../Ruppia/Env_Var_ALL.CBP_WQ.csv")

write_csv(Env_Var_nosumy1sp.CBP_WQ, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_nosumy1sp.CBP_WQ.csv")
write_csv(Env_Var_nosumy1sp.CBP_WQ, "../Ruppia/Env_Var.TSS_nosumy1sp.CBP_WQ.csv")

Env_Var_nosumy1sp.CBP_WQ

#Incorporating TSS into the Env_Var#
#NOTES: to use this .tss data (all data will have .tss at the end), you should filter OUT the TSS numbers because those are the FINAL_VALUE column that is just an NA for pre 1999. Not all stations are right for the TSS.r data, but enough are that throwing out all pre 99 TSS data isnt necessary
#DECISION TIME 1/15/21- im hashtagging all the TSS data and bringing TSS.r data only over to the analyses. If you want the FINAL_VALUE data, just filter(TSS.r > 1999)

####Station Zone with TSS data####
#is there a better way to do this than to make the whole dataset over again? im not sure.

Env_Var_year.CBP_WQ.tss <- CBP_WQ_2019.tss %>%
  # na.omit() %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(#TSS.max = max(TSS, na.rm = T), TSS.min = min(TSS, na.rm = T), TSS.me = mean(TSS, na.rm = T), TSS.ran = TSS.max - TSS.min, 
            #TSS.Dpos = max(TSS.D, na.rm = T), TSS.Dneg = min(TSS.D, na.rm = T), TSS.Dme = mean(TSS.D, na.rm = T),
            TSSr.max = max(TSSr, na.rm = T), TSSr.min = min(TSSr, na.rm = T), TSSr.me = mean(TSSr, na.rm = T), TSSr.ran = TSSr.max - TSSr.min, 
            TSSr.Dpos = max(TSSr.D, na.rm = T), TSSr.Dneg = min(TSSr.D, na.rm = T), TSSr.Dme = mean(TSSr.D, na.rm = T))

#calculate y-1 yearly
Env_Var_yearly.tss <- Env_Var_year.CBP_WQ.tss %>% group_by(STATION) %>%
  mutate(#TSS.y1max = lag(TSS.max), TSS.y1min = lag(TSS.min), TSS.y1me = lag(TSS.me), TSS.y1ran = lag(TSS.ran), TSS.y1Dpos = lag(TSS.Dpos), TSS.y1Dneg = lag(TSS.Dneg), TSS.y1Dme = lag(TSS.Dme),
         TSSr.y1max = lag(TSSr.max), TSSr.y1min = lag(TSSr.min), TSSr.y1me = lag(TSSr.me), TSSr.y1ran = lag(TSSr.ran), TSSr.y1Dpos = lag(TSSr.Dpos), TSSr.y1Dneg = lag(TSSr.Dneg), TSSr.y1Dme = lag(TSSr.Dme))


#calc growing season means (called "sum" as in summer but you get the point)
Env_Var_AprAug.CBP_WQ.tss <- CBP_WQ_2019.tss %>%
  filter(between(month, 4, 8)) %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(#TSS.summax = max(TSS, na.rm = T), TSS.summin = min(TSS, na.rm = T), TSS.summe = mean(TSS, na.rm = T), TSS.sumran = TSS.summax - TSS.summin, 
            TSSr.summax = max(TSSr, na.rm = T), TSSr.summin = min(TSSr, na.rm = T), TSSr.summe = mean(TSSr, na.rm = T), TSSr.sumran = TSSr.summax - TSSr.summin,
            #TSS.sumDpos = max(TSS.D, na.rm = T), TSS.sumDneg = min(TSS.D, na.rm = T), TSS.sumDme = mean(TSS.D, na.rm = T), 
            TSSr.sumDpos = max(TSSr.D, na.rm = T), TSSr.sumDneg = min(TSSr.D, na.rm = T), TSSr.sumDme = mean(TSSr.D, na.rm = T))

#growing season y-1
Env_Var_grow.tss <- Env_Var_AprAug.CBP_WQ.tss %>% group_by(STATION) %>%
  mutate(#TSS.sumy1max = lag(TSS.summax), TSS.sumy1min = lag(TSS.summin), TSS.sumy1me = lag(TSS.summe), TSS.sumy1ran = lag(TSS.sumran), TSS.sumy1Dpos = lag(TSS.sumDpos), TSS.sumy1Dneg = lag(TSS.sumDneg), TSS.sumy1Dme = lag(TSS.sumDme), 
         TSSr.sumy1max = lag(TSSr.summax), TSSr.sumy1min = lag(TSSr.summin), TSSr.sumy1me = lag(TSSr.summe), TSSr.sumy1ran = lag(TSSr.sumran), TSSr.sumy1Dpos = lag(TSSr.sumDpos), TSSr.sumy1Dneg = lag(TSSr.sumDneg), TSSr.sumy1Dme = lag(TSSr.sumDme))

#need a just grow season y1 because grow season of Y might be uninformative, esp compared to spring
Env_Var_sumy1.tss <- Env_Var_grow.tss %>% select(#-TSS.sumDpos, -TSS.sumDneg, -TSS.sumDme, -TSS.summax, -TSS.summin, -TSS.summe, -TSS.sumran, 
  -TSSr.sumDpos, -TSSr.sumDneg, -TSSr.sumDme, -TSSr.summax, -TSSr.summin, -TSSr.summe, -TSSr.sumran)

##calc spring means
Env_Var_MarchMay.CBP_WQ.tss <- CBP_WQ_2019.tss %>% 
  filter(between(month, 3, 5)) %>%
  group_by(year, STATION, LATITUDE, LONGITUDE, CBSEG_2003) %>% 
  summarise(#TSS.spmax = max(TSS, na.rm = T), TSS.spmin = min(TSS, na.rm = T), TSS.spme = mean(TSS, na.rm = T), TSS.spran = TSS.spmax - TSS.spmin, 
            TSSr.spmax = max(TSSr, na.rm = T), TSSr.spmin = min(TSSr, na.rm = T), TSSr.spme = mean(TSSr, na.rm = T), TSSr.spran = TSSr.spmax - TSSr.spmin, 
            #TSS.spDpos = max(TSS.D, na.rm = T), TSS.spDneg = min(TSS.D, na.rm = T), TSS.spDme = mean(TSS.D, na.rm = T), 
            TSSr.spDpos = max(TSSr.D, na.rm = T), TSSr.spDneg = min(TSSr.D, na.rm = T), TSSr.spDme = mean(TSSr.D, na.rm = T))

##Spring y-1
#note: skipped 1/14/21
#Env_Var_spring.tss <- Env_Var_MarchMay.CBP_WQ.tss %>% group_by(STATION) %>%
#  mutate(TSS.spy1max = lag(TSS.spmax), TSS.spy1min = lag(TSS.spmin), TSS.spy1me = lag(TSS.spme), TSS.spy1ran = lag(TSS.spran), TSS.spy1Dpos = lag(TSS.spDpos), TSS.spy1Dneg = lag(TSS.spDneg), TSS.spy1Dme = lag(TSS.spDme))

#master TSS dataset. TSSr is RAW_VALUE Env_Var.tss####
#tss use this data
#here is the "more informative" one with grow season and spring y1 out of it
Env_Var.tss <- full_join(Env_Var_yearly.tss, Env_Var_sumy1.tss) %>% 
  full_join(Env_Var_MarchMay.CBP_WQ.tss) %>%
  mutate(STATION = replace(STATION, STATION == "LE5.5", "LE5.5-W"))

is.na(Env_Var.tss) <- Env_Var.tss == "NaN"
is.na(Env_Var.tss) <- Env_Var.tss == "Inf"
is.na(Env_Var.tss) <- Env_Var.tss == "-Inf"

write_csv(Env_Var.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var.tss.csv")

#merge into the other dataset below, but prob want to save this code for later
Env_Var_nosumy1sp.CBP_WQ.tss <- full_join(Env_Var_nosumy1sp.CBP_WQ, Env_Var.tss)

write_csv(Env_Var_nosumy1sp.CBP_WQ.tss, "/Volumes/savshare2/Current Projects/Ruppia/Data/CBP WQ Station data/Env_Var_nosumy1sp.CBP_WQ.tss.csv")
write_csv(Env_Var_nosumy1sp.CBP_WQ.tss, "../Ruppia/Env_Var_nosumy1sp.CBP_WQ.tss.tss.csv")

write_csv(Env_Var_nosumy1sp.CBP_WQ.tss, "./_data/Baywide SAV/WaterQual_CBP2019.csv")

#
##
###
####
############Subesturay Env data#########
####
###
##
#

SubEs_Env <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/subestEnv15.csv") #this has only flow and N/P up until 2015, no 13/14 for manure and others

restime.SubE <- read.csv("/Volumes/savshare2/Current Projects/Ruppia/Data/residence_time.subest.2csv.csv")

#clean this up and merge with the restime document 
SubEs_flowmerge <- full_join(SubEs_Env, restime.SubE, by = c("SUBEST_ID")) %>%
  group_by(Year, SUBEST_ID) %>% select(-X, -Area_HA, -Area_HA1, -Area_HA3, -Area_HA5, -SalinityZone, -ew, -Richness, -Hectares1MPlusSAVComp, -log10.Area_HA, -log10.Area_HA1, -log10.Area_HA3, -log10.Area_HA5, -log10.Richness, -log10.Hectares1MPlusSAVComp) %>% ungroup() #%>%na.omit() #minus 1300 w na.omit



#create function to calculate adjusted loads, based on flow out of the estuaries.
#rate is per day
#this function can be used generally
  adjust_load_gen <- function(load, vol.m3, Flush1, Flush2) {
    adj_load <- ((load - (365 - mean(c(Flush1, Flush2))*load) / vol.m3))
    return(adj_load)
  }


  #group_by(SUBEST_ID, Year) %>%
  #select(flow, nptn, nptp, ptn, tssx, manuren_kg, manurep_kg, fertilizern_kg, fertilizerp_kg, log10.flow, log10.tssx, log10.nptp, log10.nptp, log10.ptn, log10.ptp, log10.manureN, log10.manureP, log10.fertilizerN, log10.fertilizerP) %>% na.omit() %>%
  
#adjusted load formula applied here. the y1 variables are made from adjusted load, i didnt do the .adj tho. kinda annoying i suppose 
Subestuary_WQ <- SubEs_flowmerge %>% 
  select(-log10.flow, -log10.fertilizerN, -log10.tssx, -log10.fertilizerP, -log10.nptp, -log10.nptn, -log10.WatershedHa, -log10.ptn, -log10.ptp, -log10.manureN, -log10.manureP) %>% 
  group_by(SUBEST_ID, Year) %>%
  mutate(Vol.km3 = VOLUMm3/1000000000, #km data rounds dwn to 0
 #        Flush.me = mean(TidalPFlushTime, TidalPFlushTime.RT0.1), #tidal flush time is in days
         P = (TidalPrism/12.42)*24) %>% #P = tidal exchange rate (m3/day) %>% 
  #mutate(RT = (365-Flush.me)/365) %>% #RT
  mutate(VoFlom3 = VOLUMm3 + flow, FP = flow + P) %>% 
  select(-VOLUMkm3, -TidalPFlushTime, -TidalPFlushTime.RT0.1) %>% 
  mutate(nptn.adj = (flow * nptn)/FP, nptp.adj = (flow * nptp)/FP, ptn.adj = (flow * ptn)/FP, ptp.adj = (flow * ptp)/FP, tssx.adj = (flow * tssx)/FP, manuren_kg.adj = (flow * manuren_kg)/FP, manurep_kg.adj = (flow * manurep_kg)/FP, fertilizern_kg.adj = (flow * fertilizern_kg)/FP, fertilizerp_kg.adj = (flow * fertilizerp_kg)/FP) %>% ungroup()

#new don formula: conc = (flow rate * material concentration (nptp/vol))/flow + tidal exchange rate (P)
  
#%>% 
#  mutate(nptn.adj = (nptn-(RT*nptn))/Vol.km3, nptp.adj = (nptp-(RT*nptp))/Vol.km3, ptn.adj = (ptn-(RT*ptn))/Vol.km3, ptp.adj = (ptp-(RT*ptp))/Vol.km3, tssx.adj = (tssx-(RT*tssx))/Vol.km3, manuren_kg.adj = (manuren_kg-(RT*manuren_kg))/Vol.km3, manurep_kg.adj = (manurep_kg-(RT*manurep_kg))/Vol.km3, fertilizern_kg.adj = (fertilizern_kg-(RT*fertilizern_kg))/Vol.km3, fertilizerp_kg.adj = (fertilizerp_kg-(RT*fertilizerp_kg))/Vol.km3) %>% ungroup() 

Subestuary_WQ[Subestuary_WQ$tssx.adj > 5.e+11, "tssx.adj"] <- 5.e+05

Subestuary_WQ<- Subestuary_WQ %>%
  mutate(nptn.y1 = lag(nptn.adj), nptp.y1 = lag(nptp.adj), ptn.y1 = lag(ptn.adj), ptp.y1 = lag(ptp.adj), tssx.y1 = lag(tssx.adj), manuren_kg.y1 = lag(manuren_kg.adj), manurep_kg.y1 = lag(manurep_kg.adj), fertilizern_kg.y1 = lag(fertilizern_kg.adj), fertilizerp_kg.y1 = lag(fertilizerp_kg.adj), Developed.y1 = lag(Developed), Agro.y1 = lag(Agro), flow.y1 = lag(flow))

##there are a couple tss.x outliers that need fixing

Subestuary_WQ[Subestuary_WQ$tssx.y1 > 5.e+11, "tssx.y1"] <- 5.e+05

write_csv(Subestuary_WQ, "/Volumes/savshare2/Current Projects/Ruppia/Data/Subestuary_WQ.csv")
write_csv(Subestuary_WQ, "../Ruppia/data/Subestuary_WQ.csv")


########Connowingo Dam data########
connodisc <- read.csv("../Ruppia/data/conno_discharge.csv")

conno_monthly <- 
  ggplot(data = connodisc) + 
  stat_summary(aes(x = year_nu, y = discharge_cufs, group = month_nu, color = month_nu), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9) +
  stat_summary(aes(x = year_nu, y = discharge_cufs), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = 1.5, color = "black") +
  #geom_smooth(aes(x = year, y = dens.percomp.change, group = STATION, color = STATION), method = "lm", alpha = 0.5, size = .5) +
  theme_bw(base_size=20) + 
  ylab("Mean discharge") + xlab("") +
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")

conno_Spmonthly <- 
  ggplot(data = connodisc %>% filter(between(month_nu, 3, 5))) + 
  stat_summary(aes(x = year_nu, y = discharge_cufs), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9, color = "black") +
  stat_summary(aes(x = year_nu, y = discharge_cufs), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = 1.5, color = "black") +
  #geom_smooth(aes(x = year, y = dens.percomp.change, group = STATION, color = STATION), method = "lm", alpha = 0.5, size = .5) +
  theme_bw(base_size=20) + 
  ylab("Mean Spring discharge") + xlab("") +
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "")

conno_Wimonthly <- 
  ggplot(data = connodisc %>% filter(month_nu %in% c(1,2,10,11,12))) + 
  stat_summary(aes(x = year_nu, y = discharge_cufs, group = month_nu, color = month_nu), fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), size = .9) +
  stat_summary(aes(x = year_nu, y = discharge_cufs), fun.data = mean_se, geom = "line", fun.args = list(mult = 1), size = 1.5, color = "black") +
  #geom_smooth(aes(x = year, y = dens.percomp.change, group = STATION, color = STATION), method = "lm", alpha = 0.5, size = .5) +
  theme_bw(base_size=20) + 
  ylab("Mean Winter discharge") + xlab("") +
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")








           
  