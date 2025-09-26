#attemption crawl model
packages <- c("readr","tidyverse","lubridate","tidyr","dplyr","data.table","emmeans","ggplot2","ggspatial","ggthemes","ggeffects","ggpubr","egg",
              "postGIStools","RPostgreSQL","rgdal","rgeos","sp","raster","sf","crawl","forecast","moveHMM",
              "ctmm","ctmmweb","moveVis","move","momentuHMM","mgcv","mgcViz","itsadug","visreg","gratia")

###############################################################
#######     install required packages     #####################
###############################################################
install.packages("ggpubr")
install.packages("egg")
install.packages("postGIStools") #must find in archive
install.packages("rgeos")#must find in archive
install.packages("raster")
install.packages("forecast")
install.packages("ctmm") 
install.packages("ctmmweb") #find in archive
install.packages("moveVis") #archive
install.packages("move")
install.packages("mgcViz")
install.packages("itsadug")
install.packages("visreg")
install.packages("gratia")
install.packages("momentuHMM")
library(momentuHMM)
#all that could be installed have been

lapply(packages, require, character.only = TRUE)

nrow(pos_mean_wels_5min)
#downloading cleaned data, filered to be within rimov polygon
pos_mean_wels_5min<- read_csv("./data/pos_mean_wels_filter.csv")

unique(pos_mean_wels_5min$fishid) #15 individuals


pos_mean_wels_5min$fishid <- as.factor(pos_mean_wels_5min$fishid)

#check for duplicate rows
duplicates<- duplicated(pos_mean_wels_5min)
pos_mean_wels_5min[duplicates,]

#time_range <- pos_mean_wels_5min$timestamp_5min       >= as.POSIXct("2017-08-21 10:55:00") & pos_mean_wels_5min$timestamp_5min       <= as.POSIXct("2017-08-21 11:05:00")
#subset_data <- pos_mean_wels_5min[time_range, ]
#View(subset_data)

#filter for only fish id, timestamp, and location
pos_mean_wels_5min_subset<- pos_mean_wels_5min[,c(1,2,6:7)]

#rename columns
colnames(pos_mean_wels_5min_subset)[1] <- "ID"
colnames(pos_mean_wels_5min_subset)[2] <- "timestamp"
colnames(pos_mean_wels_5min_subset)[3] <- "x"
colnames(pos_mean_wels_5min_subset)[4] <- "y"

#create dataframe from subsetted data
wels_raw <- data.frame(ID = pos_mean_wels_5min_subset$ID,
                            time = pos_mean_wels_5min_subset$timestamp,
                            x = pos_mean_wels_5min_subset$x,
                            y = pos_mean_wels_5min_subset$y)

#make ID a factor, and order rows based on ID & time
wels_raw$ID <- as.numeric(as.factor(wels_raw$ID))
wels_raw <- wels_raw[with(wels_raw, order(ID, time)), ]

# total/column number of NA's
sum(is.na(wels_raw))
colSums(is.na(wels_raw))


#############################################################################################################
########################## Calculation of tracks based on gaps in the time series ###########################
#############################################################################################################



# Set thresholds: for gaps of what length do we still propose imputation of missing values?
th_gaps <- 30 # in minutes

# Tracks of less than what duration do we throw away?
th_dur <- 60 # minutes
th_dur <- th_dur / 15 # minutes/15 (sampling frequency) (this value excludes tracks with less than 8 rows, but includes tracks with more than 9 rows (i.e. 2.30h)


# Identify all observations with gaps greater than the threshold:
tdiff <- diff(wels_raw$time, units = "mins")                                                   # Get time differences between consecutive observations
newtrack <- c(TRUE, tdiff > th_gaps)                                                            # TRUE whenever a new track starts, the first value is TRUE (first track)
newtrack[(c(0,tdiff)<0)] <- TRUE                                                                # Also TRUE when the time series of a new animal starts, because then we need a new track 
wels_raw$trackID <- cumsum(newtrack)                                                            # Cumulative sum (+1 every time a new track starts)


# Remove tracks of less than th_dur
tracklengths <- rle(wels_raw$trackID) # lengths of the tracks
short_tracks <- tracklengths$values[tracklengths$lengths < th_dur] # identify short tracks
wels_raw_tracks <- wels_raw[-which(wels_raw$trackID %in% short_tracks), ] # remove short tracks 





# Obtain a regular time series for each track:
# for the first track:
ts_start <- min(wels_raw_tracks$time[wels_raw_tracks$trackID == unique(wels_raw_tracks$trackID)[1]]) # starting time of first track
ts_end <- max(wels_raw_tracks$time[wels_raw_tracks$trackID == unique(wels_raw_tracks$trackID)[1]]) # end time of first track
ts_track <- seq.POSIXt(ts_start, ts_end, by = "5 min") # regular sequence for first track
ts_full <- ts_track
animalID_track <- rep(unique(wels_raw_tracks$ID[wels_raw_tracks$trackID == unique(wels_raw_tracks$trackID)[1]]), length(ts_track))
animalID_full<-animalID_track #MMF Addition

for (i in unique(wels_raw_tracks$trackID)[-1]) { # do it for the other tracks as well
  ts_start <- min(wels_raw_tracks$time[wels_raw_tracks$trackID == i])
  ts_end <- max(wels_raw_tracks$time[wels_raw_tracks$trackID == i])
  
  ts_track <- seq.POSIXt(ts_start, ts_end, by = "5 min")
  ts_full <- c(ts_full, ts_track)
  
  animalID_track <- rep(unique(wels_raw_tracks$ID[wels_raw_tracks$trackID == i]), length(ts_track))
  animalID_full <- c(animalID_full, animalID_track)
}
ts_full <- data.frame(time = ts_full, ID = animalID_full)

# Merge regular time series with data to create NA rows for all gaps shorter than the threshold
wels_raw_tracks <- full_join(ts_full, wels_raw_tracks, by = (c("time", "ID")))


# Let's check the frequency distribution of the track lengths after subtracting the short lengths
distr_lengths<-wels_raw_tracks %>% group_by(trackID) %>% summarise(length(time))
distr_lengths<-data.table(distr_lengths)
colnames(distr_lengths)[2]<- "lengths"
ggplot(distr_lengths, aes(x = lengths)) + geom_histogram()+xlim(0,50)


# total/column number of NA's
sum(is.na(wels_raw_tracks)) #36639
colSums(is.na(wels_raw_tracks))


##########################################################
######## Data interpolation with SSMs (Crawl) ############
######## Dataset Prep                         ############
##########################################################


# Replace NA values with previous non-NA value in track ID column
wels_raw_tracks_full <- for (i in 2:length(wels_raw_tracks$trackID)) {
  if (is.na(wels_raw_tracks$trackID[i])) {
    wels_raw_tracks$trackID[i] <- wels_raw_tracks$trackID[i-1]
  }
}


wels_raw_large_tracks <- wels_raw_tracks

unique(wels_raw_large_tracks$ID) # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15

nrow(wels_raw_large_tracks) #[1] 95780

#create separate track datasets for each individual
wels_raw_large_id_tracks_1 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "1",]
wels_raw_large_id_tracks_2 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "2",]
wels_raw_large_id_tracks_3 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "3",]
wels_raw_large_id_tracks_4 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "4",]
wels_raw_large_id_tracks_5 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "5",]
wels_raw_large_id_tracks_6 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "6",]
wels_raw_large_id_tracks_7 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "7",]
wels_raw_large_id_tracks_8 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "8",]
wels_raw_large_id_tracks_9 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "9",]
wels_raw_large_id_tracks_10 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "10",]
wels_raw_large_id_tracks_11 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "11",]
wels_raw_large_id_tracks_12 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "12",]
wels_raw_large_id_tracks_13 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "13",]
wels_raw_large_id_tracks_14 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "14",]
wels_raw_large_id_tracks_15 <- wels_raw_large_tracks[wels_raw_large_tracks$ID == "15",]

wels_raw_large_tracks_1 <- wels_raw_large_id_tracks_1
wels_raw_large_tracks_2 <- wels_raw_large_id_tracks_2
wels_raw_large_tracks_3 <- wels_raw_large_id_tracks_3
wels_raw_large_tracks_4 <- wels_raw_large_id_tracks_4
wels_raw_large_tracks_5 <- wels_raw_large_id_tracks_5
wels_raw_large_tracks_6 <- wels_raw_large_id_tracks_6
wels_raw_large_tracks_7 <- wels_raw_large_id_tracks_7
wels_raw_large_tracks_8 <- wels_raw_large_id_tracks_8
wels_raw_large_tracks_9 <- wels_raw_large_id_tracks_9
wels_raw_large_tracks_10 <- wels_raw_large_id_tracks_10
wels_raw_large_tracks_11 <- wels_raw_large_id_tracks_11
wels_raw_large_tracks_12 <- wels_raw_large_id_tracks_12
wels_raw_large_tracks_13 <- wels_raw_large_id_tracks_13
wels_raw_large_tracks_14 <- wels_raw_large_id_tracks_14
wels_raw_large_tracks_15 <- wels_raw_large_id_tracks_15


#rename ID column
wels_raw_large_tracks_1$ID <- NULL
colnames(wels_raw_large_tracks_1)[4] <- "ID"

wels_raw_large_tracks_2$ID <- NULL
colnames(wels_raw_large_tracks_2)[4] <- "ID"

wels_raw_large_tracks_3$ID <- NULL
colnames(wels_raw_large_tracks_3)[4] <- "ID"

wels_raw_large_tracks_4$ID <- NULL
colnames(wels_raw_large_tracks_4)[4] <- "ID"

wels_raw_large_tracks_5$ID <- NULL
colnames(wels_raw_large_tracks_5)[4] <- "ID"

wels_raw_large_tracks_6$ID <- NULL
colnames(wels_raw_large_tracks_6)[4] <- "ID"

wels_raw_large_tracks_7$ID <- NULL
colnames(wels_raw_large_tracks_7)[4] <- "ID"

wels_raw_large_tracks_8$ID <- NULL
colnames(wels_raw_large_tracks_8)[4] <- "ID"

wels_raw_large_tracks_9$ID <- NULL
colnames(wels_raw_large_tracks_9)[4] <- "ID"

wels_raw_large_tracks_10$ID <- NULL
colnames(wels_raw_large_tracks_10)[4] <- "ID"

wels_raw_large_tracks_11$ID <- NULL
colnames(wels_raw_large_tracks_11)[4] <- "ID"

wels_raw_large_tracks_12$ID <- NULL
colnames(wels_raw_large_tracks_12)[4] <- "ID"

wels_raw_large_tracks_13$ID <- NULL
colnames(wels_raw_large_tracks_13)[4] <- "ID"

wels_raw_large_tracks_14$ID <- NULL
colnames(wels_raw_large_tracks_14)[4] <- "ID"

wels_raw_large_tracks_15$ID <- NULL
colnames(wels_raw_large_tracks_15)[4] <- "ID"

# List of unique IDs
unique_ids <- unique(wels_raw_large_tracks$ID) #15

#summarize tracks from individual 1 across all columns
wels_raw_large_tracks_1 %>%
  summarise_all(list(n_distinct))






#################################################################
################## Running Crawl Model ##########################
################################################################### redownload data w/o day hour calculated ##detect duplicates --> try crawl model with pikeperch see if it works 


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw`_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_1 <- crawlWrap(obsData=wels_raw_large_tracks_1, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 650000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_1 <- wels_raw_large_tracks_crawl_1$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_1 <- as.data.frame(wels_raw_large_tracks_crawl.df_1)

colnames(wels_raw_large_tracks_crawl.df_1)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_1)[4] <- "y"


# It has not computed the track "21" so I'll add apart from the original dataset
wels_raw_large_tracks_1_21 <- wels_raw_large_tracks_1[wels_raw_large_tracks_1$ID == "21",]

wels_raw_large_tracks_crawl_1_21 <- crawlWrap(obsData=wels_raw_large_tracks_1_21, timeStep="5 min",
                                              fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_1_21<- wels_raw_large_tracks_crawl_1_21$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_1_21<- as.data.frame(wels_raw_large_tracks_crawl.df_1_21)

colnames(wels_raw_large_tracks_crawl.df_1_21)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_1_21)[4] <- "y"

nrow(wels_raw_large_tracks_crawl.df_1_21)
nrow(wels_raw_large_tracks_1_21)


# It has not computed the track "47" so I'll add apart from the original dataset
wels_raw_large_tracks_1_47 <- wels_raw_large_tracks_1[wels_raw_large_tracks_1$ID == "47",]

wels_raw_large_tracks_crawl_1_47 <- crawlWrap(obsData=wels_raw_large_tracks_1_47, timeStep="5 min",
                                              fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_1_47<- wels_raw_large_tracks_crawl_1_47$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_1_47<- as.data.frame(wels_raw_large_tracks_crawl.df_1_47)

colnames(wels_raw_large_tracks_crawl.df_1_47)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_1_47)[4] <- "y"

nrow(wels_raw_large_tracks_crawl.df_1_47)
nrow(wels_raw_large_tracks_1_47)


# It has not computed the track "84" so I'll add apart from the original dataset
wels_raw_large_tracks_1_84 <- wels_raw_large_tracks_1[wels_raw_large_tracks_1$ID == "84",]

wels_raw_large_tracks_crawl_1_84 <- crawlWrap(obsData=wels_raw_large_tracks_1_84, timeStep="5 min",
                                              fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_1_84<- wels_raw_large_tracks_crawl_1_84$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_1_84<- as.data.frame(wels_raw_large_tracks_crawl.df_1_84)

colnames(wels_raw_large_tracks_crawl.df_1_84)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_1_84)[4] <- "y"

nrow(wels_raw_large_tracks_crawl.df_1_84)
nrow(wels_raw_large_tracks_1_84)


# It has not computed the track "92" so I'll add apart from the original dataset
wels_raw_large_tracks_1_92 <- wels_raw_large_tracks_1[wels_raw_large_tracks_1$ID == "92",]

wels_raw_large_tracks_crawl_1_92 <- crawlWrap(obsData=wels_raw_large_tracks_1_92, timeStep="5 min",
                                              fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_1_92<- wels_raw_large_tracks_crawl_1_92$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_1_92<- as.data.frame(wels_raw_large_tracks_crawl.df_1_92)

colnames(wels_raw_large_tracks_crawl.df_1_92)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_1_92)[4] <- "y"

nrow(wels_raw_large_tracks_crawl.df_1_92)
nrow(wels_raw_large_tracks_1_92)



#bind original dataset with non-computed tracks
wels_raw_large_tracks_crawl.df_1 <- rbind(wels_raw_large_tracks_crawl.df_1_21,  wels_raw_large_tracks_crawl.df_1_47, wels_raw_large_tracks_crawl.df_1_84, wels_raw_large_tracks_crawl.df_1_92,  wels_raw_large_tracks_crawl.df_1)

wels_raw_large_tracks_crawl.df_1 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_1$ID,timestamp = wels_raw_large_tracks_crawl.df_1$time,
                                               location.long = wels_raw_large_tracks_crawl.df_1$x,location.lat = wels_raw_large_tracks_crawl.df_1$y)



# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_1 <- data.frame(ID = wels_raw_large_id_tracks_1$ID, individual.local.identifier = wels_raw_large_id_tracks_1$trackID)


wels_raw_large_tracks_crawl.df_1$others <- df_col_wels_1$ID[match(wels_raw_large_tracks_crawl.df_1$individual.local.identifier, df_col_wels_1$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_1$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_1)[4] <- "individual.local.identifier"


#compare numbers of tracks to crawl - should be = 
nrow(wels_raw_large_tracks_crawl.df_1) #[1] 8349
nrow(wels_raw_large_tracks_1) #[1] 8349
unique(wels_raw_large_tracks_crawl.df_1$individual.local.identifier) [1]



#############################################################################################################################################################
################################ Individual 2 ###################################
#############################################################################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_2 <- crawlWrap(obsData=wels_raw_large_tracks_2, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 500000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_2 <- wels_raw_large_tracks_crawl_2$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_2 <- as.data.frame(wels_raw_large_tracks_crawl.df_2)

colnames(wels_raw_large_tracks_crawl.df_2)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_2)[4] <- "y"


# It has not computed the track "268" so I'll add apart from the original dataset
wels_raw_large_tracks_2_268 <- wels_raw_large_tracks_2[wels_raw_large_tracks_2$ID == "268",]

wels_raw_large_tracks_crawl_2_268 <- crawlWrap(obsData=wels_raw_large_tracks_2_268, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_2_268<- wels_raw_large_tracks_crawl_2_268$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_2_268<- as.data.frame(wels_raw_large_tracks_crawl.df_2_268)

colnames(wels_raw_large_tracks_crawl.df_2_268)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_2_268)[4] <- "y"


nrow(wels_raw_large_tracks_crawl.df_2_268)
nrow(wels_raw_large_tracks_2_268)


# It has not computed the track "318" so I'll add apart from the original dataset
wels_raw_large_tracks_2_318 <- wels_raw_large_tracks_2[wels_raw_large_tracks_2$ID == "318",]

wels_raw_large_tracks_crawl_2_318 <- crawlWrap(obsData=wels_raw_large_tracks_2_318, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_2_318<- wels_raw_large_tracks_crawl_2_318$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_2_318<- as.data.frame(wels_raw_large_tracks_crawl.df_2_318)

colnames(wels_raw_large_tracks_crawl.df_2_318)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_2_318)[4] <- "y"


nrow(wels_raw_large_tracks_crawl.df_2_318)
nrow(wels_raw_large_tracks_2_318)


# It has not computed the track "410" so I'll add apart from the original dataset
wels_raw_large_tracks_2_410 <- wels_raw_large_tracks_2[wels_raw_large_tracks_2$ID == "410",]

wels_raw_large_tracks_crawl_2_410 <- crawlWrap(obsData=wels_raw_large_tracks_2_410, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_2_410<- wels_raw_large_tracks_crawl_2_410$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_2_410<- as.data.frame(wels_raw_large_tracks_crawl.df_2_410)

colnames(wels_raw_large_tracks_crawl.df_2_410)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_2_410)[4] <- "y"


nrow(wels_raw_large_tracks_crawl.df_2_410)
nrow(wels_raw_large_tracks_2_410)


# It has not computed the track "427" so I'll add apart from the original dataset
wels_raw_large_tracks_2_427 <- wels_raw_large_tracks_2[wels_raw_large_tracks_2$ID == "427",]

wels_raw_large_tracks_crawl_2_427 <- crawlWrap(obsData=wels_raw_large_tracks_2_427, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_2_427<- wels_raw_large_tracks_crawl_2_427$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_2_427<- as.data.frame(wels_raw_large_tracks_crawl.df_2_427)

colnames(wels_raw_large_tracks_crawl.df_2_427)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_2_427)[4] <- "y"


nrow(wels_raw_large_tracks_crawl.df_2_427)
nrow(wels_raw_large_tracks_2_427)


#combine original tracks + recomputed tracks
wels_raw_large_tracks_crawl.df_2 <- rbind(wels_raw_large_tracks_crawl.df_2_268, wels_raw_large_tracks_crawl.df_2_318, wels_raw_large_tracks_crawl.df_2_410, wels_raw_large_tracks_crawl.df_2_427,  wels_raw_large_tracks_crawl.df_2)

wels_raw_large_tracks_crawl.df_2 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_2$ID,timestamp = wels_raw_large_tracks_crawl.df_2$time,
                                               location.long = wels_raw_large_tracks_crawl.df_2$x,location.lat = wels_raw_large_tracks_crawl.df_2$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_2 <- data.frame(ID = wels_raw_large_id_tracks_2$ID, individual.local.identifier = wels_raw_large_id_tracks_2$trackID)


wels_raw_large_tracks_crawl.df_2$others <- df_col_wels_2$ID[match(wels_raw_large_tracks_crawl.df_2$individual.local.identifier, df_col_wels_2$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_2$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_2)[4] <- "individual.local.identifier"


#double check rows are = 
nrow(wels_raw_large_tracks_crawl.df_2) #[1] 10616
nrow(wels_raw_large_tracks_2) #[1] 10616
unique(wels_raw_large_tracks_crawl.df_2$individual.local.identifier)



########################################################################################################################################
################################ Individual 3 ###################################
########################################################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_3 <- crawlWrap(obsData=wels_raw_large_tracks_3, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")            # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_3 <- wels_raw_large_tracks_crawl_3$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_3 <- as.data.frame(wels_raw_large_tracks_crawl.df_3)


colnames(wels_raw_large_tracks_crawl.df_3)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_3)[4] <- "y"


# It has not computed the track "476" so I'll add apart from the original dataset
wels_raw_large_tracks_3_476 <- wels_raw_large_tracks_3[wels_raw_large_tracks_3$ID == "476",]

wels_raw_large_tracks_crawl_3_476 <- crawlWrap(obsData=wels_raw_large_tracks_3_476, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_3_476 <- wels_raw_large_tracks_crawl_3_476$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_3_476 <- as.data.frame(wels_raw_large_tracks_crawl.df_3_476)

colnames(wels_raw_large_tracks_crawl.df_3_476)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_3_476)[4] <- "y"


#combine original tracks + recomputed tracks
wels_raw_large_tracks_crawl.df_3<- rbind(wels_raw_large_tracks_crawl.df_3_476,   wels_raw_large_tracks_crawl.df_3)

wels_raw_large_tracks_crawl.df_3<- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_3$ID,timestamp = wels_raw_large_tracks_crawl.df_3$time,
                                              location.long = wels_raw_large_tracks_crawl.df_3$x,location.lat = wels_raw_large_tracks_crawl.df_3$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_3 <- data.frame(ID = wels_raw_large_id_tracks_3$ID, individual.local.identifier = wels_raw_large_id_tracks_3$trackID)


wels_raw_large_tracks_crawl.df_3$others <- df_col_wels_3$ID[match(wels_raw_large_tracks_crawl.df_3$individual.local.identifier, df_col_wels_3$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_3$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_3)[4] <- "individual.local.identifier"



nrow(wels_raw_large_tracks_crawl.df_3) #[1] 4364
nrow(wels_raw_large_tracks_3) #[1] 4364
unique(wels_raw_large_tracks_crawl.df_3$individual.local.identifier)


############################################################################################################
################################################# Individual 4 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_4 <- crawlWrap(obsData=wels_raw_large_tracks_4, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 155000, method="Nelder-Mead")                # Parameters to adjust (e.g., theta=c(6.855, -0.007))

# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_4 <- wels_raw_large_tracks_crawl_4$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_4 <- as.data.frame(wels_raw_large_tracks_crawl.df_4)

colnames(wels_raw_large_tracks_crawl.df_4)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_4)[4] <- "y"


# It has not computed the track "553" so I'll add apart from the original dataset
wels_raw_large_tracks_4_553 <- wels_raw_large_tracks_4[wels_raw_large_tracks_4$ID == "553",]

wels_raw_large_tracks_crawl_4_553 <- crawlWrap(obsData=wels_raw_large_tracks_4_553, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_4_553 <- wels_raw_large_tracks_crawl_4_553$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_4_553 <- as.data.frame(wels_raw_large_tracks_crawl.df_4_553)

colnames(wels_raw_large_tracks_crawl.df_4_553)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_4_553)[4] <- "y"


wels_raw_large_tracks_crawl.df_4 <- rbind(wels_raw_large_tracks_crawl.df_4_553, wels_raw_large_tracks_crawl.df_4)



wels_raw_large_tracks_crawl.df_4 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_4$ID,timestamp = wels_raw_large_tracks_crawl.df_4$time,
                                               location.long = wels_raw_large_tracks_crawl.df_4$x,location.lat = wels_raw_large_tracks_crawl.df_4$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_4 <- data.frame(ID = wels_raw_large_id_tracks_4$ID, individual.local.identifier = wels_raw_large_id_tracks_4$trackID)


wels_raw_large_tracks_crawl.df_4$others <- df_col_wels_4$ID[match(wels_raw_large_tracks_crawl.df_4$individual.local.identifier, df_col_wels_4$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_4$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_4)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_4) #7270
nrow(wels_raw_large_tracks_4) #7270
unique(wels_raw_large_tracks_crawl.df_4$individual.local.identifier)



############################################################################################################
################################################# Individual 5 #############################################
############################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_5 <- crawlWrap(obsData=wels_raw_large_tracks_5, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 155000, method="Nelder-Mead")                          # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_5 <- wels_raw_large_tracks_crawl_5$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_5 <- as.data.frame(wels_raw_large_tracks_crawl.df_5)

colnames(wels_raw_large_tracks_crawl.df_5)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_5)[4] <- "y"

#no track errors

wels_raw_large_tracks_crawl.df_5 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_5$ID,timestamp = wels_raw_large_tracks_crawl.df_5$time,
                                               location.long = wels_raw_large_tracks_crawl.df_5$x,location.lat = wels_raw_large_tracks_crawl.df_5$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_5 <- data.frame(ID = wels_raw_large_id_tracks_5$ID, individual.local.identifier = wels_raw_large_id_tracks_5$trackID)


wels_raw_large_tracks_crawl.df_5$others <- df_col_wels_5$ID[match(wels_raw_large_tracks_crawl.df_5$individual.local.identifier, df_col_wels_5$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_5$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_5)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_5) #11126
nrow(wels_raw_large_tracks_5) #11126
unique(wels_raw_large_tracks_crawl.df_5$individual.local.identifier)



############################################################################################################
################################################# Individual 6 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_6 <- crawlWrap(obsData=wels_raw_large_tracks_6, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 750000, method="Nelder-Mead")                                   # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_6 <- wels_raw_large_tracks_crawl_6$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_6 <- as.data.frame(wels_raw_large_tracks_crawl.df_6)

colnames(wels_raw_large_tracks_crawl.df_6)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_6)[4] <- "y"


# It has not computed the track "722" so I'll add apart from the original dataset
wels_raw_large_tracks_6_722 <- wels_raw_large_tracks_6[wels_raw_large_tracks_6$ID == "722",]

wels_raw_large_tracks_crawl_6_722 <- crawlWrap(obsData=wels_raw_large_tracks_6_722, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_6_722 <- wels_raw_large_tracks_crawl_6_722$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_6_722 <- as.data.frame(wels_raw_large_tracks_crawl.df_6_722)

colnames(wels_raw_large_tracks_crawl.df_6_722)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_6_722)[4] <- "y"


wels_raw_large_tracks_crawl.df_6 <- rbind(wels_raw_large_tracks_crawl.df_6_722 ,wels_raw_large_tracks_crawl.df_6)

wels_raw_large_tracks_crawl.df_6 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_6$ID,timestamp = wels_raw_large_tracks_crawl.df_6$time,
                                               location.long = wels_raw_large_tracks_crawl.df_6$x,location.lat = wels_raw_large_tracks_crawl.df_6$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_6 <- data.frame(ID = wels_raw_large_id_tracks_6$ID, individual.local.identifier = wels_raw_large_id_tracks_6$trackID)


wels_raw_large_tracks_crawl.df_6$others <- df_col_wels_6$ID[match(wels_raw_large_tracks_crawl.df_6$individual.local.identifier, df_col_wels_6$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_6$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_6)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_6) #2841
nrow(wels_raw_large_tracks_6) #2841
unique(wels_raw_large_tracks_crawl.df_6$individual.local.identifier)


############################################################################################################
################################################# Individual 7 #############################################
############################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_7 <- crawlWrap(obsData=wels_raw_large_tracks_7, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_7 <- wels_raw_large_tracks_crawl_7$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_7 <- as.data.frame(wels_raw_large_tracks_crawl.df_7)

colnames(wels_raw_large_tracks_crawl.df_7)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_7)[4] <- "y"



# It has not computed the track "808" so I'll add apart from the original dataset
wels_raw_large_tracks_7_808 <- wels_raw_large_tracks_7[wels_raw_large_tracks_7$ID == "808",]

wels_raw_large_tracks_crawl_7_808 <- crawlWrap(obsData=wels_raw_large_tracks_7_808, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_7_808 <- wels_raw_large_tracks_crawl_7_808$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_7_808 <- as.data.frame(wels_raw_large_tracks_crawl.df_7_808)

colnames(wels_raw_large_tracks_crawl.df_7_808)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_7_808)[4] <- "y"


# It has not computed the track "863" so I'll add apart from the original dataset
wels_raw_large_tracks_7_863 <- wels_raw_large_tracks_7[wels_raw_large_tracks_7$ID == "863",]

wels_raw_large_tracks_crawl_7_863 <- crawlWrap(obsData=wels_raw_large_tracks_7_863, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_7_863 <- wels_raw_large_tracks_crawl_7_863$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_7_863 <- as.data.frame(wels_raw_large_tracks_crawl.df_7_863)

colnames(wels_raw_large_tracks_crawl.df_7_863)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_7_863)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_7 <- rbind(wels_raw_large_tracks_crawl.df_7_808, wels_raw_large_tracks_crawl.df_7_863,  wels_raw_large_tracks_crawl.df_7)

wels_raw_large_tracks_crawl.df_7 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_7$ID,timestamp = wels_raw_large_tracks_crawl.df_7$time,
                                               location.long = wels_raw_large_tracks_crawl.df_7$x,location.lat = wels_raw_large_tracks_crawl.df_7$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_7 <- data.frame(ID = wels_raw_large_id_tracks_7$ID, individual.local.identifier = wels_raw_large_id_tracks_7$trackID)


wels_raw_large_tracks_crawl.df_7$others <- df_col_wels_7$ID[match(wels_raw_large_tracks_crawl.df_7$individual.local.identifier, df_col_wels_7$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_7$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_7)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_7) #11837
nrow(wels_raw_large_tracks_7) #11837
unique(wels_raw_large_tracks_crawl.df_7$individual.local.identifier)



############################################################################################################
################################################# Individual 8 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_8 <- crawlWrap(obsData=wels_raw_large_tracks_8, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_8 <- wels_raw_large_tracks_crawl_8$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_8 <- as.data.frame(wels_raw_large_tracks_crawl.df_8)

colnames(wels_raw_large_tracks_crawl.df_8)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_8)[4] <- "y"


# It has not computed the track "914" so I'll add apart from the original dataset
wels_raw_large_tracks_8_914 <- wels_raw_large_tracks_8[wels_raw_large_tracks_8$ID == "914",]

wels_raw_large_tracks_crawl_8_914 <- crawlWrap(obsData=wels_raw_large_tracks_8_914, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_8_914 <- wels_raw_large_tracks_crawl_8_914$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_8_914 <- as.data.frame(wels_raw_large_tracks_crawl.df_8_914)

colnames(wels_raw_large_tracks_crawl.df_8_914)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_8_914)[4] <- "y"


# It has not computed the track "947" so I'll add apart from the original dataset
wels_raw_large_tracks_8_947 <- wels_raw_large_tracks_8[wels_raw_large_tracks_8$ID == "947",]

wels_raw_large_tracks_crawl_8_947 <- crawlWrap(obsData=wels_raw_large_tracks_8_947, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_8_947 <- wels_raw_large_tracks_crawl_8_947$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_8_947 <- as.data.frame(wels_raw_large_tracks_crawl.df_8_947)

colnames(wels_raw_large_tracks_crawl.df_8_947)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_8_947)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_8 <- rbind(wels_raw_large_tracks_crawl.df_8_914, wels_raw_large_tracks_crawl.df_8_947,  wels_raw_large_tracks_crawl.df_8)

wels_raw_large_tracks_crawl.df_8 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_8$ID,timestamp = wels_raw_large_tracks_crawl.df_8$time,
                                               location.long = wels_raw_large_tracks_crawl.df_8$x,location.lat = wels_raw_large_tracks_crawl.df_8$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_8 <- data.frame(ID = wels_raw_large_id_tracks_8$ID, individual.local.identifier = wels_raw_large_id_tracks_8$trackID)


wels_raw_large_tracks_crawl.df_8$others <- df_col_wels_8$ID[match(wels_raw_large_tracks_crawl.df_8$individual.local.identifier, df_col_wels_8$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_8$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_8)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_8) #5449
nrow(wels_raw_large_tracks_8) #5449
unique(wels_raw_large_tracks_crawl.df_8$individual.local.identifier)


############################################################################################################
################################################# Individual 9 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_9 <- crawlWrap(obsData=wels_raw_large_tracks_9, timeStep="5 min",
                                           fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_9 <- wels_raw_large_tracks_crawl_9$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_9 <- as.data.frame(wels_raw_large_tracks_crawl.df_9)

colnames(wels_raw_large_tracks_crawl.df_9)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_9)[4] <- "y"


# It has not computed the track "1036" so I'll add apart from the original dataset
wels_raw_large_tracks_9_1036 <- wels_raw_large_tracks_9[wels_raw_large_tracks_9$ID == "1036",]

wels_raw_large_tracks_crawl_9_1036 <- crawlWrap(obsData=wels_raw_large_tracks_9_1036, timeStep="5 min",
                                                fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_9_1036 <- wels_raw_large_tracks_crawl_9_1036$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_9_1036 <- as.data.frame(wels_raw_large_tracks_crawl.df_9_1036)

colnames(wels_raw_large_tracks_crawl.df_9_1036)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_9_1036)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_9 <- rbind(wels_raw_large_tracks_crawl.df_9_1036, wels_raw_large_tracks_crawl.df_9)

wels_raw_large_tracks_crawl.df_9 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_9$ID,timestamp = wels_raw_large_tracks_crawl.df_9$time,
                                               location.long = wels_raw_large_tracks_crawl.df_9$x,location.lat = wels_raw_large_tracks_crawl.df_9$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_9 <- data.frame(ID = wels_raw_large_id_tracks_9$ID, individual.local.identifier = wels_raw_large_id_tracks_9$trackID)


wels_raw_large_tracks_crawl.df_9$others <- df_col_wels_9$ID[match(wels_raw_large_tracks_crawl.df_9$individual.local.identifier, df_col_wels_9$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_9$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_9)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_9) #1358
nrow(wels_raw_large_tracks_9) #1358
unique(wels_raw_large_tracks_crawl.df_9$individual.local.identifier)


############################################################################################################
################################################# Individual 10 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_10 <- crawlWrap(obsData=wels_raw_large_tracks_10, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_10 <- wels_raw_large_tracks_crawl_10$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_10 <- as.data.frame(wels_raw_large_tracks_crawl.df_10)

colnames(wels_raw_large_tracks_crawl.df_10)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_10)[4] <- "y"

#no tracks with errors


wels_raw_large_tracks_crawl.df_10 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_10$ID,timestamp = wels_raw_large_tracks_crawl.df_10$time,
                                                location.long = wels_raw_large_tracks_crawl.df_10$x,location.lat = wels_raw_large_tracks_crawl.df_10$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_10 <- data.frame(ID = wels_raw_large_id_tracks_10$ID, individual.local.identifier = wels_raw_large_id_tracks_10$trackID)


wels_raw_large_tracks_crawl.df_10$others <- df_col_wels_10$ID[match(wels_raw_large_tracks_crawl.df_10$individual.local.identifier, df_col_wels_10$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_10$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_10)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_10) #951
nrow(wels_raw_large_tracks_10) #951
unique(wels_raw_large_tracks_crawl.df_10$individual.local.identifier)


############################################################################################################
################################################# Individual 11 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_11 <- crawlWrap(obsData=wels_raw_large_tracks_11, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_11 <- wels_raw_large_tracks_crawl_11$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_11 <- as.data.frame(wels_raw_large_tracks_crawl.df_11)

colnames(wels_raw_large_tracks_crawl.df_11)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_11)[4] <- "y"

# It has not computed the track "1210" so I'll add apart from the original dataset
wels_raw_large_tracks_11_1210 <- wels_raw_large_tracks_11[wels_raw_large_tracks_11$ID == "1210",]

wels_raw_large_tracks_crawl_11_1210 <- crawlWrap(obsData=wels_raw_large_tracks_11_1210, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_11_1210 <- wels_raw_large_tracks_crawl_11_1210$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_11_1210 <- as.data.frame(wels_raw_large_tracks_crawl.df_11_1210)

colnames(wels_raw_large_tracks_crawl.df_11_1210)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_11_1210)[4] <- "y"


# It has not computed the track "1213" so I'll add apart from the original dataset
wels_raw_large_tracks_11_1213 <- wels_raw_large_tracks_11[wels_raw_large_tracks_11$ID == "1213",]

wels_raw_large_tracks_crawl_11_1213 <- crawlWrap(obsData=wels_raw_large_tracks_11_1213, timeStep="5 min",
                                               fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_11_1213 <- wels_raw_large_tracks_crawl_11_1213$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_11_1213 <- as.data.frame(wels_raw_large_tracks_crawl.df_11_1213)

colnames(wels_raw_large_tracks_crawl.df_11_1213)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_11_1213)[4] <- "y"


# It has not computed the track "1252" so I'll add apart from the original dataset
wels_raw_large_tracks_11_1252 <- wels_raw_large_tracks_11[wels_raw_large_tracks_11$ID == "1252",]

wels_raw_large_tracks_crawl_11_1252 <- crawlWrap(obsData=wels_raw_large_tracks_11_1252, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_11_1252 <- wels_raw_large_tracks_crawl_11_1252$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_11_1252 <- as.data.frame(wels_raw_large_tracks_crawl.df_11_1252)

colnames(wels_raw_large_tracks_crawl.df_11_1252)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_11_1252)[4] <- "y"



# It has not computed the track "1263" so I'll add apart from the original dataset
wels_raw_large_tracks_11_1263 <- wels_raw_large_tracks_11[wels_raw_large_tracks_11$ID == "1263",]

wels_raw_large_tracks_crawl_11_1263 <- crawlWrap(obsData=wels_raw_large_tracks_11_1263, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_11_1263 <- wels_raw_large_tracks_crawl_11_1263$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_11_1263 <- as.data.frame(wels_raw_large_tracks_crawl.df_11_1263)

colnames(wels_raw_large_tracks_crawl.df_11_1263)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_11_1263)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_11 <- rbind(wels_raw_large_tracks_crawl.df_11_1210, wels_raw_large_tracks_crawl.df_11_1213, wels_raw_large_tracks_crawl.df_11_1252, wels_raw_large_tracks_crawl.df_11_1263, wels_raw_large_tracks_crawl.df_11)

wels_raw_large_tracks_crawl.df_11 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_11$ID,timestamp = wels_raw_large_tracks_crawl.df_11$time,
                                                location.long = wels_raw_large_tracks_crawl.df_11$x,location.lat = wels_raw_large_tracks_crawl.df_11$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_11 <- data.frame(ID = wels_raw_large_id_tracks_11$ID, individual.local.identifier = wels_raw_large_id_tracks_11$trackID)


wels_raw_large_tracks_crawl.df_11$others <- df_col_wels_11$ID[match(wels_raw_large_tracks_crawl.df_11$individual.local.identifier, df_col_wels_11$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_11$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_11)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_11) #3568
nrow(wels_raw_large_tracks_11) #3568
unique(wels_raw_large_tracks_crawl.df_11$individual.local.identifier)



############################################################################################################
################################################# Individual 12 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_12 <- crawlWrap(obsData=wels_raw_large_tracks_12, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_12 <- wels_raw_large_tracks_crawl_12$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_12 <- as.data.frame(wels_raw_large_tracks_crawl.df_12)

colnames(wels_raw_large_tracks_crawl.df_12)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_12)[4] <- "y"


# It has not computed the track "1463" so I'll add apart from the original dataset
wels_raw_large_tracks_12_1463 <- wels_raw_large_tracks_12[wels_raw_large_tracks_12$ID == "1463",]

wels_raw_large_tracks_crawl_12_1463 <- crawlWrap(obsData=wels_raw_large_tracks_12_1463, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_12_1463 <- wels_raw_large_tracks_crawl_12_1463$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_12_1463 <- as.data.frame(wels_raw_large_tracks_crawl.df_12_1463)

colnames(wels_raw_large_tracks_crawl.df_12_1463)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_12_1463)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_12 <- rbind(wels_raw_large_tracks_crawl.df_12_1463, wels_raw_large_tracks_crawl.df_12)

wels_raw_large_tracks_crawl.df_12 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_12$ID,timestamp = wels_raw_large_tracks_crawl.df_12$time,
                                                location.long = wels_raw_large_tracks_crawl.df_12$x,location.lat = wels_raw_large_tracks_crawl.df_12$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_12 <- data.frame(ID = wels_raw_large_id_tracks_12$ID, individual.local.identifier = wels_raw_large_id_tracks_12$trackID)


wels_raw_large_tracks_crawl.df_12$others <- df_col_wels_12$ID[match(wels_raw_large_tracks_crawl.df_12$individual.local.identifier, df_col_wels_12$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_12$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_12)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_12) #11199
nrow(wels_raw_large_tracks_12) #11199
unique(wels_raw_large_tracks_crawl.df_12$individual.local.identifier)



############################################################################################################
################################################# Individual 13 #############################################
############################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_13 <- crawlWrap(obsData=wels_raw_large_tracks_13, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_13 <- wels_raw_large_tracks_crawl_13$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_13 <- as.data.frame(wels_raw_large_tracks_crawl.df_13)

colnames(wels_raw_large_tracks_crawl.df_13)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_13)[4] <- "y"


# It has not computed the track "1519" so I'll add apart from the original dataset
wels_raw_large_tracks_13_1519 <- wels_raw_large_tracks_13[wels_raw_large_tracks_13$ID == "1519",]

wels_raw_large_tracks_crawl_13_1519 <- crawlWrap(obsData=wels_raw_large_tracks_13_1519, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_13_1519 <- wels_raw_large_tracks_crawl_13_1519$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_13_1519 <- as.data.frame(wels_raw_large_tracks_crawl.df_13_1519)

colnames(wels_raw_large_tracks_crawl.df_13_1519)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_13_1519)[4] <- "y"


# It has not computed the track "1545" so I'll add apart from the original dataset
wels_raw_large_tracks_13_1545 <- wels_raw_large_tracks_13[wels_raw_large_tracks_13$ID == "1545",]

wels_raw_large_tracks_crawl_13_1545 <- crawlWrap(obsData=wels_raw_large_tracks_13_1545, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_13_1545 <- wels_raw_large_tracks_crawl_13_1545$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_13_1545 <- as.data.frame(wels_raw_large_tracks_crawl.df_13_1545)

colnames(wels_raw_large_tracks_crawl.df_13_1545)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_13_1545)[4] <- "y"


# It has not computed the track "1605" so I'll add apart from the original dataset
wels_raw_large_tracks_13_1605 <- wels_raw_large_tracks_13[wels_raw_large_tracks_13$ID == "1605",]

wels_raw_large_tracks_crawl_13_1605 <- crawlWrap(obsData=wels_raw_large_tracks_13_1605, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_13_1605 <- wels_raw_large_tracks_crawl_13_1605$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_13_1605 <- as.data.frame(wels_raw_large_tracks_crawl.df_13_1605)

colnames(wels_raw_large_tracks_crawl.df_13_1605)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_13_1605)[4] <- "y"



# It has not computed the track "1686" so I'll add apart from the original dataset
wels_raw_large_tracks_13_1686 <- wels_raw_large_tracks_13[wels_raw_large_tracks_13$ID == "1686",]

wels_raw_large_tracks_crawl_13_1686 <- crawlWrap(obsData=wels_raw_large_tracks_13_1686, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_13_1686 <- wels_raw_large_tracks_crawl_13_1686$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_13_1686 <- as.data.frame(wels_raw_large_tracks_crawl.df_13_1686)

colnames(wels_raw_large_tracks_crawl.df_13_1686)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_13_1686)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_13 <- rbind(wels_raw_large_tracks_crawl.df_13_1519, wels_raw_large_tracks_crawl.df_13_1545, wels_raw_large_tracks_crawl.df_13_1605,  wels_raw_large_tracks_crawl.df_13_1686,  wels_raw_large_tracks_crawl.df_13)

wels_raw_large_tracks_crawl.df_13 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_13$ID,timestamp = wels_raw_large_tracks_crawl.df_13$time,
                                                location.long = wels_raw_large_tracks_crawl.df_13$x,location.lat = wels_raw_large_tracks_crawl.df_13$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_13 <- data.frame(ID = wels_raw_large_id_tracks_13$ID, individual.local.identifier = wels_raw_large_id_tracks_13$trackID)


wels_raw_large_tracks_crawl.df_13$others <- df_col_wels_13$ID[match(wels_raw_large_tracks_crawl.df_13$individual.local.identifier, df_col_wels_13$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_13$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_13)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_13) #4047
nrow(wels_raw_large_tracks_13) #4047
unique(wels_raw_large_tracks_crawl.df_13$individual.local.identifier)



############################################################################################################
################################################# Individual 14 #############################################
############################################################################################################


# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_14 <- crawlWrap(obsData=wels_raw_large_tracks_14, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_14 <- wels_raw_large_tracks_crawl_14$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_14 <- as.data.frame(wels_raw_large_tracks_crawl.df_14)

colnames(wels_raw_large_tracks_crawl.df_14)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_14)[4] <- "y"


# It has not computed the track "1769" so I'll add apart from the original dataset
wels_raw_large_tracks_14_1769 <- wels_raw_large_tracks_14[wels_raw_large_tracks_14$ID == "1769",]

wels_raw_large_tracks_crawl_14_1769 <- crawlWrap(obsData=wels_raw_large_tracks_14_1769, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_14_1769 <- wels_raw_large_tracks_crawl_14_1769$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_14_1769 <- as.data.frame(wels_raw_large_tracks_crawl.df_14_1769)

colnames(wels_raw_large_tracks_crawl.df_14_1769)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_14_1769)[4] <- "y"


# It has not computed the track "1777" so I'll add apart from the original dataset
wels_raw_large_tracks_14_1777 <- wels_raw_large_tracks_14[wels_raw_large_tracks_14$ID == "1777",]

wels_raw_large_tracks_crawl_14_1777 <- crawlWrap(obsData=wels_raw_large_tracks_14_1777, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_14_1777 <- wels_raw_large_tracks_crawl_14_1777$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_14_1777 <- as.data.frame(wels_raw_large_tracks_crawl.df_14_1777)

colnames(wels_raw_large_tracks_crawl.df_14_1777)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_14_1777)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_14 <- rbind(wels_raw_large_tracks_crawl.df_14_1769, wels_raw_large_tracks_crawl.df_14_1777,  wels_raw_large_tracks_crawl.df_14)

wels_raw_large_tracks_crawl.df_14 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_14$ID,timestamp = wels_raw_large_tracks_crawl.df_14$time,
                                                location.long = wels_raw_large_tracks_crawl.df_14$x,location.lat = wels_raw_large_tracks_crawl.df_14$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_14 <- data.frame(ID = wels_raw_large_id_tracks_14$ID, individual.local.identifier = wels_raw_large_id_tracks_14$trackID)


wels_raw_large_tracks_crawl.df_14$others <- df_col_wels_14$ID[match(wels_raw_large_tracks_crawl.df_14$individual.local.identifier, df_col_wels_14$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_14$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_14)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_14) #3288
nrow(wels_raw_large_tracks_14) #3288
unique(wels_raw_large_tracks_crawl.df_14$individual.local.identifier)



############################################################################################################
################################################# Individual 15 #############################################
############################################################################################################

# Fit crawl (CTCRW) model and predict locations every 5-min
#colnames(wels_raw_tracks)[2] <- "time"
wels_raw_large_tracks_crawl_15 <- crawlWrap(obsData=wels_raw_large_tracks_15, timeStep="5 min",
                                            fixPar=c(NA,NA), attempts = 1550000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(6.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_15 <- wels_raw_large_tracks_crawl_15$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_15 <- as.data.frame(wels_raw_large_tracks_crawl.df_15)

colnames(wels_raw_large_tracks_crawl.df_15)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_15)[4] <- "y"


# It has not computed the track "1916" so I'll add apart from the original dataset
wels_raw_large_tracks_15_1916 <- wels_raw_large_tracks_15[wels_raw_large_tracks_15$ID == "1916",]

wels_raw_large_tracks_crawl_15_1916 <- crawlWrap(obsData=wels_raw_large_tracks_15_1916, timeStep="5 min",
                                                 fixPar=c(NA,NA), attempts = 40000000, method="Nelder-Mead")                                  # Parameters to adjust (e.g., theta=c(7.855, -0.007))


# Selection of interpolated locations and renaming of variables
wels_raw_large_tracks_crawl.df_15_1916 <- wels_raw_large_tracks_crawl_15_1916$crwPredict[, c("time", "ID", "mu.x", "mu.y")]
wels_raw_large_tracks_crawl.df_15_1916 <- as.data.frame(wels_raw_large_tracks_crawl.df_15_1916)

colnames(wels_raw_large_tracks_crawl.df_15_1916)[3] <- "x"
colnames(wels_raw_large_tracks_crawl.df_15_1916)[4] <- "y"


#combine original track + recalculated tracks
wels_raw_large_tracks_crawl.df_15 <- rbind(wels_raw_large_tracks_crawl.df_15_1916,  wels_raw_large_tracks_crawl.df_15)

wels_raw_large_tracks_crawl.df_15 <- data.frame(individual.local.identifier = wels_raw_large_tracks_crawl.df_15$ID,timestamp = wels_raw_large_tracks_crawl.df_15$time,
                                                location.long = wels_raw_large_tracks_crawl.df_15$x,location.lat = wels_raw_large_tracks_crawl.df_15$y)


# Create a new dataframe with the extracted ID and tracklD column variables
df_col_wels_15 <- data.frame(ID = wels_raw_large_id_tracks_15$ID, individual.local.identifier = wels_raw_large_id_tracks_15$trackID)


wels_raw_large_tracks_crawl.df_15$others <- df_col_wels_15$ID[match(wels_raw_large_tracks_crawl.df_15$individual.local.identifier, df_col_wels_15$individual.local.identifier)]

wels_raw_large_tracks_crawl.df_15$individual.local.identifier <- NULL
colnames(wels_raw_large_tracks_crawl.df_15)[4] <- "individual.local.identifier"


nrow(wels_raw_large_tracks_crawl.df_15) #9140
nrow(wels_raw_large_tracks_15) #9140
unique(wels_raw_large_tracks_crawl.df_15$individual.local.identifier)


#####################################################################################

#combine all individual tracks together
wels_raw_large_tracks_crawl.df <- rbind(wels_raw_large_tracks_crawl.df_1, wels_raw_large_tracks_crawl.df_2, wels_raw_large_tracks_crawl.df_3,
                                        wels_raw_large_tracks_crawl.df_4, wels_raw_large_tracks_crawl.df_5, wels_raw_large_tracks_crawl.df_6,
                                        wels_raw_large_tracks_crawl.df_7, wels_raw_large_tracks_crawl.df_8, wels_raw_large_tracks_crawl.df_9, wels_raw_large_tracks_crawl.df_10, wels_raw_large_tracks_crawl.df_11, wels_raw_large_tracks_crawl.df_12, wels_raw_large_tracks_crawl.df_13, wels_raw_large_tracks_crawl.df_14, wels_raw_large_tracks_crawl.df_15)



# Reassign the levels of individual.local.identifier
wels_raw_large_tracks_crawl.df$individual.local.identifier <- factor(
  wels_raw_large_tracks_crawl.df$individual.local.identifier,
  levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
  labels = c("T449268_1", "T449269_1", "T449270_1", "T449274_1", "T449275_1", "T449276_1", "T449278_1", "T449280_1", "T449282_1", "T449284_1", "T449286_1", "T449288_1", "T449314_1", "T449318_1", "T449319_1")
)


# Rename the columns
colnames(wels_raw_large_tracks_crawl.df) <- c("time", "x", "y", "ID")

# Display the head of the updated dataset to verify changes
head(wels_raw_large_tracks_crawl.df)


write.csv(wels_raw_large_tracks_crawl.df, "wels_raw_large_tracks_crawl.df.csv")

wels_raw_large_tracks_crawl.df <- data.table(read_csv("./wels_raw_large_tracks_crawl.df.csv"))



