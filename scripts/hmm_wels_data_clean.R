#performing HMM model with wels catfish calculated tracks
source("./scripts/libraries.R")

#read in track data which has been processed with crawl model
HMM<-read_csv("./data/wels_raw_large_tracks_crawl.df.csv")
HMM$ID<- as.factor(HMM$ID)

#calculate step lengths and turning angles
HMM_wels<-moveHMM::prepData(HMM, type = "UTM", coordNames = c("x","y"))

#plot all tracks, step lengths + turning angles
plot(HMM_wels, compact = T)

#investigate outliers
boxplot(HMM_wels$step,
        ylab = "step length")

#difficult to see where cut off of outliers is. lets subset to see the outliers more clearly
HMM_wels$week <- week(HMM_wels$time)

HMM_wels_subset_wk29 <- HMM_wels %>%
  filter(week == 29)

range(HMM_wels_subset_wk29$step)

boxplot(HMM_wels_subset_wk29$step,
        ylab = "step length")

#outliers noticed in initial plotting of step length. lets remove them

# step length
step_length <- 200

#filtering out steps larger than 200m (how far is physically possible for a fish of X size to swim [Body length / sec]) 
HMM_wels_clean <- HMM_wels %>%
  filter(step < step_length) #HMM does not work with NA's, must remove them if they are present

range(HMM_wels_clean$step)
#had to remove calculated step lengths so can convert it back into a moveObject in other scripts
HMM_wels_clean<- HMM_wels_clean[,-c(2:3)]

which(is.na(HMM_wels_clean)) #verify no NA's left

#save cleaned HMM data
write.csv(HMM_wels_clean, "./data/HMM_wels_clean.csv")
