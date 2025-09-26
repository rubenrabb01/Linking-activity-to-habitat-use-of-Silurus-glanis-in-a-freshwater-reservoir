#HMM Run - Model Refinement
################################################################################
##################   Running 3 State Models With ALL Data   ###################        
################################################################################

source("./scripts/hmm_function.R")
source("./scripts/libraries.R")

#download + filter data

#read in track data which has been processed with crawl model
HMM<-read_csv("./data/wels_raw_large_tracks_crawl.df.csv")
HMM$week <- week(HMM$time)

#convert timestamp to day_hr, and to seconds
HMM <- HMM %>% 
  mutate(day_time = hour(time)+(minute(time)/60)) 

HMM_2<- HMM %>% 
  mutate(hour = hour(time),
         minute = minute(time), 
         sec = second(time)) 

HMM_2<- HMM_2 %>% 
  mutate(day_time_sec = hour*60*60 + minute*60 + sec)

#range(HMM_2$day_time_sec)

#calculate cosine time variables to make time cyclical
HMM_2<- HMM_2 %>%
  mutate(cos_time = cos(2*pi*day_time_sec/86400))

HMM_2<- HMM_2 %>%
  mutate(sin_time = sin(2*pi*day_time_sec/86400))


#####################################################################
######### Compare old vs new diel period counts #####################
#####################################################################

#remove weeks with insufficient and unequal amounts of data (less than 100 observations / diel period / day )
#create column for diel period --> OLD 
HMM_2_dp_Test<- HMM_2 %>% 
  mutate(diel_period = case_when(day_time >= 8  & day_time < 20 ~ "day",
                                 day_time >= 20 & day_time < 24 | day_time >= 0 & day_time < 8 ~ "night"))

HMM_counts1b<- HMM_2_dp_Test %>%
  group_by(ID, week, diel_period) %>% 
  tally() 

HMM_counts2b<- HMM_counts1b%>% 
  spread(diel_period, n)

HMM_counts3b<- HMM_counts2b %>% 
  rename(ncount_day = day,
         ncount_night = night)

print(HMM_counts3b, n =120)



#update the diel period classification to be based on the sunrise and sunset timing 
getNightDay <- function(x,
                        lat = 48.8497428,
                        lon = 14.4903039,
                        xtz = "UTC",
                        returnSunrSuns = FALSE,
                        day = "day",
                        night = "night") {
  # Check input
  if (length(x) == 0) stop("length of input is 0")
  if (all(is.na(x))) return(as.character(x))
  
  # Ensure required packages are loaded
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("suncalc", quietly = TRUE)) {
    stop("Package 'suncalc' is required but not installed.")
  }
  library(data.table)
  library(suncalc)
  
  # Convert x to POSIXct using the provided time zone
  x <- as.POSIXct(x, tz = xtz)
  
  # Create a date sequence that covers the range of x (with an extra day on each side)
  from_date <- as.Date(min(x, na.rm = TRUE), tz = xtz) - 1
  to_date   <- as.Date(max(x, na.rm = TRUE), tz = xtz) + 1
  day_seq   <- seq.Date(from_date, to_date, by = "day")
  
  # Get sunrise and sunset times for each day
  sunTimes <- data.table::as.data.table(
    suncalc::getSunlightTimes(date = day_seq, lat = lat, lon = lon, tz = xtz)
  )
  
  # Determine which columns contain sunrise/sunset times.
  # Some versions of suncalc return "sunrise" and "sunset" while others return "sunriseTime" and "sunsetTime"
  sunrise_col <- if ("sunrise" %in% names(sunTimes)) "sunrise" else if ("sunriseTime" %in% names(sunTimes)) "sunriseTime" else stop("No sunrise column found")
  sunset_col  <- if ("sunset"  %in% names(sunTimes)) "sunset"  else if ("sunsetTime"  %in% names(sunTimes)) "sunsetTime"  else stop("No sunset column found")
  
  # Create a Date column for merging
  sunTimes[, Date := as.Date(date, tz = xtz)]
  
  # Subset and rename the sunrise and sunset columns in one step
  sunTimes_subset <- sunTimes[, .(Date,
                                  sunrise = get(sunrise_col),
                                  sunset  = get(sunset_col))]
  # Create a data.table for the input times
  dt <- data.table(x = x,
                   xorder = seq_along(x),
                   Date = as.Date(x, tz = xtz))
  
  # Merge input times with the sunrise/sunset times by Date
  dt <- merge(dt, sunTimes_subset, by = "Date", all.x = TRUE)
  
  # Determine the diel period: assign "day" if x is between sunrise and sunset, otherwise "night"
  dt[, Diel.period := ifelse(x > sunrise & x < sunset, day, night)]
  
  # Restore the original order
  setorder(dt, xorder)
  
  # Return either just the diel period or additional sunrise/sunset info
  if (returnSunrSuns) {
    return(dt[, .(Diel.period, sunrise, sunset)])
  } else {
    return(dt$Diel.period)
  }
}


# Start by making a copy and converting to data.table
HMM_2_dp_Test2 <- data.table(HMM_2_dp_Test)

# Create timestamp column from 'time'
HMM_2_dp_Test2[, timestamp := time]

# Apply getNightDay function
HMM_2_dp_Test2[, dp := getNightDay(timestamp)]

#compare to previous diel period counts
HMM_counts1<- HMM_2_dp_Test2 %>%
  group_by(ID, week, dp) %>% 
  tally() 

HMM_counts2<- HMM_counts1%>% 
  spread(dp, n)

HMM_counts2<- HMM_counts2 %>% 
  rename(ncount_day = day,
         ncount_night = night)

print(HMM_counts2, n =120)


HMM_counts<- full_join(HMM_counts1,HMM_counts2, by = c("ID", "week"))

HMM_full<- full_join(HMM_2_dp_Test2,HMM_counts, by = c("ID", "week", "dp"))

#transform days with NA count to 0 
HMM_full$ncount_day[is.na(HMM_full$ncount_day)] <- 0
HMM_full$ncount_night[is.na(HMM_full$ncount_night)] <- 0
range(HMM_full$ncount_day)

#filter out 
HMM_full2<- HMM_full %>% 
  group_by(ID, week) %>% 
  filter(ncount_day >= 100 & ncount_night >= 100)

#prep data for HMM
HMM_prep_copy<- moveHMM::prepData(HMM_full2, type = "UTM", coordNames = c("x","y"))

#filter out step lengths greater than 200m (physiologically difficult for fish to move faster than bodylength / sec)
step_length <- 200

HMM_prep_clean_copy <- HMM_prep_copy %>%
  filter(step < step_length) %>% 
  drop_na(angle)# HMM does not work with NA's, must remove them if they are present (usually just in 1st column of angle)


################################################################################
#####   We know that a 3 state model is preferred for wels, so let's try  ######
#####      models with varying covariates. Mainly: Time & Temperature     ######
################################################################################


################################################################################
####################### 3 state model + time covar #############################
################################################################################

stepMean_3s <- c(3,15,50) #means of gamma distribution (step length)
stepSD0_3s <- c(3,15,50) #std. dev of gamma distribution (step length)
#zeromass0_3s <- c(0.01,0.001) #step zero-mass
stepPar0_3s <- c(stepMean_3s,stepSD0_3s)#, zeromass0_3s)
angleMean_3s <- c(pi,0,pi/2) #means of von Mises dist (turning angles)
angleCon0_3s <- c(1,2,5) #precision of von Mises dist (turning angles)
anglePar0_3s <- c(angleMean_3s,angleCon0_3s)

#test HMM with time as 1 circular variable rather than 2 separate vars
HMM_prep_clean_copy_test <- HMM_prep_clean_copy
HMM_prep_clean_copy_test <- HMM_prep_clean_copy_test %>% 
  mutate(circ_time = atan2(cos_time,sin_time))

#HMM_wels_fit_copy_test1<- moveHMM::fitHMM(HMM_prep_clean_copy_test, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s, formula = ~circ_time) #, stepDist = "gamma", angleDist = "vm"

HMM_wels_fit_copy_time <- moveHMM::fitHMM(HMM_prep_clean_copy, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s, formula = ~cos_time + sin_time) #, stepDist = "gamma", angleDist = "vm"
HMM_wels_fit_copy_time <- HMM_wels_fit_copy


#view model output
plot(HMM_wels_fit_time, compact = T)
plot(HMM_wels_fit_test1, compact = T)

#AIC(HMM_wels_fit_time,HMM_wels_fit_test1,HMM_wels_fit_nocovars)#time as 1 variable is not preferable (based on AIC & log likelihood) to time as 2 vars (cos and sin)

################################################################################
####################    3 state model with no covars   #########################
################################################################################

HMM_wels_fit_copy_nocovars <- moveHMM::fitHMM(HMM_prep_clean_copy, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s) #, stepDist = "gamma", angleDist = "vm"
plot(HMM_wels_fit_nocovars, compact = T)

################################################################################
###################### 3 state model + temp covar ##############################
################################################################################
#download mean epilimnion data 
epi_mean_temp<- read_csv("./data/mean_temp_6m.csv")
epi_mean_temp<- epi_mean_temp %>% 
  rename(time = hd_timestamp_utc) %>% 
  filter(time > "2017-07-01 UTC" & time < "2017-09-01 UTC")

range(epi_mean_temp$time)

#merge will cleaned dataset
HMM_full2_copy <- left_join(HMM_full2,epi_mean_temp, by = c("time"))

#standardize covariates
#HMM_full2_copy2 <- HMM_full2_copy %>% 
#  filter(!is.na(temp_mean)) %>% 
#  mutate(stand_mean_temp = (temp_mean - mean(temp_mean))/sd(temp_mean))

#Define mu and sigma first
mu <- mean(HMM_full2_copy$temp_mean, na.rm = TRUE)
sigma <- sd(HMM_full2_copy$temp_mean, na.rm = TRUE)

# Use them to standardize
HMM_full2_copy2 <- HMM_full2_copy %>%
  filter(!is.na(temp_mean)) %>%
  mutate(stand_mean_temp = (temp_mean - mu) / sigma)

#lowT_original <- min(HMM_full2_copy2$stand_mean_temp) * sigma + mu
#highT_original <- max(HMM_full2_copy2$stand_mean_temp) * sigma + mu



#prep data for hmm
HMM_full2_copy_prep2b<- moveHMM::prepData(HMM_full2_copy2, type = "UTM", coordNames = c("x","y"))
 
#filter out step lengths greater than 200m
step_length <- 200

HMM_full2_copy_prep2 <- HMM_full2_copy_prep2b %>%
  filter(step < step_length) %>% 
  drop_na(angle)# HMM does not work with NA's, must remove them if they are present


stepMean_3s <- c(3,15,50) #means of gamma distribution (step length)
stepSD0_3s <- c(3,15,50) #std. dev of gamma distribution (step length)
#zeromass0_3s <- c(0.01,0.001) #step zero-mass
stepPar0_3s <- c(stepMean_3s,stepSD0_3s)#, zeromass0_3s)
angleMean_3s <- c(pi,0,pi/2) #means of von Mises dist (turning angles)
angleCon0_3s <- c(1,2,5) #precision of von Mises dist (turning angles)
anglePar0_3s <- c(angleMean_3s,angleCon0_3s)

HMM_wels_fit_copy_temp<- moveHMM::fitHMM(HMM_full2_copy_prep2, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s, formula = ~stand_mean_temp) #, stepDist = "gamma", angleDist = "vm"

#plotStationary & plotStates are unable to consider cosine + sine together
plotStationary(HMM_wels_fit_temp, plotCI=TRUE)
plot(HMM_wels_fit_temp, compact = T)

dev.off()


################################################################################
#################### 3 state model + time + temp covar #########################
################################################################################

#run model with 3 covars --> had to change initial parms 
stepMean_3s <- c(3,25,75) #means of gamma distribution (step length)
stepSD0_3s <- c(3,25,75) #std. dev of gamma distribution (step length)
#zeromass0_3s <- c(0.01,0.001) #step zero-mass
stepPar0_3s <- c(stepMean_3s,stepSD0_3s)#, zeromass0_3s)
angleMean_3s <- c(pi,0,pi/2) #means of von Mises dist (turning angles)
angleCon0_3s <- c(1,2,5) #precision of von Mises dist (turning angles)
anglePar0_3s <- c(angleMean_3s,angleCon0_3s)

HMM_wels_fit_copy_3covars<- moveHMM::fitHMM(HMM_full2_copy_prep2, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s, formula = ~cos_time + sin_time + stand_mean_temp) #, stepDist = "gamma", angleDist = "vm"
plot(HMM_wels_fit_3covars, compact = T)

################################################################################
#################### 3 state model + time*temp covar ###########################
################################################################################

stepMean_3s <- c(3,15,50) #means of gamma distribution (step length)
stepSD0_3s <- c(3,15,50) #std. dev of gamma distribution (step length)
#zeromass0_3s <- c(0.01,0.001) #step zero-mass
stepPar0_3s <- c(stepMean_3s,stepSD0_3s)#, zeromass0_3s)
angleMean_3s <- c(pi,0,pi/2) #means of von Mises dist (turning angles)
angleCon0_3s <- c(1,2,5) #precision of von Mises dist (turning angles)
anglePar0_3s <- c(angleMean_3s,angleCon0_3s)

HMM_wels_fit_copy_int<- moveHMM::fitHMM(HMM_full2_copy_prep2, nbStates = 3, stepPar0 = stepPar0_3s, anglePar0 = anglePar0_3s, formula = ~(cos_time + sin_time)*stand_mean_temp) #, stepDist = "gamma", angleDist = "vm"

#-436390.5 , state 1/2/3 steps = 2.34/71.32/21.06
#state 1/2/3 angle = 3.14/0.004/0.005

#plotStationary & plotStates are unable to consider cosine + sine together
plotStationary(HMM_wels_fit_copy_int, plotCI=TRUE)
plot(HMM_wels_fit_copy_int, compact = T)
plotStates(HMM_wels_fit_copy_int)

plotStationary(HMM_wels_fit_copy, plotCI=TRUE)

################################################################################
################## Model Selection by AIC & PseudoResiduals #################### 
################################################################################

plotPR(HMM_wels_fit_copy_nocovars) #no covars 
plotPR(HMM_wels_fit_copy_temp) #~mean_temp
plotPR(HMM_wels_fit_copy_int)#~(cos_time + sin_time)*mean_temp #or HMM_wels_fit_int_stand
plotPR(HMM_wels_fit_copy) #~cos_time + sin_time
plotPR(HMM_wels_fit_copy_3covars) #~cos_time + sin_time + mean_temp

AIC(HMM_wels_fit_copy_int,HMM_wels_fit_copy_3covars,HMM_wels_fit_copy_temp,HMM_wels_fit_copy_nocovars,HMM_wels_fit_copy)


#extract states from model
HMM_final_copy_output<- HMM_full2_copy_prep2
HMM_final_copy_output$state_3s<-factor(moveHMM::viterbi(HMM_wels_fit_copy_int))

covars2<- HMM_final_copy_output %>% 
  dplyr::select(sin_time, cos_time,stand_mean_temp) %>% 
  as.data.frame()

HMM_final_output_copy_prob <- stationary(HMM_wels_fit_copy_int, covs = covars2)
colnames(HMM_final_output_copy_prob)<- c("prob_s1", "prob_s2", "prob_s3")

HMM_final_copy_output2<- cbind(HMM_final_copy_output, HMM_final_output_copy_prob)

#HMM_wels_fit_int
#State 1 -> inactive (low step length high turning angle)
#State 2 -> High Activity (high step length, low turning angle)
#State 3 -> Low Activity (medium step length, low turning angle)

#export final data
write_csv(HMM_final_copy_output2, "./data/final_model_output_data/HMM_final_output_TimeTempInt_UPDATED.csv")



################################################################################
#################  Visualizations of outputs for final model  ##################
################################################################################
#########################
##### Step Lengths ######
#########################
#c("#A50021", "#F5D683","#B4B8B8")
# Colours and labels
mycols <- c("cadetblue3", "darkolivegreen", "brown3")
state_labels <- c("2" = "High Activity", "3" = "Low Activity", "1" = "Inactivity")
state_levels <- c("High Activity", "Low Activity", "Inactivity")

# Number of states
nstate <- 3
states <- viterbi(HMM_wels_fit_copy_int)

# Reorder params to match labels (2: High, 3: Low, 1: Inactive)
reorder_idx <- c(2, 3, 1)

stepMean <- HMM_wels_fit_copy_int$mle$stepPar["mean",][reorder_idx]
stepSD   <- HMM_wels_fit_copy_int$mle$stepPar["sd",][reorder_idx]

stepShape <- stepMean^2 / stepSD^2
stepRate  <- stepMean / stepSD^2

# Step grid
stepgrid <- seq(min(HMM_final_copy_output2$step, na.rm=TRUE),
                max(HMM_final_copy_output2$step, na.rm=TRUE), length = 1000)

# Gamma density data
gamma_df <- data.frame()
for (s in 1:nstate) {
  gamma_df <- bind_rows(gamma_df, data.frame(
    stepgrid = stepgrid,
    density = dgamma(stepgrid, shape = stepShape[s], rate = stepRate[s]),
    state = factor(state_levels[s], levels = state_levels)
  ))
}

# Plot data
plot_data <- HMM_final_copy_output2 %>%
  filter(step != 0) %>%
  mutate(state = factor(state_labels[as.character(states)],
                        levels = state_levels))

# Plot
step_plot <- ggplot(plot_data, aes(x = step, fill = state, color = state)) +
  #geom_histogram(aes(y = after_stat(density)), bins = 20, color = "black", alpha = 0.5) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, boundary = 0, color = "black", alpha = 0.5)+
  geom_line(data = gamma_df, aes(x = stepgrid, y = density), linewidth = 2) +
  scale_fill_manual(values = mycols, name = "State") +
  scale_color_manual(values = mycols, name = "State") +
  theme_bw(base_size = 20) +
  labs(#title = "Step Length Distribution by State",
       x = "Step Length (m)",
       y = "Density")+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),,
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20),
        legend.position = 'none')



#########################
##### Turning Angles ####
#########################

# Reorder turning angle parameters
reorder_idx <- c(2, 3, 1)
angleMean <- HMM_wels_fit_copy_int$mle$anglePar["mean",][reorder_idx]
angleCon  <- HMM_wels_fit_copy_int$mle$anglePar["concentration",][reorder_idx]

# Grid for von Mises density
anglegrid <- seq(-pi, pi, length.out = 1000)

# Von Mises density data
von_mises_df <- data.frame()
for (s in 1:nstate) {
  von_mises_df <- bind_rows(von_mises_df, data.frame(
    anglegrid = anglegrid,
    density = dvm(anglegrid, mu = angleMean[s], kappa = angleCon[s]),
    state = factor(state_levels[s], levels = state_levels)
  ))
}

# Prepare plot data
plot_data <- HMM_final_copy_output2 %>%
  mutate(state = factor(state_labels[as.character(states)],
                        levels = state_levels))

range(plot_data$angle)

# Plot
angle_plot <- ggplot(plot_data, aes(x = angle, fill = state, color = state)) +
  #geom_histogram(aes(y = after_stat(density)), bins = 20, color = "black", alpha = 0.5) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, color = "black", alpha = 0.5)+
    geom_line(data = von_mises_df, aes(x = anglegrid, y = density, color = state), linewidth = 2) +
  scale_fill_manual(values = mycols, name = "State") +
  scale_color_manual(values = mycols, name = "State") +
  theme_bw(base_size = 20) +
  labs(#title = "Turning Angle Distribution by State",
       x = "Turning Angle (radians)",
       y = "Density")+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10)),,
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        axis.title.y = element_blank(),
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20))


##arrange hmm outputs together

ggarrange(step_plot, angle_plot, ncol = 2,
          labels = c('a)','b)'), 
          align = 'hv',
          label.x = 0.17,       # Horizontal position (0 = left, 1 = right)
          label.y = 0.98)       # Vertical position (0 = bottom, 1 = top))

#########################
####   Proportions   ####
#########################
HMM_final_output2

state_proportion <- HMM_final_copy_output2 %>%
  filter(ID != "T449282_1") %>% 
  group_by(ID, state_3s) %>%
  summarise(proportion = n(),.groups = "drop")

state_proportion2 <- state_proportion %>% 
  group_by(ID) %>% 
  mutate(sum_prop = sum(proportion),
            overall_prop = proportion / sum_prop) %>%
  mutate(state = factor(state_3s, levels = c(2, 3, 1), labels = state_labels))  # Ensure state is a factor for ggplot

state_proportion3 <- state_proportion2 %>% 
  mutate(inac_sum = ifelse(state == "Inactivity", overall_prop, 0),
         short_ID = stringr::str_extract(ID, "(?<=T449)\\d{3}"))

print(state_proportion3, n=40)

state_proportion3 %>% 
  summarise(totProp = sum(proportion, na.rm = TRUE))

sum(state_proportion3$proportion)

unique(state_proportion3$short_ID)

#rename labels 
label_map <- c("275" = "1", "280" = "2", "288" = "3", "268" = "4", "270" = "5", 
               "319" = "6", "278" = "7", "269" = "8", "274" = "9", "314" = "10", 
               "286" = "11", "276" = "12", "318" = "13")


#1 overall calculation of detection
ggplot(state_proportion3, aes(x = reorder(short_ID,-inac_sum), y = overall_prop, fill = state)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Catfish #", y = "Proportion of Detections", fill = "State") +
  scale_y_continuous(labels = scales::percent) +  # Converts y-axis to percentage
  scale_x_discrete(labels = label_map) +  # Optional: apply custom labels
    theme_bw(base_size = 20) +
  scale_fill_manual(values = c("cadetblue3", "darkolivegreen4", "brown3"))+
 # geom_text(aes(x = short_ID, y = 1.05, label = sum_prop),  # Position text above bars
  #          size = 4)+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),,
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
        #axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

#  geom_text(aes(label = proportion  ),   #to add # of detections per state
#          position = position_fill(vjust = 1.05),  # Positions text above bars
#          size = 4)   # Adjust text size


################################################################################
#################  Visualize state probs. of final model  ######################
#################  under different temperature conditions ######################
################################################################################

#plot state probabilites as a function of the combination of cosine and sine
covars<- HMM_full2_copy_prep %>% 
  dplyr::select(cos_time, sin_time) %>% 
  as.data.frame() %>%
  mutate(stand_mean_temp = mean(HMM_full2_copy_prep$stand_mean_temp)) 

covars_lowT<- HMM_full2_copy_prep %>% 
  dplyr::select(cos_time, sin_time) %>% 
  as.data.frame()%>%
  mutate(stand_mean_temp = min(HMM_full2_copy_prep$stand_mean_temp)) 

covars_highT<- HMM_full2_copy_prep %>% 
  dplyr::select(cos_time, sin_time) %>% 
  as.data.frame()%>%
  mutate(stand_mean_temp = max(HMM_full2_copy_prep$stand_mean_temp ))

day_hr<- HMM_full2_copy_prep %>% 
  dplyr::select(day_time) %>% 
  as.data.frame()


# Compute stationary probabilities for each covariate pair
stationary_probs2 <- stationary(HMM_wels_fit_copy_int, covs = covars)
stationary_probs_lowT2 <- stationary(HMM_wels_fit_copy_int, covs = covars_lowT)
stationary_probs_highT2 <- stationary(HMM_wels_fit_copy_int, covs = covars_highT)

# Combine the results with the time sequence
stationary_df_lowT <- data.frame(
  time_seq = day_hr,
  state1_prob = stationary_probs_lowT2[, 1],  # Probabilities for state 1
  state2_prob = stationary_probs_lowT2[, 2],  # Probabilities for state 2
  state3_prob = stationary_probs_lowT2[, 3]   # Probabilities for state 3
)

#rename states, pivot longer, and round datapoints to 1 for smoothing
stationary_df_lowT2<- stationary_df_lowT %>% 
  rename(inactivity = "state1_prob", 
         low_activity = "state3_prob",
         high_activity = "state2_prob") %>% 
  pivot_longer(cols = c(-day_time), #pivot all columns except under_state & ID
               names_to = "state", values_to = "prob")  %>% #put the col names into a col called "over_state", and the values into "overlap" col
  mutate(day_time_r = round(day_time),
         prob_r = round(prob, 100)) %>% 
  mutate(state =  factor(state, levels = c("high_activity", "low_activity", "inactivity")))


#join probabilities and time together
stationary_df_highT <- data.frame(
  time_seq = day_hr,
  state1_prob = stationary_probs_highT2[, 1],  # Probabilities for state 1
  state2_prob = stationary_probs_highT2[, 2],  # Probabilities for state 2
  state3_prob = stationary_probs_highT2[, 3]   # Probabilities for state 3
)

#rename states, pivot longer, and round datapoints to 1 for smoothing
stationary_df_highT2<- stationary_df_highT %>% 
  rename(inactivity = "state1_prob", 
         low_activity = "state3_prob",
         high_activity = "state2_prob") %>% 
  pivot_longer(cols = c(-day_time), #pivot all columns except under_state & ID
               names_to = "state", values_to = "prob") %>% #put the col names into a col called "over_state", and the values into "overlap" col
  mutate(day_time_r = round(day_time,1),
         prob_r = round(prob, 100)) %>% 
  mutate(state =  factor(state, levels = c("high_activity", "low_activity", "inactivity")))

#mean epilimnion temperature

# Combine the results with the time sequence
stationary_df_avgT <- data.frame(
  time_seq = day_hr,
  state1_prob = stationary_probs2[, 1],  # Probabilities for state 1
  state2_prob = stationary_probs2[, 2],  # Probabilities for state 2
  state3_prob = stationary_probs2[, 3]   # Probabilities for state 3
)

#rename states, pivot longer, and round datapoints to 1 for smoothing
stationary_df_avgT2<- stationary_df_avgT %>% 
  rename(inactivity = "state1_prob", 
         low_activity = "state3_prob",
         high_activity = "state2_prob") %>% 
  pivot_longer(cols = c(-day_time), #pivot all columns except under_state & ID
               names_to = "state", values_to = "prob")  %>% #put the col names into a col called "over_state", and the values into "overlap" col
  mutate(day_time_r = round(day_time),
         prob_r = round(prob, 100)) %>% 
  mutate(state =  factor(state, levels = c("high_activity", "low_activity", "inactivity")))


#calculate non standardized temperature values for interpretation purposes

lowT_original <- min(HMM_full2_copy_prep2$stand_mean_temp) * sigma + mu
#19.31 -- 19.71
highT_original <- max(HMM_full2_copy_prep2$stand_mean_temp) * sigma + mu
#23.78 -- 22.66
avgT_original <- mean(HMM_full2_copy_prep2$stand_mean_temp) * sigma + mu
#20.98 -- 20.99

range(HMM_full2_copy_prep$temp_mean)
sd(HMM_full2_copy2$temp_mean)

# Use same dataset used in HMM model
mu <- mean(HMM_full2_copy2$temp_mean, na.rm = TRUE)
sigma <- sd(HMM_full2_copy2$temp_mean, na.rm = TRUE)

# Match with the standardized covariates used in stationary()
lowT_original <- min(HMM_full2_copy2$stand_mean_temp, na.rm = TRUE) * sigma + mu
highT_original <- max(HMM_full2_copy2$stand_mean_temp, na.rm = TRUE) * sigma + mu
avgT_original <- mean(HMM_full2_copy2$stand_mean_temp, na.rm = TRUE) * sigma + mu


library(ggpattern)

#store here for now
gg_lowT
# Plot the stationary probabilities for each state
gg_lowT2<- ggplot(stationary_df_lowT2, aes(x = day_time_r, y = prob, color = state)) +
  #geom_rect(aes(xmin = c(2.19), xmax = c(4.17), ymin = 0, ymax = 1),
   #         fill = "darkgrey", col = "darkgrey")+
  #geom_rect(aes(xmin = c(3.02), xmax = c(4.17), ymin = 0, ymax = 1),
   #         fill = "lightgrey", col = "lightgrey")+
  #geom_rect(aes(xmin = 17.09, xmax = 19.09, ymin = 0, ymax = 1),
  #          fill = "lightgrey", col = "lightgrey")+
  #geom_rect(aes(xmin =17.51, xmax = 19.09, ymin = 0, ymax = 1),
  #          fill = "darkgrey", col = "darkgrey")+
  #geom_vline(xintercept = 4.17, linetype = "dashed", color = "darkgrey", linewidth = 1) +  # Dashed vertical line
  #geom_vline(xintercept = 19.09, linetype = "dashed", color = "lightgrey", linewidth = 1) +  # Dashed vertical line
  #geom_smooth(se = T)+
  labs( y = "Probability of Being in State", color = "State") +
  #ggtitle("Under Low Epilimnion Temperature Conditions (19.7 °C)") +
  theme_bw(base_size = 15)+
  scale_x_continuous(breaks = seq(0,24,4),
                     expand = c(0,0),
                     minor_breaks = seq(0,24,1))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  theme(plot.margin = margin(1, 1.5 , 1, 1.5, "cm"))+# Adds space around the plot
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        #axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        #legend.position = 'none'
        )  

gg_highT
gg_highT2<- ggplot(stationary_df_highT2, aes(x = day_time_r, y = prob, color = state)) +
  #geom_rect(aes(xmin = c(2.19), xmax = c(4.17), ymin = 0, ymax = 1),
  #          fill = "darkgrey", col = "darkgrey")+
  #geom_rect(aes(xmin = c(3.02), xmax = c(4.17), ymin = 0, ymax = 1),
  #          fill = "lightgrey", col = "lightgrey")+
  #geom_rect(aes(xmin = 17.09, xmax = 19.09, ymin = 0, ymax = 1),
  #          fill = "lightgrey", col = "lightgrey")+
  #geom_rect(aes(xmin =17.51, xmax = 19.09, ymin = 0, ymax = 1),
  #          fill = "darkgrey", col = "darkgrey")+
  #geom_vline(xintercept = 4.17, linetype = "dashed", color = "darkgrey", linewidth = 1) +  # Dashed vertical line
  #geom_vline(xintercept = 19.09, linetype = "dashed", color = "lightgrey", linewidth = 1) +  # Dashed vertical line
  #  geom_smooth(se = T)+
  labs(x = "Time (UTC Hours)", y = "Probability of Being in State", color = "State") +
  #ggtitle("Under High Epilimnion Temperature Conditions (22.7 °C)") +
  theme_bw(base_size = 15)+
  scale_x_continuous(breaks = seq(0,24,4),
                     expand = c(0,0),
                     minor_breaks = seq(0,24,1))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  theme(plot.margin = margin(1, 1.5 , 1, 1.5, "cm"))+ # Adds space around the plot
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        #axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = 'none')  


gg_avgT
gg_avgT2<- ggplot(stationary_df_avgT2, aes(x = day_time_r, y = prob, color = state)) +
  geom_rect(aes(xmin = c(2.19), xmax = c(4.17), ymin = 0, ymax = 1),
            fill = "darkgrey", col = "darkgrey")+
  geom_rect(aes(xmin = c(3.02), xmax = c(4.17), ymin = 0, ymax = 1),
            fill = "lightgrey", col = "lightgrey")+
  geom_rect(aes(xmin = 17.09, xmax = 19.09, ymin = 0, ymax = 1),
            fill = "lightgrey", col = "lightgrey")+
  geom_rect(aes(xmin =17.51, xmax = 19.09, ymin = 0, ymax = 1),
            fill = "darkgrey", col = "darkgrey")+
  geom_vline(xintercept = 4.17, linetype = "dashed", color = "darkgrey", linewidth = 1) +  # Dashed vertical line
  geom_vline(xintercept = 19.09, linetype = "dashed", color = "lightgrey", linewidth = 1) +  # Dashed vertical line
  geom_smooth(se = T)+
  labs(x = "Hour (UTC)", y = "Probability of Being in State", color = "State") +
  #ggtitle("Under Mean Epilimnion Temperature Conditions (21 °C)") +
  theme_bw(base_size = 20)+
  scale_x_continuous(breaks = seq(0,24,4),
                     expand = c(0,0),
                     minor_breaks = seq(0,24,1))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  theme(plot.margin = margin(1, 1.5 , 1, 1.5, "cm"))+ # Adds space around the plot
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        #axis.title.x = element_text(margin = margin(t = 10)),
    #axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    legend.text = element_text(size = 20),     # Increase legend text size
    legend.title = element_text(size = 20))
        #legend.position = 'none')  


6655
combined_plots_temp
combined_plots_temp2 <- ggarrange(gg_lowT2,gg_avgT2,gg_highT2, nrow = 1, ncol = 3, common.legend = T, align = "hv")

#combined_plots_temp2 <- annotate_figure(combined_plots_temp, top = text_grob("Diurnal Prediction of Stationary State Probabilities",face = "bold", size = 14))
combined_plots_temp3 <- annotate_figure(combined_plots_temp2, bottom = text_grob("Hour (UTC)", size = 15))

ggsave("gg_temps_plot.pdf", plot = combined_plots_temp2, width = 10, height = 8)

help(package=momentuHMM)
################################################################################
######## Plot Change in Epilimnion Temperature Through Study Period ############
################################################################################

epi_mean_temp_july <- epi_mean_temp 
epi_mean_temp_july$date <- date(epi_mean_temp_july$time) 
  
range(epi_mean_temp_july$temp_mean)
sd()

epi_mean_temp_july2 <- epi_mean_temp_july %>% 
  group_by(date) %>% 
  summarise(daily_mean_epiTemp = mean(temp_mean),
            SE = sd(temp_mean) / sqrt(n()))

epi_mean_temp_july2$

epi_mean_temp_july2 %>% 
  ggplot(aes(x = date, y = daily_mean_epiTemp))+
  geom_point()+
  geom_path()+
  # Add error bars (SE) to each point
  geom_errorbar(aes(ymin = daily_mean_epiTemp - SE, ymax = daily_mean_epiTemp + SE), width = 0.2) +
    theme_bw()+
  labs(x = "Date", y = "Mean Temperature (°C)",
       title = "Change in Epilimnion Temperature over Study Period")+
theme(axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)))
################################################################################
#################       Save Model Output into CSV        ######################
################################################################################

write_csv(stationary_df_lowT, "./data/final_model_output_data/HMM_final_predict_lowT.csv")

write_csv(stationary_df_highT, "./data/final_model_output_data/HMM_final_predict_highT.csv")





##########################################################################
####### No longer needed by good to keep for future code needs ###########
##########################################################################
#############         download other Covariates          ################# 
#############                                            ################# 
##########################################################################

detection_depth<- read_csv("./data/detection_depth.csv")
#combine detection data with cleaned HMM dataset
HMM_wels_clean_depth<-left_join(HMM_prep_clean, detection_depth, by = c("ID", "time"))

write_csv(HMM_wels_clean_depth,"./data/wels_filtered_3s_depth.csv")


HMM_wels_clean_depth$ID<- as.factor(HMM_wels_clean_depth$ID)
str(HMM_wels_clean_depth)

bot_depth <- read_csv("./data/pos_mean_wels_filter.csv")
bot_depth<- bot_depth %>% rename(ID = fishid, time = timestamp_5min)

sum(is.na(bot_depth$bottom_depth))

wels_depth_filtered_covar<- left_join(HMM_wels_clean_depth,bot_depth, by = c("ID", "time"))
str(wels_depth_filtered_covar)

write_csv(wels_depth_filtered_covar, "./data/wels_filtered_3s_covars.csv")



