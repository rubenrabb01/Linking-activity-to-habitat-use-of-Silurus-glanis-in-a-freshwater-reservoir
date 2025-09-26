#final clmm code
source("./scripts/libraries.R")
library(ordinal)
library(corrr)
library(car)

#perform Cumulative Link Mixed Models (CLMMs)

#load in hmm state data, and other descriptive parameter data
wels_data<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt_UPDATED.csv")

unique(wels_data$ID)

detection_depth<- read_csv("./data/detection_depth.csv")

bot_depth <- read_csv("./data/pos_mean_wels_filter.csv")
bot_depth <- bot_depth %>%
  rename(ID = fishid, time = timestamp_5min) %>% 
  dplyr::select(-depth)

dist2shore <- read_csv("./data/dist2shore_calc.csv")
dist2shore <- dist2shore %>% 
  dplyr::select(ID, time, dist_calc)

#combine all extracted / calculated habitat parms with cleaned HMM dataset (activity states)
wels_clmm_dat<-left_join(wels_data, detection_depth, by = c("ID", "time")) %>% 
  left_join(.,bot_depth, by = c("ID", "time")) %>% 
  left_join(.,dist2shore, by = c("ID", "time"))

#select only necessary columns
wels_clmm_dat2 <- wels_clmm_dat %>% 
  dplyr::select(ID, time, dp, state_3s, depth, bottom_depth, dist_calc, day_time)


wels_clmm_dat2 <- wels_clmm_dat2 %>% 
  rename( "diel_period" = dp)

#############################################
## Estimate missing values in bottom depth ##
#############################################

clmm_full_estim <- wels_clmm_dat2 %>% 
  filter(ID != "T449282_1")

# Start by making a copy and converting to data.table
clmm_full_estim2 <- data.table(clmm_full_estim)

# Create timestamp column from 'time'
clmm_full_estim2[, timestamp := time]

#Estimate missing depth values by fish ID and with no timing gaps later than 45 minutes
# Make sure data is sorted by timestamp
setorder(clmm_full_estim2, timestamp)

# Define maximum allowed gap in seconds (45 minutes)
max_gap_seconds <- 45 * 60

# New grouped interpolation function
interpolate_column<- function(dt, groupcol, colname, timecol = "timestamp", maxgap = max_gap_seconds) {
  
  dt <- copy(dt) # To avoid modifying original by reference
  
  # Process separately for each group
  dt[, (colname) := {
    
    x <- get(timecol)
    y <- get(colname)
    
    # Find indices of non-NA values
    non_na_idx <- which(!is.na(y))
    
    if (length(non_na_idx) == 0 || length(non_na_idx) == length(y)) {
      y
    } else {
      for (i in which(is.na(y))) {
        
        prev_idx <- suppressWarnings(max(non_na_idx[non_na_idx < i], na.rm = TRUE))
        next_idx <- suppressWarnings(min(non_na_idx[non_na_idx > i], na.rm = TRUE))
        
        if (!is.infinite(prev_idx) && !is.infinite(next_idx)) {
          gap_seconds <- as.numeric(difftime(x[next_idx], x[prev_idx], units = "secs"))
          
          if (!is.na(gap_seconds) && gap_seconds <= maxgap) {
            frac <- as.numeric(difftime(x[i], x[prev_idx], units = "secs")) / gap_seconds
            interpolated_value <- y[prev_idx] + frac * (y[next_idx] - y[prev_idx])
            y[i] <- interpolated_value
          }
        }
      }
      y
    }
    
  }, by = groupcol]
  
  return(dt)
}

# Apply interpolation to both 'depth' and 'bottom_depth'
clmm_full_estim3 <- interpolate_column(clmm_full_estim2, "ID", "depth")
clmm_full_estim4 <- interpolate_column(clmm_full_estim3, "ID", "bottom_depth")

#calculate distance to bottom 
clmm_full_estim5<-clmm_full_estim4 %>% 
  mutate(dist2bot = (bottom_depth - depth))

max(na.omit(clmm_full_estim2$dist2bot))


#ensure columns are labelled correctly
clmm_full_estim5$ID<- as.factor(clmm_full_estim5$ID)
clmm_full_estim5$diel_period <- as.factor(clmm_full_estim5$diel_period)

#create new column which correctly identifies which state corresponds to what activity lvl
clmm_full_estim6 <- clmm_full_estim5 %>% 
  mutate(activity_state = case_when(
    state_3s == "3" ~ "Low_active",
    state_3s == "2" ~ "High_active",
    state_3s == "1" ~ "Inactive"))

#assign response variable as an ordered factor
clmm_full_estim6$activity_state <- factor(clmm_full_estim6$activity_state,
                                        levels = c("Inactive", "Low_active", "High_active"),
                                        ordered = TRUE)

#calculate cyclic time

#convert timestamp to day_hr, and to seconds
clmm_full_estim6 <- clmm_full_estim6 %>% 
  mutate(day_time = hour(time)+(minute(time)/60)) 

clmm_full_estim6<- clmm_full_estim6 %>% 
  mutate(hour = hour(time),
         minute = minute(time), 
         sec = second(time)) 

clmm_full_estim6<- clmm_full_estim6 %>% 
  mutate(day_time_sec = hour*60*60 + minute*60 + sec)

#calculate cosine time variables to make time cyclical
clmm_full_estim6<- clmm_full_estim6 %>%
  mutate(cos_time = cos(2*pi*day_time_sec/86400),
                    sin_time = sin(2*pi*day_time_sec/86400))



################################################################################
############## Run individual clmms based on different hypotheses ############## 
################################################################################

##################################
#### CLMM with only diel period ####
##################################
clmm_Diel<- clmm(activity_state ~ diel_period + (1|ID), data = clmm_full_estim6, Hess = TRUE)
summary(clmm_Diel)
#tab_model(clmm_Diel,clmm_H1)



##################################
#### CLMM with cos + sin time ####
##################################
#clmm_time<- clmm(activity_state ~ sin_time + cos_time + (1|ID), data = clmm_full_estim6, Hess = TRUE)
#summary(clmm_time)


#####################
######## H1 #########
#####################
#H1: catfish are more likely to be closer to shore when they are inactive during the day, than when they are active
#focuses on dist_calc

#determine range to see if there are any outlier values
range(clmm_full_estim6$dist_calc)

#remove observations greater than 210m due to error in distance calculations
clmm_full_estim7 <- clmm_full_estim6 %>% 
  filter(dist_calc < 210) 

#scale distance 2 shore due to large range in values (3.386e-05 --> 2.045e+02)
clmm_full_estim7$dist_calc_scaled = as.numeric(scale(clmm_full_estim7$dist_calc))

summary(clmm_full_estim7$dist_calc)
mean_dist2shore <- mean(clmm_full_estim7$dist_calc)
sd_dist2shore <- sd(clmm_full_estim7$dist_calc)
#means 1 unit increase in dist_calc_scaled corresponds to moving 40 m away from shore
#a value of 0 in dist_calc_scaled corresponds to 37.5 m (mean distance)

#run full model
clmm_H1 <- clmm(activity_state ~ dist_calc_scaled*diel_period + (1|ID), data = clmm_full_estim7, Hess = TRUE)
#summary(clmm_H1)

#model selection - drop non significant covars. 
#CLMM_clmm_H1_drop <- drop1(clmm_H1)
#CLMM_clmm_H1_drop[CLMM_clmm_H1_drop$AIC == min(CLMM_clmm_H1_drop$AIC), ]
#no variables to drop

#run model without diel period
clmm_H1_noDiel <- clmm(activity_state ~ dist_calc_scaled + (1|ID), data = clmm_full_estim7, Hess = TRUE)
#summary(clmm_H1_noDiel)

#run model with cyclic time
#clmm_H1_cycTime <- clmm(activity_state ~ dist_calc_scaled*diel_period + cos_time + sin_time + (1|ID), data = clmm_full_estim4, Hess = TRUE)

#run model with additive diel period
clmm_H1_addDiel <- clmm(activity_state ~ dist_calc_scaled + diel_period + (1|ID), data = clmm_full_estim7, Hess = TRUE)

tab_model(clmm_H1_noDiel,clmm_H1_addDiel,clmm_H1,show.aic = T)

############################
## Check CLMM Assumptions ##
############################
#check collinearity
library(performance)
check_collinearity(clmm_H1) #low collinearity, assumption met

#check proportional odds (parallel slope)
# Fit a non-proportional odds model (nominal = TRUE)
clmm_H1_nonprop <- clmm(activity_state ~ dist_calc_scaled * diel_period + (1 | ID), 
                        data = clmm_full_estim7, Hess = TRUE, nominal = ~ dist_calc_scaled + diel_period)

# Compare models
anova(clmm_H1, clmm_H1_nonprop)
#predictor effects consistent across response threshold: assumption met

#random effects with normal distribution 
# Extract random effects
random_effects <- ranef(clmm_H1)$ID
hist(random_effects$`(Intercept)`, main="Random Effects Distribution", xlab="Random Effect (ID)", col="lightblue")

shapiro.test(random_effects$`(Intercept)`) #NON sig = assumption met

#check model convergence 
summary(clmm_H1)  # Look for large standard errors
range(clmm_full_estim7$dist_calc)




###########################
## visualize predictions ##
###########################
#run clmm2 to visualize predictions predictions
clmm_full_estim7$ID <- as.factor(clmm_full_estim7$ID)

# Step 1: Define control parameters with a different optimizer
control_params <- clmm2.control(
  method = "nlminb",  # Change to "bfgs", "nlminb", or other optimizers
  maxit = 1000,  # Increase max iterations
  trace = 0  # Suppress trace output if it's too verbose
)

CLMM2_H1_end <- clmm2(activity_state ~ dist_calc_scaled * diel_period,
                      random = ID, 
                      Hess = TRUE, 
                      control = control_params,
                      data = clmm_full_estim7)

summary(CLMM2_H1_end)

#add cyclic time to clmm
#clmm_full_estim7 <- clmm_full_estim7 %>% 
#  mutate(circ_time = atan2(cos_time,sin_time))


#CLMM2_H1_end_time <- clmm2(activity_state ~ dist_calc_scaled * dp + cos_time + sin_time,
 #                     random = ID, 
  #                    Hess = TRUE, 
   #                   control = control_params,
    #                  data = clmm_full_estim7)

#summary(CLMM2_H1_end_time)

#acf(CLMM2_H1_end_time$location)

#try time as 1 variable
#CLMM2_H1_end_time2 <- clmm2(activity_state ~ dist_calc_scaled * dp + circ_time,
#                           random = ID, 
 #                          Hess = TRUE, 
  #                         control = control_params,
   #                        data = clmm_full_estim7)

#summary(CLMM2_H1_end_time2)
#acf(CLMM2_H1_end_time2$location)

#time as 1 variable is not preferable (based on AIC & log likelihood) to time as 2 vars (cos and sin)

# Compute Pearson residuals
#pearson_resid <- residuals(CLMM2_H1_end_time, type = "pearson")

# Check residuals
#acf(pearson_resid, na.action = na.omit)

#ranef_resid <- ranef(CLMM2_H1_end_time$location)$ID  
#acf(ranef_resid, na.action = na.omit)



###########################
## Visualize Predictions ##
###########################

# Create New Data for Prediction
newdata_dist2shore <- expand.grid(
  dist_calc_scaled = seq(from = min(clmm_full_estim7$dist_calc_scaled), 
                         to = max(clmm_full_estim7$dist_calc_scaled), 
                         length.out = 100),
  diel_period = c("day", "night")  # Adding diel period as a factor with two levels: day and night
)

newdata_dist2shore$diel_period <- factor(newdata_dist2shore$diel_period, 
                                         levels = c("day", "night"))

# Predict Probabilities for Each Activity State (Inactive, Low_active, High_active)
predict_dist2shore <- sapply(levels(clmm_full_estim7$activity_state), function(x) {
  newdata1 <- merge(newdata_dist2shore, 
                    data.frame(activity_state = factor(x, levels = levels(clmm_full_estim7$activity_state))),
                    by = NULL)  # Ensures correct mapping of diel_period
  
  predict(CLMM2_H1_end, newdata = newdata1, type = "response")  # Get probabilities using your model
})

# Combine new data with predictions
predict_dist2shore_2 <-  as.data.frame(cbind(newdata_dist2shore, predict_dist2shore))

str(predict_dist2shore)

predict_dist2shore_2$dist_calc_unscaled <- predict_dist2shore_2$dist_calc_scaled * sd_dist2shore + mean_dist2shore

# Convert to long format for easier plotting
predict_dist2shore_long <- gather(predict_dist2shore_2, 'Inactive', 'Low_active', 'High_active', key = "State", value = "Probability")

# Reorder legend labels (High_active first, then Low_active, then Inactive)
predict_dist2shore_long$State <- factor(predict_dist2shore_long$State, 
                                        levels = c("High_active", "Low_active", "Inactive"))



# Visualize the results
H1_fig <- ggplot(predict_dist2shore_long, aes(x = dist_calc_unscaled, y = Probability, color = State, linetype = dp)) +
  geom_line(linewidth = 1.5) +
  labs(x = "Distance from Shore (m)", 
       y = "Probability of being in State",
       color = "State",
       linetype = "Diel Period") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  ylim(0,1)+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20))  



#####################
######## H2 #########
#####################
#H2: catfish are more likely to in deeper habitats when they are active, than when they are inactive. 
#they are more likely to be in deeper habitats during the night than during the day
#based on bottom depth

#determine range to see if there are any outlier values
range(clmm_full_estim7$bottom_depth)

#remove observations less than 0m due to error in bottom depth interpolation calculations
clmm_full_estim8 <- clmm_full_estim7 %>% 
  filter(bottom_depth > 0,
         !is.na(bottom_depth)) 

#scale distance 2 shore due to large range in values (3.386e-05 --> 2.045e+02)
clmm_full_estim8$botdepth_scaled = as.numeric(scale(clmm_full_estim8$bottom_depth))

summary(clmm_full_estim8$bottom_depth)
sd(clmm_full_estim8$bottom_depth)

mean_botDepth <- mean(clmm_full_estim8$bottom_depth)
sd_botDepth <- sd(clmm_full_estim8$bottom_depth)
#means 1 unit increase in botdepth_scaled corresponds to moving 8.8 m away from shore
#a value of 0 in botdepth_scaled corresponds to 9.8 m (mean distance)


#run full model
clmm_H2 <- clmm(activity_state ~ botdepth_scaled*diel_period + (1|ID), data = clmm_full_estim8, Hess = TRUE)
#summary(clmm_H2)

#run model with no diel period 
clmm_H2_noDiel <- clmm(activity_state ~ botdepth_scaled + (1|ID), data = clmm_full_estim8, Hess = TRUE)


#run model with cyclic time
#clmm_H2_cycTime <- clmm(activity_state ~ botdepth_scaled*dp + cos_time + sin_time + (1|ID), data = clmm_full_estim8, Hess = TRUE)

#run model with additive diel
clmm_H2_addDiel <- clmm(activity_state ~ botdepth_scaled + diel_period + (1|ID), data = clmm_full_estim8, Hess = TRUE)

tab_model(clmm_H2_noDiel,clmm_H2_addDiel,clmm_H2,show.aic = T)

############################
## Check CLMM Assumptions ##
############################
#check collinearity
library(performance)
check_collinearity(clmm_H2) #low collinearity, assumption met

#check proportional odds (parallel slope)
# Fit a non-proportional odds model (nominal = TRUE)
clmm_H2_nonprop <- clmm(activity_state ~ botdepth_scaled * dp + (1 | ID), 
                        data = clmm_full_estim8, Hess = TRUE, nominal = ~ botdepth_scaled + dp)

# Compare models
anova(clmm_H2, clmm_H2_nonprop)
#predictor effects consistent across response threshold: assumption met

#random effects with normal distribution 
# Extract random effects
random_effects <- ranef(clmm_H2)$ID
hist(random_effects$`(Intercept)`, main="Random Effects Distribution", xlab="Random Effect (ID)", col="lightblue")

shapiro.test(random_effects$`(Intercept)`) #NON sig = assumption met




###########################
######   run clmm2    #####
###########################
#run clmm2 to visualize predictions predictions
clmm_full_estim8$ID <- as.factor(clmm_full_estim8$ID)

# Step 1: Define control parameters with a different optimizer
control_params <- clmm2.control(
  method = "nlminb",  # Change to "bfgs", "nlminb", or other optimizers
  maxit = 1000,  # Increase max iterations
  trace = 0  # Suppress trace output if it's too verbose
)

CLMM2_H2_end <- clmm2(activity_state ~ botdepth_scaled * diel_period , 
                      random = ID, 
                      Hess = TRUE, 
                      control = control_params,
                      data = clmm_full_estim8)

summary(CLMM2_H2_end)

###########################
## Visualize Predictions ##
###########################

# Create New Data for Prediction
newdata_botdepth_scaled <- expand.grid(
  botdepth_scaled = seq(from = min(clmm_full_estim8$botdepth_scaled), 
                         to = max(clmm_full_estim8$botdepth_scaled), 
                         length.out = 100),
  diel_period = c("day", "night")  # Adding diel period as a factor with two levels: day and night
)

newdata_botdepth_scaled$dp <- factor(newdata_botdepth_scaled$diel_period, 
                                         levels = c("day", "night"))

# Predict Probabilities for Each Activity State (Inactive, Low_active, High_active)
predict_botdepth_scaled <- sapply(levels(clmm_full_estim8$activity_state), function(x) {
  newdata1 <- merge(newdata_botdepth_scaled, 
                    data.frame(activity_state = factor(x, levels = levels(clmm_full_estim8$activity_state))),
                    by = NULL)  # Ensures correct mapping of diel_period
  
  predict(CLMM2_H2_end, newdata = newdata1, type = "response")  # Get probabilities using your model
})

# Combine new data with predictions
predict_botdepth_scaled_2 <-  as.data.frame(cbind(newdata_botdepth_scaled, predict_botdepth_scaled))

predict_botdepth_scaled_2$botdepth_unscaled <- predict_botdepth_scaled_2$botdepth_scaled * sd_botDepth + mean_botDepth


# Convert to long format for easier plotting
predict_botdepth_scaled_long <- gather(predict_botdepth_scaled_2, 'Inactive', 'Low_active', 'High_active', key = "State", value = "Probability")

# Reorder legend labels (High_active first, then Low_active, then Inactive)
predict_botdepth_scaled_long$State <- factor(predict_botdepth_scaled_long$State, 
                                        levels = c("High_active", "Low_active", "Inactive"))

# Visualize the results
H2_fig <- ggplot(predict_botdepth_scaled_long, aes(x = botdepth_unscaled, y = Probability, color = State, linetype = dp)) +
  geom_line(linewidth = 1.5) +
  labs(x = "Reservoir Bottom Depth (m)", 
       #y = "Probability of being in State",
       color = "State",
       linetype = "Diel Period") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  ylim(0,1)+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20))





#####################
######## H3 #########
#####################
#H3: catfish are more likely to be closer to the bottom substrate when they are inactive, than when they are active
#their distance to the bottom substrate will change with diel period

#determine range to see if there are any outlier values
range(clmm_full_estim7$dist2bot)

#remove observations that are less than 0 due to mismatch in bottom depth and depth measures
clmm_full_estim_dist2bot <- clmm_full_estim7 %>% 
  filter(dist2bot > 0, 
         !is.na(dist2bot)) 

#scale distance 2 shore due to large hessian value after performing clmm with non scaled data
clmm_full_estim_dist2bot$dist2bot_scaled = as.numeric(scale(clmm_full_estim_dist2bot$dist2bot))
range(clmm_full_estim_dist2bot$dist2bot_scaled)


summary(clmm_full_estim_dist2bot$dist2bot)

mean_dist2bot <- mean(clmm_full_estim_dist2bot$dist2bot)
sd_dist2bot <- sd(clmm_full_estim_dist2bot$dist2bot)
#a value of 0 in dist2bot_scaled corresponds to 8.4 m (mean distance from bottom)

sd(clmm_full_estim_dist2bot$dist2bot)
#means 1 unit increase in dist_calc_scaled corresponds to vertically moving 8.9m 
#away from the bottom of the reservoir


#run full model
clmm_H3 <- clmm(activity_state ~ dist2bot_scaled*diel_period + (1|ID), data = clmm_full_estim_dist2bot, Hess = TRUE)
#summary(clmm_H3)


#run model without diel period
clmm_H3_noDiel <- clmm(activity_state ~ dist2bot_scaled + (1|ID), data = clmm_full_estim_dist2bot, Hess = TRUE)


#run model with cyclic time
#clmm_H3_cycTime <- clmm(activity_state ~ dist2bot_scaled*dp + cos_time + sin_time + (1|ID), data = clmm_full_estim_dist2bot, Hess = TRUE)

#run model with additive diel
clmm_H3_addDiel <- clmm(activity_state ~ dist2bot_scaled + diel_period + (1|ID), data = clmm_full_estim_dist2bot, Hess = TRUE)

#AIC(clmm_H3_addDiel)

tab_model(clmm_H3_noDiel,clmm_H3_addDiel,clmm_H3,show.aic = T)

tab_model(clmm_H3)
############################
## Check CLMM Assumptions ##
############################
#check collinearity
library(performance)
check_collinearity(clmm_H3) #low collinearity, assumption met

#check proportional odds (parallel slope)
# Fit a non-proportional odds model (nominal = TRUE)
clmm_H3_nonprop <- clmm(activity_state ~ dist2bot*dp + (1|ID), data = clmm_full_estim_dist2bot, Hess = TRUE,
                        nominal = ~ dist2bot + dp)

# Compare models
anova(clmm_H3, clmm_H3_nonprop)
#predictor effects consistent across response threshold: assumption met

#random effects with normal distribution 
# Extract random effects
random_effects <- ranef(clmm_H3)$ID
hist(random_effects$`(Intercept)`, main="Random Effects Distribution", xlab="Random Effect (ID)", col="lightblue")

shapiro.test(random_effects$`(Intercept)`) #NON sig = assumption met

#check model convergence 
summary(clmm_H3)  # Look for large standard errors

###########################
## visualize predictions ##
###########################
#run clmm2 to visualize predictions predictions
clmm_full_estim_dist2bot$ID <- as.factor(clmm_full_estim_dist2bot$ID)

# Step 1: Define control parameters with a different optimizer
control_params <- clmm2.control(
  method = "nlminb",  # Change to "bfgs", "nlminb", or other optimizers
  maxit = 1000,  #  max iterations
  trace = 0  # Suppress trace output if it's too verbose
)

CLMM2_H3_end <- clmm2(activity_state ~ dist2bot_scaled * diel_period , 
                      random = ID, 
                      Hess = TRUE, 
                      control = control_params,
                      data = clmm_full_estim_dist2bot)

summary(CLMM2_H3_end)


###########################
## Visualize Predictions ##
###########################

# Create New Data for Prediction
newdata_dist2bot_scaled <- expand.grid(
  dist2bot_scaled = seq(from = min(clmm_full_estim_dist2bot$dist2bot_scaled), 
                         to = max(clmm_full_estim_dist2bot$dist2bot_scaled), 
                         length.out = 100),
  diel_period = c("day", "night")  # Adding diel period as a factor with two levels: day and night
)

newdata_dist2bot_scaled$diel_period <- factor(newdata_dist2bot_scaled$diel_period, 
                                         levels = c("day", "night"))

# Predict Probabilities for Each Activity State (Inactive, Low_active, High_active)
predict_dist2bot_scaled <- sapply(levels(clmm_full_estim_dist2bot$activity_state), function(x) {
  newdata1 <- merge(newdata_dist2bot_scaled, 
                    data.frame(activity_state = factor(x, levels = levels(clmm_full_estim_dist2bot$activity_state))),
                    by = NULL)  # Ensures correct mapping of diel_period
  
  predict(CLMM2_H3_end, newdata = newdata1, type = "response")  # Get probabilities using your model
})

# Combine new data with predictions
predict_dist2bot_2 <-  as.data.frame(cbind(newdata_dist2bot_scaled, predict_dist2bot_scaled))

predict_dist2bot_2$dist2bot_unscaled <- predict_dist2bot_2$dist2bot_scaled * sd_dist2bot + mean_dist2bot


# Convert to long format for easier plotting
predict_dist2bot_long <- gather(predict_dist2bot_2, 'Inactive', 'Low_active', 'High_active', key = "State", value = "Probability")

# Reorder legend labels (High_active first, then Low_active, then Inactive)
predict_dist2bot_long$State <- factor(predict_dist2bot_long$State, 
                                        levels = c("High_active", "Low_active", "Inactive"))

# Visualize the results
H3_fig <- ggplot(predict_dist2bot_long, aes(x = dist2bot_unscaled, y = Probability, color = State, linetype = diel_period)) +
  geom_line(linewidth = 1.5) +
  labs(x = "Distance to Reservoir Bottom (m)",
       y = "Probability of being in State",
       color = "State",
       linetype = "Diel Period") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  ylim(0,1)+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20))  


#####################
######## H4 #########
#####################
#H4: A catfish's depth probability will change based on the combined influence of
# their activity state and diel period

#determine range to see if there are any outlier values
range(clmm_full_estim7$depth)

clmm_full_estim9 <- clmm_full_estim7 %>% 
  filter(!is.na(depth)) 


#run full model
clmm_H4 <- clmm(activity_state ~ depth*diel_period + (1|ID), data = clmm_full_estim9, Hess = TRUE)
summary(clmm_H4)

#run model with no Diel period
clmm_H4_noDiel <- clmm(activity_state ~ depth + (1|ID), data = clmm_full_estim9, Hess = TRUE)


#run model with cyclic time
#clmm_H4_cycTime <- clmm(activity_state ~ depth*diel_period + cos_time + sin_time + (1|ID), data = clmm_full_estim3, Hess = TRUE)

#run model with additive diel
clmm_H4_addDiel <- clmm(activity_state ~ depth + diel_period + (1|ID), data = clmm_full_estim9, Hess = TRUE)

tab_model(clmm_H4_noDiel,clmm_H4_addDiel,clmm_H4,show.aic = T)



############################
## Check CLMM Assumptions ##
############################
#check collinearity
check_collinearity(clmm_H4) #low collinearity, assumption met

#check proportional odds (parallel slope)
# Fit a non-proportional odds model (nominal = TRUE)
clmm_H4_nonprop <- clmm(activity_state ~ depth*dp + (1|ID), data = clmm_full_estim9, Hess = TRUE,
                        nominal = ~ depth + dp)

# Compare models
anova(clmm_H4, clmm_H4_nonprop)
#predictor effects consistent across response threshold: assumption met

#random effects with normal distribution 
# Extract random effects
random_effects <- ranef(clmm_H4)$ID
hist(random_effects$`(Intercept)`, main="Random Effects Distribution", xlab="Random Effect (ID)", col="lightblue")

shapiro.test(random_effects$`(Intercept)`) #NON sig = assumption met

#check model convergence 
summary(clmm_H4)  # Look for large standard errors


###########################
## visualize predictions ##
###########################
#run clmm2 to visualize predictions predictions
clmm_full_estim9$ID <- as.factor(clmm_full_estim9$ID)

# Step 1: Define control parameters with a different optimizer
control_params <- clmm2.control(
  method = "nlminb",  # Change to "bfgs", "nlminb", or other optimizers
  maxit = 1000,  #  max iterations
  trace = 0  # Suppress trace output if it's too verbose
)

CLMM2_H4_end <- clmm2(activity_state ~ depth * diel_period , 
                      random = ID, 
                      Hess = TRUE, 
                      control = control_params,
                      data = clmm_full_estim9)

summary(CLMM2_H4_end)


###########################
## Visualize Predictions ##
###########################

# Create New Data for Prediction
newdata_depth <- expand.grid(
  depth = seq(from = min(clmm_full_estim9$depth), 
                        to = max(clmm_full_estim9$depth), 
                        length.out = 100),
  diel_period = c("day", "night")  # Adding diel period as a factor with two levels: day and night
)

newdata_depth$dp <- factor(newdata_depth$diel_period, 
                                              levels = c("day", "night"))


# Predict Probabilities for Each Activity State (Inactive, Low_active, High_active)
predict_depth <- sapply(levels(clmm_full_estim9$activity_state), function(x) {
  newdata1 <- merge(newdata_depth, 
                    data.frame(activity_state = factor(x, levels = levels(clmm_full_estim9$activity_state))),
                    by = NULL)  # Ensures correct mapping of diel_period
  
  predict(CLMM2_H4_end, newdata = newdata1, type = "response")  # Get probabilities using your model
})

# Combine new data with predictions
predict_depth_2 <-  as.data.frame(cbind(newdata_depth, predict_depth))

# Convert to long format for easier plotting
predict_depth_long <- gather(predict_depth_2, 'Inactive', 'Low_active', 'High_active', key = "State", value = "Probability")

# Reorder legend labels (High_active first, then Low_active, then Inactive)
predict_depth_long$State <- factor(predict_depth_long$State, 
                                      levels = c("High_active", "Low_active", "Inactive"))

# Visualize the results
H4_fig <- ggplot(predict_depth_long, aes(x = depth, y = Probability)) +
  geom_line(aes(color = State, linetype = diel_period), linewidth = 1.5) +
  labs(x = "Fish Depth (m)", 
       #y = "Probability of being in State",
       color = "State",
       linetype = "Diel Period") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("cadetblue3", "darkolivegreen","brown3"), labels = c("high activity","low activity","inactivity"))+
  ylim(0,1)+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Increase legend text size
        legend.title = element_text(size = 20))     # Increase legend title size)  

##F5FACD
###FEF9B5
##DCE7D7
###F5D683
##FF8E32

#arrange all plots in  1 fig
combined_plot <- ggarrange(H1_fig, H2_fig, H3_fig, H4_fig, 
                           ncol = 2, nrow = 2,
                           align = "hv",
                           labels = c("a)", "b)", "c)", "d)"),  # Add panel labels
                           label.x = 0.17,       # Horizontal position (0 = left, 1 = right)
                           label.y = 0.98,       # Vertical position (0 = bottom, 1 = top)
                           common.legend = TRUE, legend = "top")

# Add common y-axis label
final_plot <- annotate_figure(combined_plot,
                              left = text_grob("Probability of Being in State", 
                                               rot = 90, vjust = 1, size = 20))

# Print the final plot
print(final_plot)

### CHECK AIC

AIC(clmm_time,clmm_Diel,clmm_H1,clmm_H1_noDiel,clmm_H1_cycTime,clmm_H2,clmm_H2_noDiel,clmm_H2_cycTime,clmm_H3, clmm_H3_noDiel, clmm_H3_cycTime, clmm_H4, clmm_H4_noDiel,clmm_H4_cycTime)

tab_model(clmm_H1,clmm_H2,clmm_H3, clmm_H4, show.aic = T)
tab_model(CLMM2_H1_end, CLMM2_H2_end, CLMM2_H3_end, CLMM2_H4_end) #clmm2 and clmm have the same model output

summary(clmm_H1)
summary(CLMM2_H1_end)

#test
# Arrange the plots with a shared legend
combined_plot <- ggarrange(H1_fig, H2_fig, H3_fig, H4_fig, 
                           ncol = 2, nrow = 2,align = "hv",
                           common.legend = TRUE, legend = "top")

# Add common y-axis label
final_plot <- annotate_figure(combined_plot,
                              left = text_grob("Probability of Being in State", 
                                               rot = 90, vjust = 1, size = 15))

# Print the final plot
print(final_plot)


###################
###### Extra ######
###################
###############################################
####    CLMM with Day_Time  & Epilim temp  ####
###############################################
#download mean epilimnion data 
epi_mean_temp<- read_csv("./data/mean_temp_6m.csv")
epi_mean_temp<- epi_mean_temp %>% 
  rename(time = hd_timestamp_utc) %>% 
  filter(time > "2017-07-01 UTC" & time < "2017-09-01 UTC")


#merge will cleaned dataset
clmm_full_estim3_epiTemp<- full_join(clmm_full_estim3,epi_mean_temp, by = c("time"))

#standardize covariates
clmm_full_estim3_epiTemp2<- clmm_full_estim3_epiTemp %>% 
  na.omit(temp_mean) %>% 
  mutate(stand_mean_temp = (temp_mean - mean(temp_mean))/sd(temp_mean))

clmm_day_hr_temp <- clmm(activity_state ~ stand_mean_temp*day_time + (1|ID), data = clmm_full_estim3_epiTemp2, Hess = TRUE)
summary(clmm_day_hr_temp)

# Create new data frame with time range (e.g., 0 to 24 hours)
time_seq <- seq(0, 24, by = 1)

# Generate prediction data for low & high temp
new_data_lowT <- data.frame(day_time = time_seq, stand_mean_temp = min(clmm_full_estim3_epiTemp2$stand_mean_temp))
new_data_highT <- data.frame(day_time = time_seq, stand_mean_temp = max(clmm_full_estim3_epiTemp2$stand_mean_temp))

# Predict probabilities for each activity state
preds_lowT <- predict(clmm_day_hr_temp, new_data_lowT, type = "prob") #predict function doesnt work for clmm models, resort back to direct HMM predictions
preds_highT <- predict(clmm_day_hr_temp, new_data_highT, type = "prob")



#Scale fish ID variance
resid_var <- 3.29 #fixed residual variance for a logit model

#H1 Intra-CLass Correlation variance (ICC)
fish_var <- 0.2977
ICC <- fish_var/(fish_var + resid_var)
#0.08297795
#only 8.3% of total variance for distance to shore*diel period is due to differences between fish 
#Thus, most of the variation is within individuals 

#H2 ICC
fish_var <- 0.274
ICC <- fish_var/(fish_var + resid_var)
#0.07687991
#only 7.7% of total variance for bottom depth*diel period is due to differences between fish 
#Thus, most of the variation is within individuals 

#H3 ICC
fish_var <- 0.3868
ICC <- fish_var/(fish_var + resid_var)
#only 10.5% of total variance for distance to bottom*diel period is due to differences between fish 
#Thus, most of the variation is within individuals 

#H4 ICC
fish_var <- 0.1571
ICC <- fish_var/(fish_var + resid_var)
#only 4.5% of total variance for fish depth*diel period is due to differences between fish 
#Thus, most of the variation is within individuals 


##Create plots of Estimated Marginal Means for each habitat variable
library(emmeans)
emm_H1 <- emmeans(clmm_H1, ~ dist_calc_scaled * diel_period, type = 'response') 
emm_H1_df <- as.data.frame(emm_H1)

ggplot(emm_H1_df, aes(x = dist_calc_scaled, y = emmean, color = diel_period, group = diel_period)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  geom_line(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = SE, ymax = SE), 
                width = 0.1, position = position_dodge(0.2)) +
  labs(x = "Variable 1",
       y = "Estimated Marginal Mean (Probability)",
       color = "Variable 2",
       title = "Estimated Marginal Means with Interaction") +
  theme_minimal()

#Create plots of Odds Ratio for each habitat variable

library(broom.mixed)

# Tidy and transform each model
oddsH1 <- tidy(clmm_H1, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 1")
oddsH2 <- tidy(clmm_H2, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 2")
oddsH3 <- tidy(clmm_H3, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 3")
oddsH4 <- tidy(clmm_H4, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 4")

# Combine into one data frame
odds_all <- bind_rows(oddsH1, oddsH2, oddsH3, oddsH4)

odds_all$term<- factor(odds_all$term, levels = c(
  "Inactive|Low_active",
  "Low_active|High_active",
  "diel_periodnight",
  "dist_calc_scaled",
  "dist_calc_scaled:diel_periodnight",
  "botdepth_scaled",
  "botdepth_scaled:diel_periodnight",
  "dist2bot_scaled",
  "dist2bot_scaled:diel_periodnight",
  "depth",
  "depth:diel_periodnight"
))

ggplot(odds_all, aes(x = term, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Comparison of Odds Ratios Across Models",
       y = "Odds Ratio", x = "Term", 
       color = 'Hypothesis') +
  scale_color_manual(values = c("Model 1" = "steelblue", 
                                "Model 2" = "forestgreen", 
                                "Model 3" = "darkorange", 
                                "Model 4" = "purple"),
                     labels = c("Model 1" = "Model 1: Distance to Shore",
                                "Model 2" = "Model 2: Reservoir Bottom Depth",
                                "Model 3" = "Model 3: Distance to Bottom",
                                "Model 4" = "Model 4: Fish Depth")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

