################################################################################
####### Calculating Variance of Observations Across States & Diel Period ####### 
################################################################################

source("./scripts/libraries.R")

#import data 
wels_data<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt.csv")
wels_data$state_3s <- as.factor(wels_data$state_3s)

wels_3s<- wels_data %>% 
  group_by(ID, state_3s) %>% 
  count() %>% 
  pivot_wider(names_from = state_3s, values_from = n) %>% 
  rename("Inactive" = "1",
         "High Activity" = "2", 
         "Low Activity" = "3") 

wels_3s[is.na(wels_3s)] <- 0

wels_3s<- wels_3s %>% 
  column_to_rownames(var = "ID")


#########################################################
##### Cluster Individuals based on Activity Levels ######
#########################################################
library(vegan)
d<- vegdist(wels_3s, method = "bray")
clustering<- hclust(d,method = "complete")
plot(clustering)


#**** Add more variables to the clusters --> inclu. dist2bot, bottom_depth, body_temp, & dist2shore



#########################################################
#### Cluster Based on Habitat Use & Activity Levels #####
#########################################################

detection_depth<- read_csv("./data/detection_depth.csv")
#combine detection data with cleaned HMM dataset
wels_depth<-left_join(wels_data, detection_depth, by = c("ID", "time"))

wels_depth$ID<- as.factor(wels_depth$ID)
wels_depth2<- wels_depth %>% 
  dplyr::select(ID, time, x,y,state_3s, diel_period, temp_mean, state_3s, depth, body_temp)

colnames(wels_depth)
head(wels_depth$depth)
#import + join bottom depth data

bot_depth <- read_csv("./data/pos_mean_wels_filter.csv")
bot_depth<- bot_depth %>% rename(ID = fishid, time = timestamp_5min)
bot_depth2<- bot_depth %>% 
  dplyr::select(ID,time,bottom_depth)

sum(is.na(bot_depth$bottom_depth))

wels_depth3<- left_join(wels_depth2,bot_depth2, by = c("ID", "time"))

colnames(wels_depth3)
library(dplyr)

wels_cluster<- wels_depth3 %>% 
  mutate(day = day(time),
         diel_period = dplyr::recode(diel_period, "day" = 1, "night" = 2)) %>% 
  dplyr::select(ID, day, state_3s, diel_period, temp_mean, state_3s, depth, body_temp, bottom_depth) %>% 
  na.omit()

wels_cluster2<- wels_cluster %>% 
  group_by(ID, state_3s) %>% 
  summarise(cluster_temp = mean(temp_mean), 
         cluster_depth = mean(depth), 
         cluster_bodytemp = mean(body_temp), 
         cluster_botdepth = mean(bottom_depth))

wels_cluster3<- wels_cluster2 %>% 
  group_by(ID, state_3s, cluster_temp, cluster_depth, cluster_bodytemp, cluster_botdepth) %>% 
  count() %>% 
  na.omit() %>% 
  unite(ID_state,ID,state_3s) %>% 
  column_to_rownames(var = "ID_state")

#clustering based on activity & env state avgs
d2<- vegdist(wels_cluster3, method = "euclidean")
clustering2<- hclust(d2,method = "complete")
plot(clustering2)



######
#try pca clustering
######

wels_cluster2 %>% 
 ungroup() %>% 
  group_by(ID)

#pca on activity amt & env vars  
results<- prcomp(wels_cluster3, scale = TRUE)
biplot(results)

# Aggregate data by ID
wels_cluster_agg <- wels_cluster %>% 
  group_by(ID, state_3s) %>% 
  summarise( mean_temp = mean(temp_mean), 
             mean_depth = mean(depth), 
             mean_body_temp = mean(body_temp), 
             mean_bottom_depth = mean(bottom_depth), 
             mean_diel_period = mean(diel_period)) %>% 
  mutate(state = state_3s) %>% 
  unite(ID_state,ID,state_3s) %>% 
  na.omit() %>% 
  column_to_rownames(var = "ID_state") %>% 
  mutate(across(where(is.factor), as.numeric))

# Perform PCA w/ env vars, diel period & state
results3 <- prcomp(wels_cluster_agg, scale = TRUE) 
# Create biplot 
biplot(results3)


#pca with just env vars 
wels_cluster_agg2<- wels_cluster_agg %>% 
  dplyr::select(mean_temp, mean_depth, mean_body_temp, mean_bottom_depth)

results4<- prcomp(wels_cluster_agg2, scale = T)
biplot(results4)

