#visualizing 3d state overlap 


library(rKIN)

#rim<- read_csv("./data/wels_filtered_3s_covars.csv")
#rim<- rim %>% 
#  dplyr::select(ID, step,angle,x,y,time,diel_period,state_3s,body_temp, bottom_depth)

rim<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt_UPDATED.csv")
rim2<- rim %>% 
      dplyr::select(ID, step,angle,x,y,time,dp,state_3s)


                 
bot_depth <- read_csv("./data/pos_mean_wels_filter.csv")
bot_depth2<- bot_depth %>% rename(ID = fishid, time = timestamp_5min)

rim_depth<- left_join(rim2,bot_depth2, by = c("ID", "time"))#, "bottom_depth"))

rim_depth_test <- rim_depth

#############################################
## Estimate missing values in bottom depth ##
#############################################

rim_depth_test2 <- rim_depth_test %>% 
  filter(ID != "T449282_1")

# Start by making a copy and converting to data.table
rim_depth_test2 <- data.table(rim_depth_test2)

# Create timestamp column from 'time'
rim_depth_test2[, timestamp := time]

#Estimate missing depth values by fish ID and with no timing gaps later than 45 minutes
# Make sure data is sorted by timestamp
setorder(rim_depth_test2, timestamp)

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
rim_depth_test3 <- interpolate_column(rim_depth_test2, "ID", "depth")
rim_depth_test4 <- interpolate_column(rim_depth_test3, "ID", "bottom_depth")

#calculate distance to bottom 
rim_depth_test5<-rim_depth_test4 %>% 
  mutate(dist2bot = (bottom_depth - depth))

head(rim_depth_test5)

rim_depth_test3 <- as_tibble(rim_depth_test2) #data with non-interpolated depth + bottom depth values
rim_depth_test6 <- as_tibble(rim_depth_test5) #data with estimated depths 
#changing the diel period classification has extremely large influence on output visuals

rim_depth2<- rim_depth_test6 %>% 
  na.omit(depth, bottom_depth) %>% 
  filter(bottom_depth>0) %>% 
  dplyr::select(depth,bottom_depth,state_3s, dp) %>% 
  as.data.frame() #must be in dataframe to work

hval_set <- c(2,2)
rim_depth2$depth<- -rim_depth2$depth


####For comparing day vs night --> will have to make 2 separate datasets for day vs night###

rim_depth2_day <- rim_depth2 %>% 
  filter(dp == "day")#filter(diel_period == "day")#
rim_depth2_night <- rim_depth2 %>% 
  filter(dp == "night") #filter(diel_period == "night")#

# import rimov thermocline
rimov_th <- data.table(read_csv("./data/rimov_thermo.csv"))
rimov_th <- rimov_th[ date >= "2017-06-01" & date <= "2017-08-31", .(date, thermo.depth) ]
mean(rimov_th$thermo.depth)
range(rimov_th$thermo.depth)



library(ggrepel)
library(scales)
library(ggpubr)

S1_col <- "brown3"
S2_col <- "cadetblue3"
S3_col <- "darkolivegreen"
pol_col <- "grey"
lab_state<- c("Inactivity", "High activity", "Low activity")
alpha_vec <- c(0.8,0.4,0.4 ) #c(0.6,0.3,0.1 )
ratio_size <- 4


############################################################################################
############ DAY VS NIGHT ####################
  
    
# calculation of overlay day v night
rim.kin_day<- estKIN(data=rim_depth2_day, x="bottom_depth", y="depth", group="state_3s", levels=c(50,95), hval = hval_set)
# Extract the area of each polygon
rim.area_day<- getArea(rim.kin_day)

rim.kin_night<- estKIN(data=rim_depth2_night, x="bottom_depth", y="depth", group="state_3s", levels=c(50,95), hval = hval_set)
# Extract the area of each polygon
rim.area_night<- getArea(rim.kin_night)

# determine polygon overlap for all polygons
rim.olp_day<- calcOverlap(rim.kin_day)
rim_line_day <- fortify(rim.kin_day$estObj)
rim_line_day <- st_as_sf(rim_line_day)

rim.olp_night<- calcOverlap(rim.kin_night)
rim_line_night <- fortify(rim.kin_night$estObj)
rim_line_night <- st_as_sf(rim_line_night)

#############################################
### create table of overlap for core area ###
#############################################
#day table
overlap_2d_day<- rim.olp_day %>% 
  pivot_longer(!OverlapID, names_to = "overlapping_state", values_to = "percent_overlap")

overlap_2d_day<- overlap_2d_day %>% 
  rename("Overlapped_state" = OverlapID)

#night
overlap_2d_night<- rim.olp_night %>% 
  pivot_longer(!OverlapID, names_to = "overlapping_state", values_to = "percent_overlap")

overlap_2d_night<- overlap_2d_night %>% 
  rename("Overlapped_state" = OverlapID)

overlap_2d_day2 <- overlap_2d_day %>% 
  mutate(dp = "day")

overlap_2d_night2 <- overlap_2d_night %>% 
  mutate(dp = "night")

overlap_2d<- full_join(overlap_2d_day2,overlap_2d_night2, by = c("Overlapped_state", "overlapping_state", "percent_overlap", "dp"))

overlap_2d_a <- overlap_2d %>% 
  mutate(percent_overlap = percent_overlap*100)

###create table of overlap for core area (50%)
overlap_50_tab <-overlap_2d_a %>%
  filter(str_detect(Overlapped_state, "_50") & str_detect(overlapping_state, "_50")) %>% 
  mutate(
    Overlapped_state = recode_factor(Overlapped_state,
                                     '1_50' = "Inactivity",
                                     '2_50' = "High Activity",
                                     '3_50' = "Low Activity"),
    overlapping_state = recode_factor(overlapping_state,
                                      '1_50' = "Inactivity",
                                      '2_50' = "High Activity",
                                      '3_50' = "Low Activity")
  )

overlap_50_tab2 <- overlap_50_tab %>% 
  filter(Overlapped_state != overlapping_state)

# Group the table by diel_period
overlap_50_tab_grouped <- overlap_50_tab2 %>%
  spread(key = dp, value = percent_overlap)

# Now create a table with the custom header and grouped by diel_period
tab_df(overlap_50_tab_grouped, digits = 1,
       col.header = c("Reference State", "Overlapping State", "Day Period", "Night Period"),title = "Percent Overlap of Core Area (50% Confidence Interval)")

#########################################
### create table for Home Range (95%) ###
#########################################

overlap_95_tab <-overlap_2d_a %>%
  filter(str_detect(Overlapped_state, "_95") & str_detect(overlapping_state, "_95")) %>% 
  mutate(
    Overlapped_state = recode_factor(Overlapped_state,
                                     '1_95' = "Inactivity",
                                     '2_95' = "High Activity",
                                     '3_95' = "Low Activity"),
    overlapping_state = recode_factor(overlapping_state,
                                      '1_95' = "Inactivity",
                                      '2_95' = "High Activity",
                                      '3_95' = "Low Activity")
  )

overlap_95_tab2 <- overlap_95_tab %>% 
  filter(Overlapped_state != overlapping_state)

# Group the table by diel_period
overlap_95_tab_grouped <- overlap_95_tab2 %>%
  spread(key = dp, value = percent_overlap)

# Now create a table with the custom header and grouped by diel_period
tab_df(overlap_95_tab_grouped,digits = 1,
       col.header = c("Reference State", "Overlapping State", "Day Period", "Night Period"),title = "Percent Overlap of Home Range (95% Confidence Interval)")



######### VISUALIZE #################

library(scales)
library(ggpubr)

###day###

# Define colors with explicit mapping to Group values
group_colors1 <- c(
  "1_50" = "brown3",       # Inactivity
  "2_50" = "cadetblue3",   # High activity
  "3_50" = "darkolivegreen" # Low activity
)

# Define labels explicitly
group_labels1 <- c(
  "1_50" = "Inactivity",
  "2_50" = "High activity",
  "3_50" = "Low activity"
)


# Define colors with explicit mapping to Group values
group_colors2 <- c(
  "1_95" = "brown3",       # Inactivity
  "2_95" = "cadetblue3",   # High activity
  "3_95" = "darkolivegreen" # Low activity
)

# Define labels explicitly
group_labels2 <- c(
  "1_95" = "Inactivity",
  "2_95" = "High activity",
  "3_95" = "Low activity"
)

rim50_gg_day <- rim_line_day %>%
  mutate(Group_ConfInt = factor(paste(Group, ConfInt, sep = "_"),
         levels = c("2_50", "2_95", "3_50", "3_95", "1_50", "1_95"))) %>%
  filter(ConfInt == 50) %>% 
  arrange(Group) %>%
  ggplot() +
  geom_sf(aes(fill = Group_ConfInt)) +
  scale_fill_manual("State",
                    values = alpha(group_colors1,alpha_vec),
                    labels = group_labels1,
                    na.value = "grey") +  # Ensures NA values are handled
  scale_y_continuous(limits = c(-16, 0), expand = c(0, 0), labels = abs) +
  scale_x_continuous(limits = c(0, 45), expand = c(0, 0), breaks = seq(0,45,5)) +
  #ggtitle("Spatial Overlay during Day") +
  theme_bw(base_size = 20) +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.8, 0.4),
        legend.direction = "vertical",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Legend label text size
        legend.title = element_text(size = 20),     # Legend title siz
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(y = "Depth (m)", fill = "States") + # subtitle = "Core Area (50% Confidence Interval)"
  geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "#333333", linewidth = 1, linetype = "dashed")

rim50_gg_day_polygon <- rim50_gg_day + 
  annotate(geom = 'polygon', x = c(0, 0, 15), y = c(-16, 0, -16), 
           color = "#666666", fill = "#666666")


rim95_gg_day <- rim_line_day %>%
  mutate(Group_ConfInt = factor(paste(Group, ConfInt, sep = "_"),
                                levels = c("2_50", "2_95", "3_50", "3_95", "1_50", "1_95")))%>%
  filter(ConfInt==95) %>% 
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill =  Group_ConfInt))+
  scale_fill_manual("State",
                    values = alpha(group_colors2,alpha_vec),
                    labels = group_labels2,
                    na.value = "grey") + 
  scale_y_continuous(limits = c(-16, 0), expand = c(0, 0), labels = abs)+
  scale_x_continuous(limits = c(0, 45), expand = c(0, 0), breaks = seq(0,45,5))+
  #ggtitle("Spatial Overlay during Day")+
  theme_bw(base_size = 20)+
  theme(legend.position = "inside", 
        legend.position.inside=c(0.8,0.4),
        legend.direction = "vertical",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Legend label text size
        legend.title = element_text(size = 20),     # Legend title siz
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(color = "States")+#, subtitle = "Home Range (95% Confidence Interval)")+
  #facet_wrap(~diel_period)+
  geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "#333333", linewidth  = 1, linetype = "dashed")

rim95_gg_day_polygon <- rim95_gg_day + annotate(geom = 'polygon', x = c(0,0,15), y= c(-16,0,-16), 
                                                color = "#666666", fill = "#666666")#alpha=0.2, 



###night###
rim50_gg_night <- rim_line_night %>%
  mutate(Group_ConfInt = factor(paste(Group, ConfInt, sep = "_"),
                                levels = c("2_50", "2_95", "3_50", "3_95", "1_50", "1_95"))) %>%
  filter(ConfInt == 50) %>% 
  arrange(Group) %>%
  ggplot() +
  geom_sf(aes(fill = Group_ConfInt)) +
  scale_fill_manual("State",
                    values = alpha(group_colors1,alpha_vec),
                    labels = group_labels1,
                    na.value = "grey") +  # Ensures NA values are handled
  scale_y_continuous(limits = c(-16, 0), expand = c(0, 0), labels = abs) +
  scale_x_continuous(limits = c(0, 45),expand = c(0, 0), breaks = seq(0,45,5)) +
  #ggtitle("Spatial Overlay during Night") +
  theme_bw(base_size = 20) +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.8, 0.4),
        legend.direction = "vertical",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Legend label text size
        legend.title = element_text(size = 20),     # Legend title siz
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "States")+ #(x = "Bottom Depth (m)", y = "Depth (m)", fill = "States")+#,subtitle = "Core Area (50% Confidence Interval)") +
  geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "#333333", linewidth = 1, linetype = "dashed")

rim50_gg_night_polygon <- rim50_gg_night + 
  annotate(geom = 'polygon', x = c(0, 0, 15), y = c(-16, 0, -16), 
           color = "#666666", fill = "#666666")


rim95_gg_night <- rim_line_night %>%
  mutate(Group_ConfInt = factor(paste(Group, ConfInt, sep = "_"),
                                levels = c("2_50", "2_95", "3_50", "3_95", "1_50", "1_95")))%>%
  filter(ConfInt==95) %>% 
  arrange(Group)%>%
  ggplot() +
  #geom_hline(yintercept = 0, col = "black", linewidth  = 0.3, linetype = "solid")+
  geom_sf(aes(fill =  Group_ConfInt))+
  scale_fill_manual("State",
                    values = alpha(group_colors2,alpha_vec),
                    labels = group_labels2,
                    na.value = "grey") + 
  scale_y_continuous(limits = c(-16, 0), expand = c(0, 0), labels = abs)+
  scale_x_continuous(limits = c(0, 45), expand = c(0, 0), breaks = seq(0,45,5))+
  #ggtitle("Spatial Overlay during Night")+
  theme_bw(base_size = 20)+
  theme(legend.position = "inside", 
        legend.position.inside=c(0.8,0.4),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.text = element_text(size = 20),     # Legend label text size
        legend.title = element_text(size = 20),     # Legend title siz
        legend.direction = "vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(fill = "State")+#x = "Bottom Depth (m)", color = "States", subtitle = "Home Range (95% Confidence Interval)")+
  #facet_wrap(~diel_period)+
  geom_hline(yintercept = -mean(rimov_th$thermo.depth), col = "#333333", linewidth  = 1, linetype = "dashed")

rim95_gg_night_polygon <- rim95_gg_night + annotate(geom = 'polygon', x = c(0,0,15), y= c(-16,0,-16), 
                                                color = "#666666", fill = "#666666")


# Arrange the plots

combined_plot <- ggarrange(rim50_gg_day_polygon,rim95_gg_day_polygon, rim50_gg_night_polygon, rim95_gg_night_polygon, 
                            nrow = 2,ncol = 2, 
                            align = "hv", 
                           labels = c("a)", "b)", "c)", "d)"),
                            common.legend = T, 
          legend = 'top')

# Add common y-axis label
final_plot <- annotate_figure(combined_plot,
                              left = text_grob("Depth (m)", rot = 90, hjust = 1, size = 15),
                              bottom = text_grob("Bottom Depth (m)", size = 15))

# Print the final plot
print(final_plot)

################################################################################

# calcualte max reservoir + fish depth for each state for each fig 
library(sf)
library(dplyr)

# Define a function that correctly extracts values, feature-by-feature
extract_max_values <- function(sf_object, diel_period) {
  
  sf_filtered <- sf_object %>%
    st_cast("MULTIPOLYGON")   # ensure clean geometry
  
  # Create a data frame: one row per feature
  results <- lapply(seq_len(nrow(sf_filtered)), function(i) {
    
    coords <- st_coordinates(sf_filtered[i,])
    
    data.frame(
      Group = sf_filtered$Group[i],
      ConfInt = sf_filtered$ConfInt[i],
      max_depth = min(coords[,2]),
      max_resDepth = max(coords[,1]),
      Diel = diel_period
    )
  }) %>%
    bind_rows()   # Combine into one dataframe
  
  return(results)
}

# Apply function to day and night datasets
day_results <- extract_max_values(rim_line_day, "Day")
night_results <- extract_max_values(rim_line_night, "Night")

# Combine
rim_summary <- bind_rows(day_results, night_results) %>%
  arrange(Group, Diel, ConfInt)

# View
rim_summary
