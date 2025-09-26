source("./scripts/libraries.R")
source("./scripts/data_load.R")
#cleaning out outliers (points outside of rimov polygon)

pos_mean<- read_csv("./data/pos_mean2.csv")
rimov_pol<-st_read("./data/shp_files")

#convert fish positions to sf point
pos_mean_sf<- st_as_sf(pos_mean, coords = c("easting", "northing"))

#check reference systems
st_crs(pos_mean_sf)
st_crs(rimov_pol)

#merge fish positions to match rimov CRS
pos_mean_sf<-pos_mean_sf %>% 
  st_set_crs(32633)

#remove points which fall outside of rimov polygon (gps errors)
test<-st_intersects(pos_mean_sf,rimov_pol, sparse = TRUE)
pos_mean_sf$within_poly<-apply(test, 1, any) #show true / false values, and create a column in original sf dataset

nrow(pos_mean_sf)

#filter for positions that are within polygon
pos_mean_filter<- pos_mean_sf %>% 
  filter(within_poly == TRUE) 

#convert back to easting/northing data (not sure why i did it this way but it works :D) 
pos_mean_filter<- as.data.frame(st_coordinates(pos_mean_sf$geometry))
pos_mean_filter <- cbind(pos_mean_sf[,-c(1)],pos_mean_filter)

colnames(pos_mean_filter)[1] ="fishid"

pos_mean_filter <- pos_mean_filter %>%
  arrange(fishid, timestamp_5min) %>% 
  filter(within_poly == TRUE) 

head(pos_mean_filter)
nrow(pos_mean_filter)
#export filtered data to csv
write_csv(pos_mean_filter, "pos_mean_wels_filter.csv")


################################################
####### visualize filtered data ################
################################################

#plot rimov polygon, and before + after filtering
plot(st_geometry(rimov_pol), axes = TRUE, ylim = c(5407000,5408000)) 
plot(st_geometry(pos_mean_sf), add = TRUE, pch= 19, col = "red")
plot(st_geometry(pos_mean_filter), add = TRUE, pch= 19, col = "lightblue")


#add 10m buffer from rimov shore
#ggplot() +
#  geom_sf(data = st_intersection(pos_mean_filter, st_buffer(rimov_pol, 10)))