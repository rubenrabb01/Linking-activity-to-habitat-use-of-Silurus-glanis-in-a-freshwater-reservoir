library(raster)
library(adehabitatHR)
library(ggspatial)
# 2d home range ####
shape.rim.utm <- st_read("~/catfish_rimov/data/shp_files/shoreline/Rimov_pol_469m_UTM33.shp")

shape_points_rim <- data.frame(st_coordinates(shape.rim.utm)[,c(1,2)])
names(shape_points_rim) <- c("x", "y")
shape_points_rim <- shape_points_rim[1:(length(shape_points_rim$x)),]

# simplified lake, needed for kernel density estimate
rim.simple <- st_read("~/catfish_rimov/data/shp_files/simple_shoreline/rimov_pol_simple2_UTM33.shp")
# shape file buffered in order make it larger, necessary for function "getverticeshr"
rim.buffer <- st_buffer(shape.rim.utm, dist = 500)
raster.rim <- raster(rim.buffer,resolution = 10)
pixel.rim.b<-SpatialPixels(SpatialPoints(coordinates(raster.rim)))

# conversion of simpled shapes to appropriate class
sh.sim.rim<-as_Spatial(st_cast(rim.simple,"LINESTRING"))


wels_3s_covar<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt_UPDATED.csv")
unique(wels_3s_covar$ID)

###############################################################################################################
###############################################################################################################

wels_268 <- wels_3s_covar %>% 
  filter(ID == "T449268_1")
sp_268 <- st_as_sf(wels_268, coords = c("x","y"),crs = 32633)
sp_268_in <- st_filter(sp_268, rim.simple)
sp_268_in <- as_Spatial(sp_268_in)

# Cut only positions within reservoir shape
kud_268<- kernelUD(sp_268_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_268 <- getverticeshr(kud_268, 95)
area_268 <- kernel.area(kud_268,percent=95) # area in hectares
area_268$ID <- unique(area_268$ID)
overlap_268 <- kerneloverlaphr(kud_268, meth="HR", conditional=TRUE)
line_268_sf <- st_as_sf(line_268)

#ID is character, change to factor and specify levels 
line_268_sf$id<-as.factor(line_268_sf$id)

state3_268<- line_268_sf %>% 
  filter(id == 2 )
state2_268 <- line_268_sf %>% 
  filter(id == 3)
state1_268 <- line_268_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg268 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_268, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_268, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_268, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  annotation_scale(location = "bl")+
  ggtitle("268")+
  theme(#legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

gg268a <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_268, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_268, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_268, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  annotation_scale(location = "bl", text_cex = 1.5)+
  ggtitle("268")+
  theme(#legend.position='none',
    panel.border = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size = 13)
    
  )


#############################################################################################
###############################################################################################################
wels_269 <- wels_3s_covar %>% 
  filter(ID == "T449269_1")
sp_269 <- st_as_sf(wels_269, coords = c("x","y"),crs = 32633)
sp_269_in <- st_filter(sp_269, rim.simple)
sp_269_in <- as_Spatial(sp_269_in)

# Cut only positions within reservoir shape
kud_269<- kernelUD(sp_269_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_269 <- getverticeshr(kud_269, 95)
area_269 <- kernel.area(kud_269,percent=95) # area in hectares
area_269$ID <- unique(area_269$ID)
overlap_269 <- kerneloverlaphr(kud_269, meth="HR", conditional=TRUE)
line_269_sf <- st_as_sf(line_269)

#ID is character, change to factor and specify levels 
line_269_sf$id<-as.factor(line_269_sf$id)

state3_269<- line_269_sf %>% 
  filter(id == 2 )
state2_269 <- line_269_sf %>% 
  filter(id == 3)
state1_269 <- line_269_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg269 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_269, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_269, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_269, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("269")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

##################################################################
###############################################################################################################
wels_270 <- wels_3s_covar %>% 
  filter(ID == "T449270_1")
sp_270 <- st_as_sf(wels_270, coords = c("x","y"),crs = 32633)
sp_270_in <- st_filter(sp_270, rim.simple)
sp_270_in <- as_Spatial(sp_270_in)

# Cut only positions within reservoir shape
kud_270<- kernelUD(sp_270_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_270 <- getverticeshr(kud_270, 95)
area_270 <- kernel.area(kud_270,percent=95) # area in hectares
area_270$ID <- unique(area_270$ID)
overlap_270 <- kerneloverlaphr(kud_270, meth="HR", conditional=TRUE)
line_270_sf <- st_as_sf(line_270)

#ID is character, change to factor and specify levels 
line_270_sf$id<-as.factor(line_270_sf$id)

state3_270<- line_270_sf %>% 
  filter(id == 2 )
state2_270 <- line_270_sf %>% 
  filter(id == 3)
state1_270 <- line_270_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg270 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_270, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_270, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_270, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("270")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################
wels_274 <- wels_3s_covar %>% 
  filter(ID == "T449274_1")
sp_274 <- st_as_sf(wels_274, coords = c("x","y"),crs = 32633)
sp_274_in <- st_filter(sp_274, rim.simple)
sp_274_in <- as_Spatial(sp_274_in)

# Cut only positions within reservoir shape
kud_274<- kernelUD(sp_274_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_274 <- getverticeshr(kud_274, 95)
area_274 <- kernel.area(kud_274,percent=95) # area in hectares
area_274$ID <- unique(area_274$ID)
overlap_274 <- kerneloverlaphr(kud_274, meth="HR", conditional=TRUE)
line_274_sf <- st_as_sf(line_274)

#ID is character, change to factor and specify levels 
line_274_sf$id<-as.factor(line_274_sf$id)

state3_274<- line_274_sf %>% 
  filter(id == 2 )
state2_274 <- line_274_sf %>% 
  filter(id == 3)
state1_274 <- line_274_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg274 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_274, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_274, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_274, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("274")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )




###############################################################################################################
###############################################################################################################

wels_275 <- wels_3s_covar %>% 
  filter(ID == "T449275_1")
sp_275 <- st_as_sf(wels_275, coords = c("x","y"),crs = 32633)
sp_275_in <- st_filter(sp_275, rim.simple)
sp_275_in <- as_Spatial(sp_275_in)

# Cut only positions within reservoir shape
kud_275<- kernelUD(sp_275_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_275 <- getverticeshr(kud_275, 95)
area_275 <- kernel.area(kud_275,percent=95) # area in hectares
area_275$ID <- unique(area_275$ID)
overlap_275 <- kerneloverlaphr(kud_275, meth="HR", conditional=TRUE)
line_275_sf <- st_as_sf(line_275)

#ID is character, change to factor and make subsets of each state
line_275_sf$id<-as.factor(line_275_sf$id)

state3_275<- line_275_sf %>% 
  filter(id == 2 )
state2_275 <- line_275_sf %>% 
  filter(id == 3)
state1_275 <- line_275_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg275 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_275, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_275, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_275, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("275")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################


wels_276 <- wels_3s_covar %>% 
  filter(ID == "T449276_1")
sp_276 <- st_as_sf(wels_276, coords = c("x","y"),crs = 32633)
sp_276_in <- st_filter(sp_276, rim.simple)
sp_276_in <- as_Spatial(sp_276_in)

# Cut only positions within reservoir shape
kud_276<- kernelUD(sp_276_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_276 <- getverticeshr(kud_276, 95)
area_276 <- kernel.area(kud_276,percent=95) # area in hectares
area_276$ID <- unique(area_276$ID)
overlap_276 <- kerneloverlaphr(kud_276, meth="HR", conditional=TRUE)
line_276_sf <- st_as_sf(line_276)

#ID is character, change to factor and specify levels 
line_276_sf$id<-as.factor(line_276_sf$id)

state3_276<- line_276_sf %>% 
  filter(id == 2 )
state2_276 <- line_276_sf %>% 
  filter(id == 3)
state1_276 <- line_276_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg276a <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_276, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_276, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_276, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl", text_cex = 1.5)+
  ggtitle("276")+
  theme(#legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 13)
        
  )

gg276 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_276, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_276, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_276, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl", text_cex = 1.5)+
  ggtitle("276")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )
###############################################################################################################
###############################################################################################################


wels_278 <- wels_3s_covar %>% 
  filter(ID == "T449278_1")
sp_278 <- st_as_sf(wels_278, coords = c("x","y"),crs = 32633)
sp_278_in <- st_filter(sp_278, rim.simple)
sp_278_in <- as_Spatial(sp_278_in)

# Cut only positions within reservoir shape
kud_278<- kernelUD(sp_278_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_278 <- getverticeshr(kud_278, 95)
area_278 <- kernel.area(kud_278,percent=95) # area in hectares
area_278$ID <- unique(area_278$ID)
overlap_278 <- kerneloverlaphr(kud_278, meth="HR", conditional=TRUE)
line_278_sf <- st_as_sf(line_278)

#is is character, change to factor and specify levels 
line_278_sf$id<-as.factor(line_278_sf$id)

state3_278<- line_278_sf %>% 
  filter(id == 2 )
state2_278 <- line_278_sf %>% 
  filter(id == 3)
state1_278 <- line_278_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg278 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_278, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_278, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_278, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("278")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################
wels_280 <- wels_3s_covar %>% 
  filter(ID == "T449280_1")
sp_280 <- st_as_sf(wels_280, coords = c("x","y"),crs = 32633)
sp_280_in <- st_filter(sp_280, rim.simple)
sp_280_in <- as_Spatial(sp_280_in)

# Cut only positions within reservoir shape
kud_280<- kernelUD(sp_280_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_280 <- getverticeshr(kud_280, 95)
area_280 <- kernel.area(kud_280,percent=95) # area in hectares
area_280$ID <- unique(area_280$ID)
overlap_280 <- kerneloverlaphr(kud_280, meth="HR", conditional=TRUE)
line_280_sf <- st_as_sf(line_280)


#is is character, change to factor and specify levels 
line_280_sf$id<-as.factor(line_280_sf$id)

state3_280<- line_280_sf %>% 
  filter(id == 2 )
state2_280 <- line_280_sf %>% 
  filter(id == 3)
state1_280 <- line_280_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg280a <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_280, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_280, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_280, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("280")+
  theme(#legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 13)
        
  )

#visualize, plot in order of stacked appearance
gg280 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_280, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_280, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_280, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("280")+
  theme(#legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )





###############################################################################################################
###############################################################################################################


wels_282 <- wels_3s_covar %>% 
  filter(ID == "T449282_1")
sp_282 <- st_as_sf(wels_282, coords = c("x","y"),crs = 32633)
sp_282_in <- st_filter(sp_282, rim.simple)
sp_282_in <- as_Spatial(sp_282_in)

# Cut only positions within reservoir shape
kud_282<- kernelUD(sp_282_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_282 <- getverticeshr(kud_282, 95)
area_282 <- kernel.area(kud_282,percent=95) # area in hectares
area_282$ID <- unique(area_282$ID)
overlap_282 <- kerneloverlaphr(kud_282, meth="HR", conditional=TRUE)
line_282_sf <- st_as_sf(line_282)

#is is character, change to factor and specify levels 
line_282_sf$id<-as.factor(line_282_sf$id)

state3_282<- line_282_sf %>% 
  filter(id == 2 )
state2_282 <- line_282_sf %>% 
  filter(id == 3)
state1_282 <- line_282_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg282 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_282, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_282, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_282, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("282")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )


###############################################################################################################
###############################################################################################################
wels_286 <- wels_3s_covar %>% 
  filter(ID == "T449286_1")
sp_286 <- st_as_sf(wels_286, coords = c("x","y"),crs = 32633)
sp_286_in <- st_filter(sp_286, rim.simple)
sp_286_in <- as_Spatial(sp_286_in)

# Cut only positions within reservoir shape
kud_286<- kernelUD(sp_286_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_286 <- getverticeshr(kud_286, 95)
area_286 <- kernel.area(kud_286,percent=95) # area in hectares
area_286$ID <- unique(area_286$ID)
overlap_286 <- kerneloverlaphr(kud_286, meth="HR", conditional=TRUE)
line_286_sf <- st_as_sf(line_286)

#check class of variable    
class(line_286_sf$id) 



#BIA_type is character, change to factor and specify levels 
line_286_sf$id<-as.factor(line_286_sf$id)

#reorder features so migration comes first  --> didnt work to solve issue
state_select_ordered <- line_286_sf %>% 
  arrange(id = c(1,3,2))

state3_286<- line_286_sf %>% 
  filter(id == 2 )
state2_286 <- line_286_sf %>% 
  filter(id == 3)
state1_286 <- line_286_sf %>% 
  filter(id == 1)




#visualize
gg286 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_286, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_286, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_286, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("286")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################
wels_288 <- wels_3s_covar %>% 
  filter(ID == "T449288_1")
sp_288 <- st_as_sf(wels_288, coords = c("x","y"),crs = 32633)
sp_288_in <- st_filter(sp_288, rim.simple)
sp_288_in <- as_Spatial(sp_288_in)

# Cut only positions within reservoir shape
kud_288<- kernelUD(sp_288_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_288 <- getverticeshr(kud_288, 95)
area_288 <- kernel.area(kud_288,percent=95) # area in hectares
area_288$ID <- unique(area_288$ID)
overlap_288 <- kerneloverlaphr(kud_288, meth="HR", conditional=TRUE)
line_288_sf <- st_as_sf(line_288)
vol_288<-  getvolumeUD(kud_288)

#is is character, change to factor and specify levels 
line_288_sf$id<-as.factor(line_288_sf$id)

state3_288<- line_288_sf %>% 
  filter(id == 2 )
state2_288 <- line_288_sf %>% 
  filter(id == 3)
state1_288 <- line_288_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg288 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_288, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_288, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_288, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("288")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )



###############################################################################################################
###############################################################################################################

wels_314 <- wels_3s_covar %>% 
  filter(ID == "T449314_1")
sp_314 <- st_as_sf(wels_314, coords = c("x","y"),crs = 32633)
sp_314_in <- st_filter(sp_314, rim.simple)
sp_314_in <- as_Spatial(sp_314_in)

# Cut only positions within reservoir shape
kud_314<- kernelUD(sp_314_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_314 <- getverticeshr(kud_314, 95)
area_314 <- kernel.area(kud_314,percent=95) # area in hectares
area_314$ID <- unique(area_314$ID)
overlap_314 <- kerneloverlaphr(kud_314, meth="HR", conditional=TRUE)
line_314_sf <- st_as_sf(line_314)


#is is character, change to factor and specify levels 
line_314_sf$id<-as.factor(line_314_sf$id)

state3_314<- line_314_sf %>% 
  filter(id == 2 )
state2_314 <- line_314_sf %>% 
  filter(id == 3)
state1_314 <- line_314_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg314 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_314, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_314, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_314, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("314")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################
wels_318 <- wels_3s_covar %>% 
  filter(ID == "T449318_1")
sp_318 <- st_as_sf(wels_318, coords = c("x","y"),crs = 32633)
sp_318_in <- st_filter(sp_318, rim.simple)
sp_318_in <- as_Spatial(sp_318_in)

# Cut only positions within reservoir shape
kud_318<- kernelUD(sp_318_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_318 <- getverticeshr(kud_318, 95)
area_318 <- kernel.area(kud_318,percent=95) # area in hectares
area_318$ID <- unique(area_318$ID)
overlap_318 <- kerneloverlaphr(kud_318, meth="HR", conditional=TRUE)
line_318_sf <- st_as_sf(line_318)

#is is character, change to factor and specify levels 
line_318_sf$id<-as.factor(line_318_sf$id)

state3_318<- line_318_sf %>% 
  filter(id == 2 )
state2_318 <- line_318_sf %>% 
  filter(id == 3)
state1_318 <- line_318_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg318 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_318, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_318, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_318, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "bl")+
  ggtitle("318")+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )
###############################################################################################################
###############################################################################################################
wels_319 <- wels_3s_covar %>% 
  filter(ID == "T449319_1")
sp_319 <- st_as_sf(wels_319, coords = c("x","y"),crs = 32633)
sp_319_in <- st_filter(sp_319, rim.simple)# Cut only positions within reservoir shape
sp_319_in <- as_Spatial(sp_319_in)

#calculate kernel density
kud_319<- kernelUD(sp_319_in[,"state_3s"], h=50 , grid=pixel.rim.b, boundary = sh.sim.rim)  # h=50 seemed to give good results accosrding to visual inspection
line_319 <- getverticeshr(kud_319, 95)
area_319 <- kernel.area(kud_319,percent=95) # area in hectares
area_319$ID <- unique(area_319$ID)
overlap_319 <- kerneloverlaphr(kud_319, meth="HR", conditional=TRUE)
line_319_sf <- st_as_sf(line_319)

#is is character, change to factor and specify levels 
line_319_sf$id<-as.factor(line_319_sf$id)

state3_319<- line_319_sf %>% 
  filter(id == 2 )
state2_319 <- line_319_sf %>% 
  filter(id == 3)
state1_319 <- line_319_sf %>% 
  filter(id == 1)


#visualize, plot in order of stacked appearance
gg319 <- ggplot() +
  geom_sf(data = shape.rim.utm, fill =NA ) +
  geom_sf(data = state3_319, aes(col = id, fill = id), color = "cadetblue3",fill = "cadetblue3",alpha = 0.5)+
  geom_sf(data = state2_319, aes(col = id, fill = id), color = "darkolivegreen",fill = "darkolivegreen",alpha = 0.5)+
  geom_sf(data = state1_319, aes(col = id, fill = id), color = "brown1",fill = "brown1",alpha = 0.5)+
  #geom_sf(data = swim_sp, size = 0.01, alpha  = 0.3, col = "black")+
  #annotate("text", x = 404400, y = 5601200, label = paste("2D HR",round(area_swim[,1],0), "ha",sep = " "), size = 7, col = "red")+
  #annotation_scale(location = "tr")+
  ggtitle("319", )+
  theme(legend.position='none',
        panel.border = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)
        
  )

###############################################################################################################
###############################################################################################################

unique(wels_3s_covar$ID)

ggarrange(gg268,gg269,gg270,gg274,gg275,gg276,gg278,gg280,gg286,gg288, gg314,gg318,gg319, nrow = 5, ncol = 3) #removed gg282 because only 2 states detected


ggarrange(gg268a,gg276a,gg280a, ncol = 3, nrow = 1, common.legend = T)


ggarrange(gg269,gg270,gg274, ncol = 3, nrow = 1, common.legend = T)

ggarrange(gg275,gg278,gg286, ncol = 3, nrow = 1, common.legend = T)

ggarrange(gg288,gg314,gg319, ncol = 3, nrow = 1, common.legend = T)


###############################################################################################################
############################## Overall States Area Boxplot #######################################
###############################################################################################################

area_268<- area_268 %>% 
  mutate(ID = "T449268_1")
area_269<- area_269 %>% 
  mutate(ID = "T449269_1")
area_270<- area_270 %>% 
  mutate(ID = "T449270_1")
area_274<- area_274 %>% 
  mutate(ID = "T449274_1")
area_276<- area_276 %>% 
  mutate(ID = "T449276_1")
area_275<- area_275 %>% 
  mutate(ID = "T449275_1")
area_278<- area_278 %>% 
  mutate(ID = "T449278_1")
area_280<- area_280 %>% 
  mutate(ID = "T449280_1")
area_288<- area_288 %>% 
  mutate(ID = "T449288_1")
area_286<- area_286 %>% 
  mutate(ID = "T449286_1")
area_314<- area_314 %>% 
  mutate(ID = "T449314_1")
area_318<- area_318 %>% 
  mutate(ID = "T449318_1")
area_319<- area_319 %>% 
  mutate(ID = "T449319_1")

area_282 <-  area_282 %>% 
  mutate(X2 = NA, 
         ID = "T449282_1")

area_282<- area_282 %>% 
  relocate(X3, .after = X2)

#combine rows together
area_tot<- rbind(area_268,area_269,area_270,area_274,area_276,area_275,area_278,area_280,
                 area_286,area_288,area_314,area_318,area_319)#,area_282

area_tot_df<- tibble(area_tot)

area_tot_df1<- area_tot_df %>% 
  rename(S1 = X1, S2 = X3, S3 = X2) #this is the correct order of activity states

area_tot_df2<- area_tot_df1 %>% 
pivot_longer(cols = c(S1,S2,S3), names_to = 'state') %>% 
  rename(area = value)

library(EnvStats)

colorder<- c("Inactivity", "Low Activity", "High Activity")

area_tot_df3<- area_tot_df2 %>% 
  mutate(state = recode_factor(state, S1 = "Inactivity", S2 = "Low Activity", S3 = "High Activity"))



area_tot_df3 %>% 
  filter(ID != "T449282_1") %>% 
  ggplot(aes(x = reorder(state,area), y = area, fill = state))+
  geom_boxplot()+
  scale_fill_manual(values = c("brown3","darkolivegreen","cadetblue3"), guide = "none")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,125,25))+
  labs(x = "State", y = "2D Area (Ha)")+
  stat_n_text(geom = "text", y.pos = 125) #from EnvStats

area_tot_df3$area <- as.numeric(area_tot_df3$area)

area_anova <- aov(area ~ state, data = area_tot_df3)

summary(area_anova)

TukeyHSD(area_anova)

#check anova assumptions
shapiro.test(residuals(area_anova))#non sig, assumption met 
leveneTest(area ~ state, data = area_tot_df3) #sig, assumption not met

leveneTest(log(area) ~ state, data = area_tot_df3)

#perform test with logged area
area_anova2 <- aov(log(area) ~ state, data = area_tot_df3)
summary(area_anova2)

TukeyHSD(area_anova2) #all sig
shapiro.test(residuals(area_anova2))#non sig, assumption met 


area_tot_df3 %>% 
  filter(ID != "T449282_1") %>% 
  ggplot(aes(x = reorder(state,area), y = log(area), fill = state))+
  geom_boxplot()+
  #geom_violin(adjust = 0.5, draw_quantiles = 0.5)+
  scale_fill_manual(values = c("brown3","darkolivegreen","cadetblue3"), guide = "none")+
  theme_bw(base_size = 15)+
  #scale_y_continuous(breaks = seq(0,125,25))+
  labs(x = "State", y = "Log Area (Ha)")+
  #stat_n_text(geom = "text", y.pos = 0.5)+ #from EnvStats
  stat_compare_means(
    method = "t.test", 
    comparisons = list(c("Inactivity", "Low Activity"), 
                       c("Inactivity", "High Activity"), 
                       c("Low Activity", "High Activity")),
    label = "p.signif",
    label.y = c(4.8, 5.1, 5.5)
  )+
  theme(#axis.title = element_text(size = 13),
        #axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = 'none')  


  

# Calculate mean and standard error
summary_stats <- area_tot_df3 %>%
  filter(ID != "T449282_1") %>% 
  group_by(state) %>% 
  na.omit(area) %>% 
  summarize(
    mean = mean(area),
    se = sd(area) / sqrt(n())
  )


tab_df(summary_stats,
       col.header = c("State", "Mean Area (Ha)", "SE"))

print(summary_stats) #TOTAL mean+SE OF STATE 2D AREA in ha ?



hist(sqrt(area_tot_df$area)) #best transformation so far for all states together 

area_tot_df %>% 
  filter(state == "S1") %>% 
  ggplot(aes(x=sqrt(area)))+
  geom_histogram(bins = 10)


state_area_aov<- aov(sqrt(area)~state, data = area_tot_df)
summary(state_area_aov)
TukeyHSD(state_area_aov) #all states are sig. diff. from each other

#visualize assumptions
plot(state_area_aov)

#testing assumptions

#normality of residuals
#by groups 
area_tot_df %>% 
  group_by(state) %>% 
  mutate(sqr_area = sqrt(area)) %>% 
  shapiro_test(sqr_area)
#non sig = normality met

#overall
shapiro_test(residuals(state_area_aov))
#non sig = normality met

#homogeneity of variances
area_tot_df %>% 
  mutate(sqr_area = sqrt(area), 
         state = as.factor(state)) %>%
  levene_test(sqr_area~state)
#non-sig = homogeneity of variance met (also looked good in plotting)
#only sightly non sig. tho... why??

####################################################################

#check HR area b/w states by percentile
hrs <- mcp.area(area_tot_df[,3], percent=seq(50, 100, by = 5))
?mcp.area





#######Check difference in home range based on percentiles
hrs_test <- mcp.area(sp_268_in[,"state_3s"], percent = seq(50, 100, by = 5))

#varies a lot particularly b/w 95-100 in state 1 + State 2 (cols 1, 3) --> but not the same for all individuals

###############################################################################################################
##############################   Overall States Area Overlap    ###############################################
###############################################################################################################



colnames(overlap_268) <- c("S1", "S3", "S2")
rownames(overlap_268) <- c("S1", "S3", "S2")

overlap_268_df<- overlap_268 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449268_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_268_long<- overlap_268_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_268_long<- overlap_268_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)



#re-run overlap code from above across all individuals to not duplicate columns / rows
colnames(overlap_269) <- c("S1", "S3", "S2")
rownames(overlap_269) <- c("S1", "S3", "S2")

overlap_269_df<- overlap_269 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449269_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_269_long<- overlap_269_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_269_long<- overlap_269_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#######################################################################
colnames(overlap_270) <- c("S1", "S3", "S2")
rownames(overlap_270) <- c("S1", "S3", "S2")

overlap_270_df<- overlap_270 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449270_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_270_long<- overlap_270_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_270_long<- overlap_270_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#################################################################
colnames(overlap_274) <- c("S1", "S3", "S2")
rownames(overlap_274) <- c("S1", "S3", "S2")

overlap_274_df<- overlap_274 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449274_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_274_long<- overlap_274_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_274_long<- overlap_274_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
####################################################################
colnames(overlap_276) <- c("S1", "S3", "S2")
rownames(overlap_276) <- c("S1", "S3", "S2")

overlap_276_df<- overlap_276 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449276_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_276_long<- overlap_276_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_276_long<- overlap_276_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#########################################################################
colnames(overlap_275) <- c("S1", "S3", "S2")
rownames(overlap_275) <- c("S1", "S3", "S2")

overlap_275_df<- overlap_275 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449275_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_275_long<- overlap_275_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_275_long<- overlap_275_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
###########################################################################
colnames(overlap_278) <- c("S1", "S3", "S2")
rownames(overlap_278) <- c("S1", "S3", "S2")

overlap_278_df<- overlap_278 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449278_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_278_long<- overlap_278_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_278_long<- overlap_278_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#########################################################################
colnames(overlap_280) <- c("S1", "S3", "S2")
rownames(overlap_280) <- c("S1", "S3", "S2")

overlap_280_df<- overlap_280 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449280_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_280_long<- overlap_280_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_280_long<- overlap_280_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
######################################################################
colnames(overlap_288) <- c("S1", "S3", "S2")
rownames(overlap_288) <- c("S1", "S3", "S2")

overlap_288_df<- overlap_288 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449288_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_288_long<- overlap_288_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_288_long<- overlap_288_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#######################################################################
colnames(overlap_286) <- c("S1", "S3", "S2")
rownames(overlap_286) <- c("S1", "S3", "S2")

overlap_286_df<- overlap_286 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449286_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_286_long<- overlap_286_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_286_long<- overlap_286_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
######################################################################
colnames(overlap_314) <- c("S1", "S3", "S2")
rownames(overlap_314) <- c("S1", "S3", "S2")

overlap_314_df<- overlap_314 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449314_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_314_long<- overlap_314_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_314_long<- overlap_314_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
#######################################################################
colnames(overlap_318) <- c("S1", "S3", "S2")
rownames(overlap_318) <- c("S1", "S3", "S2")

overlap_318_df<- overlap_318 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449318_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_318_long<- overlap_318_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_318_long<- overlap_318_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
########################################################################
colnames(overlap_319) <- c("S1", "S3", "S2")
rownames(overlap_319) <- c("S1", "S3", "S2")

overlap_319_df<- overlap_319 %>% 
  as.data.frame() %>% 
  mutate(ID = "T449319_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_319_long<- overlap_319_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_319_long<- overlap_319_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)

###########################################################################################
overlap_282<- rbind(overlap_282, c(NA,NA))
overlap_282<- cbind(overlap_282, c(NA,NA,NA))

colnames(overlap_282) <- c("S1", "S3", "S2")
rownames(overlap_282) <- c("S1", "S3", "S2")


overlap_282_df<- overlap_282 %>% 
  as.data.frame() %>% 
  mutate( ID = "T449282_1") %>%  #add a column with the fish ID
  rownames_to_column(var = "under_state") #create a column with the row names

overlap_282_long<- overlap_282_df %>% 
  pivot_longer(cols = c(-under_state, -ID), #pivot all columns except under_state & ID
               names_to = "over_state", values_to = "overlap") #put the col names into a col called "over_state", and the values into "overlap" col

overlap_282_long<- overlap_282_long %>% 
  mutate(state_overlay= paste(under_state, over_state, sep = "_")) %>% 
  dplyr::select(ID, state_overlay, overlap)
###############################################################################################################
###############################################################################################################

#combine all datasets
overlap_tot<-rbind(overlap_268_long,overlap_269_long,overlap_270_long,overlap_274_long,overlap_275_long,overlap_276_long,overlap_278_long,overlap_280_long,overlap_286_long,overlap_288_long,overlap_314_long,overlap_318_long,overlap_319_long,overlap_282_long)

#visualize... not sure this is the best way though...
overlap_tot %>% 
  ggplot(aes(x = state_overlay, y = overlap, fill = state_overlay))+
  geom_boxplot()+
  theme_bw()+
  xlab("States (2 overlaps 1)")+
  ylab("% overlap (m2)")+
  ylim(c(0,1))+
  stat_n_text(geom = "text", y.pos = 0) #from EnvStats package


overlap_mean_2d<- overlap_tot %>% 
  filter(ID != "T449282_1") %>% 
  na.omit() %>% 
  group_by(state_overlay) %>% 
  summarise_at(vars(overlap), list(mean_overlap = mean))

overlap_tot %>% 
  filter(ID == "T449270_1")



# Calculate mean and standard error
summary_stats_overlap <- overlap_tot %>%
  filter(ID != "T449282_1",
         !state_overlay %in% c("S1_S1","S2_S2","S3_S3")) %>% 
  group_by(state_overlay ) %>% 
  na.omit(overlap) %>% 
  summarize(
    mean = mean(overlap),
    se = sd(overlap) / sqrt(n())
  )

summary_stats_overlap2 <- summary_stats_overlap %>% 
  separate(state_overlay, sep = "_", into = c("Reference_State", "Overlapping_State"))%>%
  mutate(mean = mean*100,
         se = se*100,
    Reference_State = recode_factor(Reference_State,
                             "S1" = "Inactivity",
                             "S2" = "Low Activity",
                             "S3" = "High Activity"),
    Overlapping_State = recode_factor(Overlapping_State,
                               "S1" = "Inactivity",
                               "S2" = "Low Activity",
                               "S3" = "High Activity"))

tab_df(summary_stats_overlap2,
       col.header = c("Reference State", "Overlapping State", "Mean Overlap (%)", "SE"))

###############################################################################################################
###############################################################################################################


################# Work on for loop later #####################
ID_list<- unique(wels_3s_covar$ID)
list_of_df<- list()

for (i in length(ID_list)){
    x<- wels_3s_covar %>% 
      filter(ID == ID_list[i])
    #run test
    sp <- st_as_sf(x, coords = c("x","y"),crs = 32633)
    sp_in <- st_filter(sp, rim.simple)
    sp_in <- as_Spatial(sp_in)
    
    return(list_of_df)
  }
########################################################

 
 
 
 
 
 
 
 
 
 
