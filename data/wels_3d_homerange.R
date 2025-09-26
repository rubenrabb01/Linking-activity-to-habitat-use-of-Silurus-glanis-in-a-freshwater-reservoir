#3D homerange



##########################################################################
##########################################################################
##########################################################################

#Loads necessary libraries

# https://academic.oup.com/auk/article/131/4/681/5149381?login=true
# Nathan W. Cooper, Thomas W. Sherry, Peter P. Marra, Modeling three-dimensional space use and overlap in birds, The Auk, Volume 131, Issue 4, 1 October 2014, Pages 681–693, https://doi.org/10.1642/AUK-14-17.1
library(ks)
library(MASS)
library(KernSmooth)
library(plot3D)

#data_3d_t <- data.table(read_csv( "~/3D_track/telemetry_study-3D-master/data/yaps/pike3D_positions_yaps.csv"))
#data_3d_t <- data.table(read_csv( "~/MACFISH movement_visualization/data_for_3Dpaper/position_full_period.csv"))
#data_3d_t[,month := month(timestamp_utc)]


wels_3D<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt.csv")

detection_depth<- read_csv("./data/detection_depth.csv")
#combine detection data with cleaned HMM dataset
wels_depth_3D<-left_join(wels_3D, detection_depth, by = c("ID", "time"))


wels_3D_270 <- wels_depth_3D %>% 
  filter(ID == "T449270_1")
  
  
# Diver vs lazyboy ####
#Reads files for Animal A and Animal B into R workspace, where “Animal a.csv” is a file with X,Y,Z coordinates
test 
#a<-as.data.frame(data_3d_t[nickname == "diver",.(X = x, Y = y, Z = depth)])
#a<-as.data.frame(wels_3D_subset[state_3s == 2 &!is.na(depth),.(X = x, Y = y, Z = depth)])


cols_to_select<- c("x","y","depth", "state_3s")

a_270<-wels_3D_270[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_270<- a_270[,-4]

#b<-as.data.frame(data_3d_t[nickname == "sedentary",.(X = x, Y = y, Z = depth)])
#b<-as.data.frame(wels_3D_subset[state_3s == 3 &!is.na(depth),.(X = x, Y = y, Z = depth)])

b_270<-wels_3D_270[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

b_270<- b_270[,-4]


c_270<-wels_3D_270[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_270<- c_270[,-4]
#Calls the plug-in bandwidth estimator for dataset “a” and “b"
BandA <- Hpi(a)
BandB <- Hpi(b)
#Joins the 2 datasets and then determines the minimum and maximum boundaries for each dimension. It is important #that overlapping territories are evaluated in the same physical space. Also adds small buffer so that no part of the #territories are cut off.

ab_270<-rbind(a_270,b_270)

ac_270<- rbind(a_270,c_270)

bc_270<- rbind(b_270,c_270)

#min / max of ab
minX_ab_270<-min(ab_270$X)-25

minY_ab_270<-min(ab_270$Y)-25

minZ_270<-0

maxX_ab_270<-max(ab_270$X)+25

maxY_ab_270<-max(ab_270$Y)+25

maxZ_ab_270<-max(ab_270$Z)+5

#min / max of ac
minX_a_270c<-min(ac_270$X)-25

minY_ac_270<-min(ac_270$Y)-25

minZ_270<-0

maxX_ac_270<-max(ac_270)+25

maxY_ac_270<-max(ac_270$Y)+25

maxZ_ac_270<-max(ac_270$Z)+5

#min / max of bc
minX_ab_270<-min(bc_270$X)-25

minY_ab_270<-min(bc_270$Y)-25

minZ_270<-0

maxX_bc_270<-max(bc_270$X)+25

maxY_bc_270<-max(bc_270$Y)+25

maxZ_bc_270<-max(bc_270$Z)+5

#Runs the kernel density analysis

fhata_270 <- kde(x=a, H=BandA, binned=FALSE, xmin=c(minX,minY,minZ), xmax=c(maxX,maxY,maxZ),gridsize = 151)

fhatb_270 <- kde(x=b, H=BandB, binned=FALSE, xmin=c(minX,minY,minZ), xmax=c(maxX,maxY,maxZ),gridsize = 151)

fhatc_270 <- kde(x=c, H=BandB, binned=FALSE, xmin=c(minX,minY,minZ), xmax=c(maxX,maxY,maxZ),gridsize = 151)


#Defines 95th isopleth

CL95a <- contourLevels(fhata, cont=95, approx=FALSE)

CL95b <- contourLevels(fhatb, cont=95, approx=FALSE)

CL95c<- contourLevels(fhatc, cont=95, approx=FALSE)

#Calculates the volume at the 95th isopleth

Vol95a<-contourSizes(fhata, cont=95)

Vol95b<-contourSizes(fhatb, cont=95)

Vol95c<-contourSizes(fhatc, cont=95)

#Allows for use of other isopleths when calculating VI and UDOI (see eqn 5)

fhata$estimate <- ifelse(fhata$estimate>=CL95a,fhata$estimate/.95,0)

fhatb$estimate <- ifelse(fhatb$estimate>=CL95b,fhatb$estimate/.95,0)

fhatc$estimate <- ifelse(fhatc$estimate>=CL95c,fhatc$estimate/.95,0)

#Computes individual spatial overlap at 95th isopleth


### a on b, b on a ###
fhat.overlap <- fhata

fhat.overlap$estimate <- fhata$estimate>CL95a & fhatb$estimate>CL95b

#Quantifies volume of overlapped space

Overlap95<-contourSizes(fhat.overlap, abs.cont=.5)

#Computes percentage of territory bird “a” that is overlapped by bird “b” and vice versa

PercentOverlap95ab_270<-(Overlap95/Vol95a)*100

PercentOverlap95ba_270<-(Overlap95/Vol95b)*100


### a on c, c on a ###
fhat.overlap$estimate <- fhata$estimate>CL95a & fhatc$estimate>CL95c

#Quantifies volume of overlapped space

Overlap95<-contourSizes(fhat.overlap, abs.cont=.5)

#Computes percentage of territory bird “a” that is overlapped by bird “b” and vice versa

PercentOverlap95ac_270<-(Overlap95/Vol95a)*100

PercentOverlap95ca_270<-(Overlap95/Vol95c)*100


### c on b, b on c ###
fhat.overlap <- fhatc

fhat.overlap$estimate <- fhatc$estimate>CL95c & fhatb$estimate>CL95b

#Quantifies volume of overlapped space

Overlap95_270<-contourSizes(fhat.overlap, abs.cont=.5)

#Computes percentage of territory bird “a” that is overlapped by bird “b” and vice versa

PercentOverlap95cb_270<-(Overlap95/Vol95c)*100

PercentOverlap95bc_270<-(Overlap95/Vol95b)*100


#Computes voxel size

vol.cell <- prod(sapply(fhata$eval.points, diff)[1,])

#Calculates 3D versions of VI (eqn 2) and UDOI (eqn 4)

VI3D <- sum(pmin(fhata$estimate,fhatb$estimate)*vol.cell)

UDOI3D<-Overlap95*sum(fhata$estimate*fhatb$estimate*vol.cell)



#Plots 95th isopleth of neighboring animals on same graph

plot(fhata,cont=95,col="red",drawpoints=F,xlab="", ylab="", zlab="",size=2, ptcol="black", zlim = c(6,0),
     add=FALSE, box=TRUE, axes=TRUE,
     theta=25, phi=180,  d=3,  xaxt='n')
#plot(fhatb,cont=95,col="green",drawpoints=F,xlab="", ylab="", zlab="",size=2, ptcol="black", zlim = c(6,0),
#     add=T, box=TRUE, axes=TRUE,
#     theta=25, phi=180,  d=3,  xaxt='n')
plot(fhatc,cont=95,col="blue",drawpoints=F,xlab="", ylab="", zlab="",size=2, ptcol="black", zlim = c(6,0),
     add=T, box=TRUE, axes=TRUE,
     theta=25, phi=180,  d=3,  xaxt='n')


#################################################################################
############## Compare 2d vs 3d Overlay #####################################
#################################################################################

#look at all simultaneously
overlap_tot_3d<- rbind(PercentOverlap95ab,PercentOverlap95ba,PercentOverlap95ac,PercentOverlap95ca,PercentOverlap95bc,PercentOverlap95cb)

rownames(overlap_tot_3d) <- c("S1_S3", "S3_S1", "S1_S2", "S2_S1", "S3_S2", "S2_S3")
colnames(overlap_tot_3d) <- "overlap_3d"
overlap_tot_3d<- overlap_tot_3d %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "state_overlay")

#add in 2d overlap data (will have to ensure that overlap_subset is in environment from Rscript wels_2d_homerange.R)
overlap_subset<- overlap_tot %>% 
  filter(ID == "T449270_1" & !overlap == 1) %>% 
  rename(overlap_2d = overlap)

overlay_2d_3d<- merge(overlap_subset,overlap_tot_3d, by = "state_overlay")
overlay_2d_3d<- overlay_2d_3d %>% 
  mutate(overlap_2d = overlap_2d*100)

overlay_table<- nice_table(
  overlay_2d_3d,
  title = c("Table 1", "Comparing 2D and 3D Percent Spatial Overlay"),
    )

flextable::save_as_image(overlay_table, "./figures/overlay_2d_3d_270.png")












#play3d(spin3d(axis = c(0,0,1)))

library(rgl)


#plot3d(fhata,cont=95,col="red",drawpoints=F,xlab="", ylab="", zlab="",size=2, ptcol="black", zlim = c(10,0),
#       add=FALSE, box=TRUE, axes=TRUE,
#       theta=25, phi=195,  d=3,  xaxt='n')


#################################################################################
################ Attempt function to speed up ###################################
#################################################################################

wels_3D<- read_csv("./data/final_model_output_data/HMM_final_output_TimeTempInt.csv")

detection_depth<- read_csv("./data/detection_depth.csv")
#combine detection data with cleaned HMM dataset
wels_depth_3D<-left_join(wels_3D, detection_depth, by = c("ID", "time"))
cols_to_select<- c("x","y","depth", "state_3s")
unique(wels_depth_3D$ID)

wels_3D_270 <- wels_depth_3D %>% 
  filter(ID == "T449270_1")

wels_3D_270<-wels_3D_270[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_268 <- wels_depth_3D %>% 
  filter(ID == "T449268_1")

wels_3D_268<-wels_3D_268[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_269 <- wels_depth_3D %>% 
  filter(ID == "T449269_1")

wels_3D_269<-wels_3D_269[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_274 <- wels_depth_3D %>% 
  filter(ID == "T449274_1")

wels_3D_274<-wels_3D_274[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_275 <- wels_depth_3D %>% 
  filter(ID == "T449275_1")

wels_3D_275<-wels_3D_275[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_276 <- wels_depth_3D %>% 
  filter(ID == "T449276_1")

wels_3D_276<-wels_3D_276[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_278 <- wels_depth_3D %>% 
  filter(ID == "T449278_1")

wels_3D_278<-wels_3D_278[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_280 <- wels_depth_3D %>% 
  filter(ID == "T449280_1")

wels_3D_280<-wels_3D_280[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_282 <- wels_depth_3D %>% 
  filter(ID == "T449282_1")

unique(wels_3D_282$state_3s)

wels_3D_282<-wels_3D_282[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_286 <- wels_depth_3D %>% 
  filter(ID == "T449286_1")

wels_3D_286<-wels_3D_286[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_288 <- wels_depth_3D %>% 
  filter(ID == "T449288_1")

wels_3D_288<-wels_3D_288[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_314 <- wels_depth_3D %>% 
  filter(ID == "T449314_1")

wels_3D_314<-wels_3D_314[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_318 <- wels_depth_3D %>% 
  filter(ID == "T449318_1")

wels_3D_318<-wels_3D_318[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

wels_3D_319 <- wels_depth_3D %>% 
  filter(ID == "T449319_1")

wels_3D_319<-wels_3D_319[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

# Load necessary libraries
library("ks")
library("MASS")
library("KernSmooth")

# Define a function to perform 3D kernel density estimation and spatial overlap analysis, to extract 3d volume & overlap
analyze_movement <- function(data_list) {
  # Initialize a list to store results for all groupings
  results <- list()
  
  for (state in names(data_list)) {
    # Extract data for the current state
    data <- data_list[[state]]
    
    # Estimate bandwidth for the dataset
    bandwidth <- Hpi(data)
    
    # Determine min/max boundaries with a buffer
    minX <- min(data$X) - 25
    minY <- min(data$Y) - 25
    minZ <- 0
    maxX <- max(data$X) + 25
    maxY <- max(data$Y) + 25
    maxZ <- max(data$Z) + 5
    
    # Run kernel density analysis
    kde_result <- kde(x = data, H = bandwidth, binned = FALSE, 
                      xmin = c(minX, minY, minZ), xmax = c(maxX, maxY, maxZ), gridsize = 151)
    
    # Define 95th isopleth
    CL95 <- contourLevels(kde_result, cont = 95, approx = FALSE)
    
    # Volume at the 95th isopleth
    volume_95 <- contourSizes(kde_result, cont = 95)
    
    # Normalize estimate for 95% isopleth
    kde_result$estimate <- ifelse(kde_result$estimate >= CL95, kde_result$estimate / 0.95, 0)
    
    # Save results for the current state
    results[[state]] <- list(
      volume_95 = volume_95,
      kde_result = kde_result
    )
  }
  
  return(results)
}

# Example data structure for all individuals
# Each individual's data is a list of data frames, one for each state (1, 2, "r")
# Pre-loaded data assumed here as `individual_data_list`
# Example: individual_data_list[[1]] = list("1" = df_state1, "2" = df_state2, "r" = df_state_r)

# Loop through individuals and apply the function for each
all_results <- lapply(individual_data_list, analyze_movement)

wels_3D_270_1<- wels_3D_270 %>% 
  filter(state_3s == 1)

wels_3D_270_1 <- wels_3D_270_1[,-4]

wels_3D_270_2<- wels_3D_270 %>% 
  filter(state_3s == 2)

wels_3D_270_2<- wels_3D_270_2[,-4]

wels_3D_270_3<- wels_3D_270 %>% 
  filter(state_3s == 3)

wels_3D_270_3 <- wels_3D_270_3[,-4]



# Create a list for one individual
individual_1 <- list(
  "1" = wels_3D_270_1,  # Data for state 1
  "2" = wels_3D_270_2,  # Data for state 2
  "r" = wels_3D_270_3   # Data for state "r"
)


wels_3D_268_1<- wels_3D_268 %>% 
  filter(state_3s == 1)

wels_3D_268_1 <- wels_3D_268_1[,-4]

wels_3D_268_2<- wels_3D_268 %>% 
  filter(state_3s == 2)

wels_3D_268_2<- wels_3D_268_2[,-4]

wels_3D_268_3<- wels_3D_268 %>% 
  filter(state_3s == 3)

wels_3D_268_3 <- wels_3D_268_3[,-4]


# Create a list for one individual
individual_2 <- list(
  "1" = wels_3D_268_1,  # Data for state 1
  "2" = wels_3D_268_2,  # Data for state 2
  "r" = wels_3D_268_3   # Data for state "r"
)

wels_3D_269_1<- wels_3D_269 %>% 
  filter(state_3s == 1)

wels_3D_269_1 <- wels_3D_269_1[,-4]

wels_3D_269_2<- wels_3D_269 %>% 
  filter(state_3s == 2)

wels_3D_269_2<- wels_3D_269_2[,-4]

wels_3D_269_3<- wels_3D_269 %>% 
  filter(state_3s == 3)

wels_3D_269_3 <- wels_3D_269_3[,-4]


# Create a list for one individual
individual_3 <- list(
  "1" = wels_3D_269_1,  # Data for state 1
  "2" = wels_3D_269_2,  # Data for state 2
  "r" = wels_3D_269_3   # Data for state "r"
)


wels_3D_274_1<- wels_3D_274 %>% 
  filter(state_3s == 1)

wels_3D_274_1 <- wels_3D_274_1[,-4]

wels_3D_274_2<- wels_3D_274 %>% 
  filter(state_3s == 2)

wels_3D_274_2<- wels_3D_274_2[,-4]

wels_3D_274_3<- wels_3D_274 %>% 
  filter(state_3s == 3)

wels_3D_274_3 <- wels_3D_274_3[,-4]


# Create a list for one individual
individual_4 <- list(
  "1" = wels_3D_274_1,  # Data for state 1
  "2" = wels_3D_274_2,  # Data for state 2
  "r" = wels_3D_274_3   # Data for state "r"
)


wels_3D_275_1<- wels_3D_275 %>% 
  filter(state_3s == 1)

wels_3D_275_1 <- wels_3D_275_1[,-4]

wels_3D_275_2<- wels_3D_275 %>% 
  filter(state_3s == 2)

wels_3D_275_2<- wels_3D_275_2[,-4]

wels_3D_275_3<- wels_3D_275 %>% 
  filter(state_3s == 3)

wels_3D_275_3 <- wels_3D_275_3[,-4]


# Create a list for one individual
individual_5 <- list(
  "1" = wels_3D_275_1,  # Data for state 1
  "2" = wels_3D_275_2,  # Data for state 2
  "r" = wels_3D_275_3   # Data for state "r"
)


wels_3D_276_1<- wels_3D_276 %>% 
  filter(state_3s == 1)

wels_3D_276_1 <- wels_3D_276_1[,-4]

wels_3D_276_2<- wels_3D_276 %>% 
  filter(state_3s == 2)

wels_3D_276_2<- wels_3D_276_2[,-4]

wels_3D_276_3<- wels_3D_276 %>% 
  filter(state_3s == 3)

wels_3D_276_3 <- wels_3D_276_3[,-4]


# Create a list for one individual
individual_6 <- list(
  "1" = wels_3D_276_1,  # Data for state 1
  "2" = wels_3D_276_2,  # Data for state 2
  "r" = wels_3D_276_3   # Data for state "r"
)


wels_3D_278_1<- wels_3D_278 %>% 
  filter(state_3s == 1)

wels_3D_278_1 <- wels_3D_278_1[,-4]

wels_3D_278_2<- wels_3D_278 %>% 
  filter(state_3s == 2)

wels_3D_278_2<- wels_3D_278_2[,-4]

wels_3D_278_3<- wels_3D_278 %>% 
  filter(state_3s == 3)

wels_3D_278_3 <- wels_3D_278_3[,-4]


# Create a list for one individual
individual_7 <- list(
  "1" = wels_3D_278_1,  # Data for state 1
  "2" = wels_3D_278_2,  # Data for state 2
  "r" = wels_3D_278_3   # Data for state "r"
)


wels_3D_280_1<- wels_3D_280 %>% 
  filter(state_3s == 1)

wels_3D_280_1 <- wels_3D_280_1[,-4]

wels_3D_280_2<- wels_3D_280 %>% 
  filter(state_3s == 2)

wels_3D_280_2<- wels_3D_280_2[,-4]

wels_3D_280_3<- wels_3D_280 %>% 
  filter(state_3s == 3)

wels_3D_280_3 <- wels_3D_280_3[,-4]


# Create a list for one individual
individual_8 <- list(
  "1" = wels_3D_280_1,  # Data for state 1
  "2" = wels_3D_280_2,  # Data for state 2
  "r" = wels_3D_280_3   # Data for state "r"
)


wels_3D_282_1<- wels_3D_282 %>% 
  filter(state_3s == 1)

wels_3D_282_1 <- wels_3D_282_1[,-4]

wels_3D_282_2<- wels_3D_282 %>% 
  filter(state_3s == 2)

wels_3D_282_2<- wels_3D_282_2[,-4]

wels_3D_282_3<- wels_3D_282 %>% 
  filter(state_3s == 3)

wels_3D_282_3 <- wels_3D_282_3[,-4]


# Create a list for one individual
individual_9 <- list(
  "1" = wels_3D_282_1,  # Data for state 1
  "2" = wels_3D_282_2,  # Data for state 2
  "r" = wels_3D_282_3   # Data for state "r"
)

#Individual 9 does not have depth data, hence was removed from analysis

wels_3D_286_1<- wels_3D_286 %>% 
  filter(state_3s == 1)

wels_3D_286_1 <- wels_3D_286_1[,-4]

wels_3D_286_2<- wels_3D_286 %>% 
  filter(state_3s == 2)

wels_3D_286_2<- wels_3D_286_2[,-4]

wels_3D_286_3<- wels_3D_286 %>% 
  filter(state_3s == 3)

wels_3D_286_3 <- wels_3D_286_3[,-4]


# Create a list for one individual
individual_10 <- list(
  "1" = wels_3D_286_1,  # Data for state 1
  "2" = wels_3D_286_2,  # Data for state 2
  "r" = wels_3D_286_3   # Data for state "r"
)


wels_3D_288_1<- wels_3D_288 %>% 
  filter(state_3s == 1)

wels_3D_288_1 <- wels_3D_288_1[,-4]

wels_3D_288_2<- wels_3D_288 %>% 
  filter(state_3s == 2)

wels_3D_288_2<- wels_3D_288_2[,-4]

wels_3D_288_3<- wels_3D_288 %>% 
  filter(state_3s == 3)

wels_3D_288_3 <- wels_3D_288_3[,-4]


# Create a list for one individual
individual_11 <- list(
  "1" = wels_3D_288_1,  # Data for state 1
  "2" = wels_3D_288_2,  # Data for state 2
  "r" = wels_3D_288_3   # Data for state "r"
)


wels_3D_314_1<- wels_3D_314 %>% 
  filter(state_3s == 1)

wels_3D_314_1 <- wels_3D_314_1[,-4]

wels_3D_314_2<- wels_3D_314 %>% 
  filter(state_3s == 2)

wels_3D_314_2<- wels_3D_314_2[,-4]

wels_3D_314_3<- wels_3D_314 %>% 
  filter(state_3s == 3)

wels_3D_314_3 <- wels_3D_314_3[,-4]


# Create a list for one individual
individual_12 <- list(
  "1" = wels_3D_314_1,  # Data for state 1
  "2" = wels_3D_314_2,  # Data for state 2
  "r" = wels_3D_314_3   # Data for state "r"
)


wels_3D_318_1<- wels_3D_318 %>% 
  filter(state_3s == 1)

wels_3D_318_1 <- wels_3D_318_1[,-4]

wels_3D_318_2<- wels_3D_318 %>% 
  filter(state_3s == 2)

wels_3D_318_2<- wels_3D_318_2[,-4]

wels_3D_318_3<- wels_3D_318 %>% 
  filter(state_3s == 3)

wels_3D_318_3 <- wels_3D_318_3[,-4]


# Create a list for one individual
individual_13 <- list(
  "1" = wels_3D_318_1,  # Data for state 1
  "2" = wels_3D_318_2,  # Data for state 2
  "r" = wels_3D_318_3   # Data for state "r"
)


wels_3D_319_1<- wels_3D_319 %>% 
  filter(state_3s == 1)

wels_3D_319_1 <- wels_3D_319_1[,-4]

wels_3D_319_2<- wels_3D_319 %>% 
  filter(state_3s == 2)

wels_3D_319_2<- wels_3D_319_2[,-4]

wels_3D_319_3<- wels_3D_319 %>% 
  filter(state_3s == 3)

wels_3D_319_3 <- wels_3D_319_3[,-4]


# Create a list for one individual
individual_14 <- list(
  "1" = wels_3D_319_1,  # Data for state 1
  "2" = wels_3D_319_2,  # Data for state 2
  "r" = wels_3D_319_3   # Data for state "r"
)


# Combine all individuals into a single list
individual_data_list <- list(
  "individual_1" = individual_1,
  "individual_2" = individual_2,
  "individual_3" = individual_3,
  "individual_4" = individual_4,
  "individual_5" = individual_5,
  "individual_6" = individual_6,
  "individual_7" = individual_7,
  "individual_8" = individual_8,
 # "individual_9" = individual_9,
  "individual_10" = individual_10,
  "individual_11" = individual_11,
  "individual_12" = individual_12,
  "individual_13" = individual_13,
  "individual_14" = individual_14
  # Add more individuals as needed
)

#Individual 9 does not have depth data, hence was removed from analysis


# Loop through individuals and apply the function for each
all_results <- lapply(individual_data_list, analyze_movement)



#View Volume of all individuals
all_results#$individual_1
for (i in seq_along(all_results)) {
  cat("\nVolumes for Individual", i, ":\n")
  if (!is.null(all_results[[i]])) {
    for (state in names(all_results[[i]])) {
      if (!is.null(all_results[[i]][[state]]$volume_95)) {
        volume <- all_results[[i]][[state]]$volume_95
        cat("  State", state, ": Volume_95 =", volume, "\n")
      } else {
        cat("  State", state, ": Volume_95 not available.\n")
      }
    }
  } else {
    cat("  No results available for this individual.\n")
  }
}


#View Overlap of all
# Initialize an empty list to store overlap data
overlap_results <- list()

# Loop through individuals and extract overlaps
for (i in seq_along(all_results)) {
  # Check if results for the individual exist
  if (!is.null(all_results[[i]])) {
    # Initialize a matrix to store overlaps for the current individual
    state_names <- names(all_results[[i]])
    n_states <- length(state_names)
    overlap_matrix <- matrix(NA, nrow = n_states, ncol = n_states,
                             dimnames = list(state_names, state_names))
    
    # Loop through state pairs to fill the overlap matrix
    for (state1 in state_names) {
      for (state2 in state_names) {
        if (state1 != state2) {
          overlap_matrix[state1, state2] <- all_results[[i]][[state1]]$overlap[[state2]]
        }
      }
    }
    
    # Store the overlap matrix for this individual
    overlap_results[[i]] <- overlap_matrix
  }
}

# View overlap results for a specific individual
overlap_results[[1]]  # Overlap matrix for individual 1
