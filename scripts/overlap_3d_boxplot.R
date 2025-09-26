source("./scripts/libraries.R")
source("./scripts/overlap_3d_function.R")


#Reads files for Animal A and Animal B into R workspace, where “Animal a.csv” is a file with X,Y,Z coordinates
wels_3D<- read_csv("./data/wels_filtered_3s_depth.csv")
wels_3D$state_3s<- as.factor(wels_3D$state_3s)

wels_3D %>% 
  na.omit() %>% 
  ggplot(aes(x=state_3s, y = depth, fill = state_3s))+
  geom_boxplot()

wels_3D %>% 
  filter(ID == "T449280_1") %>% 
  ggplot(aes(x=x, y =y, col = state_3s))+
    geom_point()

wels_3D$
unique(wels_3D$ID)

cols_to_select<- c("x","y","depth", "state_3s")

#########################################################################
################### ------ 269 ------- ##################################
#########################################################################
#subset by individual
wels_3D_269 <- wels_3D %>% 
  filter(ID == "T449269_1")

#then subset by state
a_269<-wels_3D_269[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_269<- a_269[,-4]

b_269<-wels_3D_269[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

b_269<- b_269[,-4]

c_269<-wels_3D_269[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_269<- c_269[,-4]

#apply overlap function
overlap_269_ab<- calculate_overlap(a_269,b_269)
overlap_269_ac<- calculate_overlap(a_269,c_269)
overlap_269_bc<- calculate_overlap(b_269,c_269)

#combine just the 1way overlaps together
overlap_269_tot<- rbind(overlap_269_ab$PercentOverlap95a, overlap_269_ac$PercentOverlap95a,overlap_269_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_269_tot<- rbind(overlap_269_ab$Vol95a, overlap_269_ab$Vol95b, overlap_269_ab$Vol50a, overlap_269_ab$Vol50b, overlap_269_ac$Vol95b,overlap_269_ac$Vol50b)


#rename columns + add ID column for both overlap + volume
#assign a = S1, b = S2, c = S3
overlap_269_tot<- overlap_269_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449269_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_269_tot<- vol_269_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449269_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 270 ------- ##################################
#########################################################################
#subset by individual
wels_3D_270 - wels_3D %>% 
  filter(ID == "T449270_1")

#then subset by state
a_270<-wels_3D_270[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_270<- a_270[,-4]

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

#apply overlap function to all combinations of states
overlap_270_ab<- calculate_overlap(a_270,b_270)
overlap_270_ac<- calculate_overlap(a_270,c_270)
overlap_270_bc<- calculate_overlap(b_270,c_270)

#combine just the 1way overlaps together
overlap_270_tot<- rbind(overlap_270_ab$PercentOverlap95a, overlap_270_ac$PercentOverlap95a,overlap_270_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_270_tot<- rbind(overlap_270_ab$Vol95a, overlap_270_ab$Vol95b, overlap_270_ab$Vol50a, overlap_270_ab$Vol50b, overlap_270_ac$Vol95b,overlap_270_ac$Vol50b)

#rename columns + add ID column 
overlap_270_tot<- overlap_270_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449270_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_270_tot<- vol_270_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449270_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 268 ------- ##################################
#########################################################################
#then subset by individual
wels_3D_268 <- wels_3D %>% 
  filter(ID == "T449268_1")

#then subset by state
a_268<-wels_3D_268[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_268<- a_268[,-4]

b_268<-wels_3D_268[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_268<- b_268[,-4]

c_268<-wels_3D_268[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_268<- c_268[,-4]

#apply overlap function to all combinations of states
overlap_268_ab<- calculate_overlap(a_268,b_268)
overlap_268_ac<- calculate_overlap(a_268,c_268)
overlap_268_bc<- calculate_overlap(b_268,c_268)

#combine just the 1way overlaps together
overlap_268_tot<- rbind(overlap_268_ab$PercentOverlap95a, overlap_268_ac$PercentOverlap95a,overlap_268_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_268_tot<- rbind(overlap_268_ab$Vol95a, overlap_268_ab$Vol95b, overlap_268_ab$Vol50a, overlap_268_ab$Vol50b, overlap_268_ac$Vol95b,overlap_268_ac$Vol50b)

#rename columns + add ID column 
overlap_268_tot<- overlap_268_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449268_1", 
         state = c("S1_S2","S1_S3", "S2_S3")) %>% 
  rename(overlap = "5%")

vol_268_tot<- vol_268_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449268_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")

#########################################################################
################### ------ 274 ------- ##################################
#########################################################################
#then subset by individual
wels_3D_274 <- wels_3D %>% 
  filter(ID == "T449274_1")

#then subset by state
a_274<-wels_3D_274[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_274<- a_274[,-4]

b_274<-wels_3D_274[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_274<- b_274[,-4]

c_274<-wels_3D_274[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_274<- c_274[,-4]

#apply overlay function to all combinations of states
overlap_274_ab<- calculate_overlap(a_274,b_274)
overlap_274_ac<- calculate_overlap(a_274,c_274)
overlap_274_bc<- calculate_overlap(b_274,c_274)


#combine just the 1way overlaps together
overlap_274_tot<- rbind(overlap_274_ab$PercentOverlap95a, overlap_274_ac$PercentOverlap95a,overlap_274_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_274_tot<- rbind(overlap_274_ab$Vol95a, overlap_274_ab$Vol95b, overlap_274_ab$Vol50a, overlap_274_ab$Vol50b, overlap_274_ac$Vol95b,overlap_274_ac$Vol50b)

#rename columns + add ID column 
overlap_274_tot<- overlap_274_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449274_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_274_tot<- vol_274_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449274_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 275 ------- ##################################
#########################################################################
#subset by individual
wels_3D_275 <- wels_3D %>% 
  filter(ID == "T449275_1")

#subset by state    
a_275<-wels_3D_275[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_275<- a_275[,-4]

b_275<-wels_3D_275[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_275<- b_275[,-4]

c_275<-wels_3D_275[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_275<- c_275[,-4]

#apply overlay function to all combinations of states
overlap_275_ab<- calculate_overlap(a_275,b_275)
overlap_275_ac<- calculate_overlap(a_275,c_275)
overlap_275_bc<- calculate_overlap(b_275,c_275)


#combine just the 1way overlaps together
overlap_275_tot<- rbind(overlap_275_ab$PercentOverlap95a, overlap_275_ac$PercentOverlap95a,overlap_275_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_275_tot<- rbind(overlap_275_ab$Vol95a, overlap_275_ab$Vol95b, overlap_275_ab$Vol50a, overlap_275_ab$Vol50b, overlap_275_ac$Vol95b,overlap_275_ac$Vol50b)

#rename columns + add ID column 
overlap_275_tot<- overlap_275_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449275_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_275_tot<- vol_275_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449275_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")

#########################################################################
################### ------ 276 ------- ##################################
#########################################################################
#subset by individual
wels_3D_276 <- wels_3D %>% 
  filter(ID == "T449276_1")

#subset by state    
a_276<-wels_3D_276[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_276<- a_276[,-4]

b_276<-wels_3D_276[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_276<- b_276[,-4]

c_276<-wels_3D_276[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_276<- c_276[,-4]

#apply overlay function to all combinations of states
overlap_276_ab<- calculate_overlap(a_276,b_276)
overlap_276_ac<- calculate_overlap(a_276,c_276)
overlap_276_bc<- calculate_overlap(b_276,c_276)

#combine just the 1way overlaps together
overlap_276_tot<- rbind(overlap_276_ab$PercentOverlap95a, overlap_276_ac$PercentOverlap95a,overlap_276_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_276_tot<- rbind(overlap_276_ab$Vol95a, overlap_276_ab$Vol95b, overlap_276_ab$Vol50a, overlap_276_ab$Vol50b, overlap_276_ac$Vol95b,overlap_276_ac$Vol50b)

#rename columns + add ID column 
overlap_276_tot<- overlap_276_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449276_1", 
         state = c("S1_S2","S1_S3", "S2_S3")) %>% 
  rename(overlap = "5%")

vol_276_tot<- vol_276_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449276_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 278 ------- ##################################
#########################################################################
#subset by individual
wels_3D_278 <- wels_3D %>% 
  filter(ID == "T449278_1")

#subset by state    
a_278<-wels_3D_278[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_278<- a_278[,-4]

b_278<-wels_3D_278[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_278<- b_278[,-4]

c_278<-wels_3D_278[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_278<- c_278[,-4]

#apply overlay function to all combinations of states
overlap_278_ab<- calculate_overlap(a_278,b_278)
overlap_278_ac<- calculate_overlap(a_278,c_278)
overlap_278_bc<- calculate_overlap(b_278,c_278)


#combine just the 1way overlaps together
overlap_278_tot<- rbind(overlap_278_ab$PercentOverlap95a, overlap_278_ac$PercentOverlap95a,overlap_278_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_278_tot<- rbind(overlap_278_ab$Vol95a, overlap_278_ab$Vol95b, overlap_278_ab$Vol50a, overlap_278_ab$Vol50b, overlap_278_ac$Vol95b,overlap_278_ac$Vol50b)

#rename columns + add ID column 
overlap_278_tot<- overlap_278_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449278_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_278_tot<- vol_278_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449278_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")



#########################################################################
################### ------ 280 ------- ##################################
#########################################################################

#subset by individual
wels_3D_280 <- wels_3D %>% 
  filter(ID == "T449280_1")

#subset by state    
a_280<-wels_3D_280[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_280<- a_280[,-4]

b_280<-wels_3D_280[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_280<- b_280[,-4]

c_280<-wels_3D_280[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_280<- c_280[,-4]

#apply overlay function to all combinations of states
overlap_280_ab<- calculate_overlap(a_280,b_280)
overlap_280_ac<- calculate_overlap(a_280,c_280)
overlap_280_bc<- calculate_overlap(b_280,c_280)


#combine just the 1way overlaps together
overlap_280_tot<- rbind(overlap_280_ab$PercentOverlap95a, overlap_280_ac$PercentOverlap95a,overlap_280_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_280_tot<- rbind(overlap_280_ab$Vol95a, overlap_280_ab$Vol95b, overlap_280_ab$Vol50a, overlap_280_ab$Vol50b, overlap_280_ac$Vol95b,overlap_280_ac$Vol50b)

#rename columns + add ID column 
overlap_280_tot<- overlap_280_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449280_1", 
         state = c("S1_S2","S1_S3", "S2_S3")) %>% 
  rename(overlap = "5%")

vol_280_tot<- vol_280_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449280_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 282 ------- ##################################
#########################################################################
#subset by individual
wels_3D_282 <- wels_3D %>% 
  filter(ID == "T449282_1")

unique(wels_3D_282$state_3s)
unique(wels_3D_282$depth)

#individual 282 excluded from 3D spatial analysis due to lack of depth data


#########################################################################
################### ------ 286 ------- ##################################
#########################################################################
#subset by individual
wels_3D_286 <- wels_3D %>% 
  filter(ID == "T449286_1")

#subset by state    
a_286<-wels_3D_286[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_286<- a_286[,-4]

b_286<-wels_3D_286[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_286<- b_286[,-4]

c_286<-wels_3D_286[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_286<- c_286[,-4]

#apply overlay function to all combinations of states
overlap_286_ab<- calculate_overlap(a_286,b_286)
overlap_286_ac<- calculate_overlap(a_286,c_286)
overlap_286_bc<- calculate_overlap(b_286,c_286)


#combine just the 1way overlaps together
overlap_286_tot<- rbind(overlap_286_ab$PercentOverlap95a, overlap_286_ac$PercentOverlap95a,overlap_286_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_286_tot<- rbind(overlap_286_ab$Vol95a, overlap_286_ab$Vol95b, overlap_286_ab$Vol50a, overlap_286_ab$Vol50b, overlap_286_ac$Vol95b,overlap_286_ac$Vol50b)

#rename columns + add ID column 
overlap_286_tot<- overlap_286_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449286_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_286_tot<- vol_286_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449286_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 288 ------- ##################################
#########################################################################
#subset by individual
wels_3D_288 <- wels_3D %>% 
  filter(ID == "T449288_1")

#subset by state    
a_288<-wels_3D_288[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_288<- a_288[,-4]

b_288<-wels_3D_288[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_288<- b_288[,-4]

c_288<-wels_3D_288[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_288<- c_288[,-4]

#apply overlay function to all combinations of states
overlap_288_ab<- calculate_overlap(a_288,b_288)
overlap_288_ac<- calculate_overlap(a_288,c_288)
overlap_288_bc<- calculate_overlap(b_288,c_288)


#combine just the 1way overlaps together
overlap_288_tot<- rbind(overlap_288_ab$PercentOverlap95a, overlap_288_ac$PercentOverlap95a,overlap_288_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_288_tot<- rbind(overlap_288_ab$Vol95a, overlap_288_ab$Vol95b, overlap_288_ab$Vol50a, overlap_288_ab$Vol50b, overlap_288_ac$Vol95b,overlap_288_ac$Vol50b)

#rename columns + add ID column 
overlap_288_tot<- overlap_288_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449288_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_288_tot<- vol_288_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449288_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 314 ------- ##################################
#########################################################################
#subset by individual
wels_3D_314 <- wels_3D %>% 
  filter(ID == "T449314_1")

#subset by state    
a_314<-wels_3D_314[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_314<- a_314[,-4]

b_314<-wels_3D_314[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_314<- b_314[,-4]

c_314<-wels_3D_314[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_314<- c_314[,-4]

#apply overlay function to all combinations of states
overlap_314_ab<- calculate_overlap(a_314,b_314)
overlap_314_ac<- calculate_overlap(a_314,c_314)
overlap_314_bc<- calculate_overlap(b_314,c_314)


#combine just the 1way overlaps together
overlap_314_tot<- rbind(overlap_314_ab$PercentOverlap95a, overlap_314_ac$PercentOverlap95a,overlap_314_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_314_tot<- rbind(overlap_314_ab$Vol95a, overlap_314_ab$Vol95b, overlap_314_ab$Vol50a, overlap_314_ab$Vol50b, overlap_314_ac$Vol95b,overlap_314_ac$Vol50b)

#rename columns + add ID column 
overlap_314_tot<- overlap_314_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449314_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_314_tot<- vol_314_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449314_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 318 ------- ##################################
#########################################################################
#subset by individual
wels_3D_318 <- wels_3D %>% 
  filter(ID == "T449318_1")

#subset by state    
a_318<-wels_3D_318[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_318<- a_318[,-4]

b_318<-wels_3D_318[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_318<- b_318[,-4]

c_318<-wels_3D_318[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_318<- c_318[,-4]

#apply overlay function to all combinations of states
overlap_318_ab<- calculate_overlap(a_318,b_318)
overlap_318_ac<- calculate_overlap(a_318,c_318)
overlap_318_bc<- calculate_overlap(b_318,c_318)


#combine just the 1way overlaps together
overlap_318_tot<- rbind(overlap_318_ab$PercentOverlap95a, overlap_318_ac$PercentOverlap95a,overlap_318_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_318_tot<- rbind(overlap_318_ab$Vol95a, overlap_318_ab$Vol95b, overlap_318_ab$Vol50a, overlap_318_ab$Vol50b, overlap_318_ac$Vol95b,overlap_318_ac$Vol50b)

#rename columns + add ID column 
overlap_318_tot<- overlap_318_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449318_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_318_tot<- vol_318_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449318_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
################### ------ 319 ------- ##################################
#########################################################################

#subset by individual
wels_3D_319 <- wels_3D %>% 
  filter(ID == "T449319_1")

#subset by state    
a_319<-wels_3D_319[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 1 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

a_319<- a_319[,-4]

b_319<-wels_3D_319[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 2 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)


b_319<- b_319[,-4]

c_319<-wels_3D_319[cols_to_select] %>% 
  as.data.frame() %>% 
  filter(state_3s == 3 &!is.na(depth)) %>% 
  rename(X = x, 
         Y = y, 
         Z = depth)

c_319<- c_319[,-4]

#apply overlay function to all combinations of states
overlap_319_ab<- calculate_overlap(a_319,b_319)
overlap_319_ac<- calculate_overlap(a_319,c_319)
overlap_319_bc<- calculate_overlap(b_319,c_319)


#combine just the 1way overlaps together
overlap_319_tot<- rbind(overlap_319_ab$PercentOverlap95a, overlap_319_ac$PercentOverlap95a,overlap_319_bc$PercentOverlap95a)

#combine state volumes into 1 dataset
vol_319_tot<- rbind(overlap_319_ab$Vol95a, overlap_319_ab$Vol95b, overlap_319_ab$Vol50a, overlap_319_ab$Vol50b, overlap_319_ac$Vol95b,overlap_319_ac$Vol50b)

#rename columns + add ID column 
overlap_319_tot<- overlap_319_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449319_1", 
         state = c("S1_S2","S1_S3", "S2_S3" )) %>% 
  rename(overlap = "5%")

vol_319_tot<- vol_319_tot %>% 
  as.data.frame() %>% 
  mutate(ID = "T449319_1", 
         state_iso = c("S1_95", "S2_95", "S1_50", "S2_50", "S3_95", "S3_50" )) %>% 
  rename(volume = "5%")


#########################################################################
############ combine all overlay datasets together ######################
#########################################################################

overlap_3d_tot<- rbind(overlap_319_tot, overlap_318_tot, overlap_314_tot,overlap_288_tot,
                    overlap_286_tot,overlap_280_tot,overlap_278_tot,
                    overlap_276_tot,overlap_275_tot,overlap_274_tot,overlap_270_tot, 
                    overlap_269_tot,overlap_268_tot)

volume_tot<- rbind(vol_319_tot, vol_318_tot, vol_314_tot,vol_288_tot,
        vol_286_tot,vol_280_tot,vol_278_tot,
        vol_276_tot,vol_275_tot,vol_274_tot,vol_270_tot, 
        vol_269_tot,vol_268_tot)


overlap_3d_gg<- overlap_3d_tot %>% 
  ggplot(aes(x = state, y = overlap, fill = state))+
  geom_boxplot()+
  theme_bw()+
  xlab("States (2 overlaps 1)")+
  ylab("3D Overlap (%)")+
  labs(caption = "In A_B, overlap is the % of A which is overlapped / covered by B")

volume_50_gg<-volume_tot %>% 
  filter(state_iso %in% c("S1_50","S2_50","S3_50")) %>% 
  ggplot(aes(x = state_iso, y = volume, fill = state_iso))+
  geom_boxplot()+
  theme_bw()+
  xlab("States ")+
  ylab("Volume (m3)")

volume_95_gg<-volume_tot %>% 
  filter(state_iso %in% c("S1_95","S2_95","S3_95")) %>% 
  ggplot(aes(x = state_iso, y = volume, fill = state_iso))+
  geom_boxplot()+
  theme_bw()+
  xlab("States ")+
  ylab("Volume (m3)")

ggarrange(volume_50_gg,volume_95_gg, ncol = 2)
