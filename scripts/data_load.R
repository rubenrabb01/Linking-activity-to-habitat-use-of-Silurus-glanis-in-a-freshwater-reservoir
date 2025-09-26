# Script for catfish report 2017
# showing different aspect of catfish behaviour in Rimov reservoir
library(RPostgreSQL)
library(data.table)
#library(maptools)#not available - archive install did not work
library(lubridate)
library(rgdal)
#library(postGIStools) #not available - archive install did not work
library(devtools)

# connect to database
con <-  dbConnect(drv = PostgreSQL(), dbname ="teridb", host="172.21.3.20", user= "teriuser", password = "t3r1us3r!")

#select species
species <- "wels"

# make selection
select.fish.qu <- 
  paste("SELECT fi_fishid
FROM teri.fish
WHERE
fi_species = '", species,"' ;", sep = "")

# getting tag numbers
tagid.list <- data.table(dbGetQuery(con,select.fish.qu ))
tagid.list <- tagid.list[,fi_fishid]

# query for positions
select.position.qu <- paste("WITH tmp as (SELECT  * FROM teri.positions 
                            where 
                            up_fishvalid AND
                            up_timestamp_utc BETWEEN '2017-07-01 00:00:00' AND '2017-08-31 23:59:59' AND 
                            tu_tagmode = 'fish' AND
                            fi_fishid IN ('", paste(tagid.list, collapse = "','"),"')
                            )

                            SELECT up_timestamp_utc, b.lt_lake, up_easting, up_northing, ST_X(ST_Transform(up_geom, 4326)) as lon,ST_Y(ST_Transform(up_geom, 4326)) as lat, up_egam, up_ngam, up_depth, up_validpos, fi_fishid, up_gamres, up_fishvalid,  up_bottom_depth,
                                   ST_LineLocatePoint(lcl_centerline, up_geom_mid) * ST_length(lcl_centerline) distfromdam 
                                FROM tmp b INNER JOIN teri.lakecenterline a ON b.lt_lake = a.lt_lake ", sep = "")


# upload positions from database
positions <- data.table(dbGetQuery(con, select.position.qu))
#positions[, week := week(up_timestamp_utc)]
#positions[, hour := hour(up_timestamp_utc)] --> had to remove because it double counted columns on every hour
#positions[, date := as.Date(up_timestamp_utc)] --> had to remove the rest because more duplicate columns
#positions[, month := month(up_timestamp_utc)]
positions[, timestamp_5min := round_date(up_timestamp_utc, "5 minutes")]

positions <- positions[up_gamres < 100, ]

pos_mean <- positions[,.(easting = mean(up_easting), northing = mean(up_northing),
                                        depth = mean(up_depth, na.rm = T), bottom_depth = mean(up_bottom_depth)),
                         by = .(fi_fishid, timestamp_5min)]

#save as csv to working directory
write.csv(pos_mean,"pos_mean2.csv") #renamed to pos_mean2 to separate from original dataset we had downloaded which included hour calculation

head(pos_mean)
#########################################################################################
############## Download depth + body temperature detection data ##################################
#########################################################################################

library(RPostgreSQL)
library(data.table)
library(lubridate)

# open connection
con <-  dbConnect(drv = PostgreSQL(), dbname ="teridb", host="172.21.3.20", user= "teriuser", password = "t3r1us3r!")
# define species
species <- "wels"

# define time interval
start.time <- "2017-07-01 00:00:00"
end.time <-  "2017-08-31 23:59:59"

# make selection
select.fish.qu <-
  paste("SELECT fi_fishid
        FROM teri.fish
        WHERE
        fi_species = '", species,"' ;", sep = "")

# getting tag numbers
tagid.list <- data.table(dbGetQuery(con,select.fish.qu ))
tagid.list <- tagid.list[,fi_fishid]


# getting detections from single receivers in Velesinska bay (rec no. 1500117), tributary (1500)
select.detections.qu <- paste("
                              SELECT  dd_timestamp_utc, ht_hsn, fi_fishid, dd_depth, dd_bodytemp FROM teri.detsdepth_rec_fish
                              where
                              NOT dead AND
                              dd_timestamp_utc BETWEEN ","'",start.time,"'", "AND","'", end.time,"'", " AND
                              fi_fishid IN ('", paste(tagid.list, collapse = "','"),"');", sep = "")

detection <- data.table(dbGetQuery(con, select.detections.qu ))
# round_time
detection[, ts_5min := round_date(dd_timestamp_utc, " 5 mins")]
# calculate mena depth and body temperture for 5 minutes interval
detection_depth <- detection[,.(depth = mean(dd_depth, na.rm = T), body_temp = mean(dd_bodytemp, na.rm = T)), by = .(fi_fishid, ts_5min)]

summary(detection_depth)
summary(HMM_wels_clean)

detection_depth$fi_fishid<- as.factor(detection_depth$fi_fishid)
detection_depth<- detection_depth %>% 
  rename( ID = fi_fishid,
          time = ts_5min)


#save as csv to working directory
write.csv(detection_depth,"./data/detection_depth.csv")


head(detection_depth)
