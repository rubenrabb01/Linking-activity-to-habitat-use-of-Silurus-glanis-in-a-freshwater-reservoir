
#function to compute swimming distance using easting/northing
#input can be changed to use egam/ngam to calculate egam swim distance
comp.dist <- 
  function (easting, northing, depth, timestamp){
    diff.time <- c( NA, diff(as.numeric(timestamp), lag = 1))   # difference in time
    diff.east <- c(NA, diff(easting, lag = 1))                  # difference of easting
    diff.north <- c(NA, diff(northing, lag = 1))                # difference of northing
    diff.depth <- c(NA,diff(depth, lag=1))                      # difference in depth
    swim.dist2D <- sqrt(diff.east^2+diff.north^2)               # swimmed distance in 2D not including depth 
    swim.dist3D <- sqrt(swim.dist2D^2+diff.depth^2)             # swimmed distance in 3D not including depth 
    return(list(diff_depth = diff.depth, diff_time = diff.time, swim_dist_2D = swim.dist2D, swim_dist_3D = swim.dist3D))
  }

