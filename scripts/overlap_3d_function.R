library(ks) # for Hpi, kde functions

# Define a function that takes three datasets and processes them
calculate_overlap<- function(dataset_a, dataset_b, buffer_xy = 25, buffer_z = 5, grid_size = 151) {
  # Calculate plug-in bandwidth estimators for both datasets
  BandA <- Hpi(dataset_a)
  BandB <- Hpi(dataset_b)
  
  # Combine the datasets to determine boundaries
  ab <- rbind(dataset_a, dataset_b)
  # Determine min and max boundaries for each dimension (with buffer)
  minX <- min(ab$X) - buffer_xy
  minY <- min(ab$Y) - buffer_xy
  minZ <- 0
  maxX <- max(ab$X) + buffer_xy
  maxY <- max(ab$Y) + buffer_xy
  maxZ <- max(ab$Z) + buffer_z
  
  # Run kernel density analysis
  fhata <- kde(x = dataset_a, H = BandA, binned = FALSE, xmin = c(minX, minY, minZ), xmax = c(maxX, maxY, maxZ), gridsize = grid_size)
  fhatb <- kde(x = dataset_b, H = BandB, binned = FALSE, xmin = c(minX, minY, minZ), xmax = c(maxX, maxY, maxZ), gridsize = grid_size)
  
  # Define 95th isopleth
  CL95a <- contourLevels(fhata, cont = 95, approx = FALSE)
  CL95b <- contourLevels(fhatb, cont = 95, approx = FALSE)
  
  # Calculate the volume at the 95th isopleth
  Vol95a <- contourSizes(fhata, cont = 95)
  Vol95b <- contourSizes(fhatb, cont = 95)
  
  # Calculate the volume at the 50th isopleth
  Vol50a <- contourSizes(fhata, cont = 50)
  Vol50b <- contourSizes(fhatb, cont = 50)
  
  # Adjust estimates to reflect 95th isopleth
  fhata$estimate <- ifelse(fhata$estimate >= CL95a, fhata$estimate / .95, 0)
  fhatb$estimate <- ifelse(fhatb$estimate >= CL95b, fhatb$estimate / .95, 0)
  
  # Compute individual spatial overlap at 95th isopleth
  fhat.overlap <- fhata
  fhat.overlap$estimate <- fhata$estimate > CL95a & fhatb$estimate > CL95b
  
  # Quantify the volume of overlapped space
  Overlap95 <- contourSizes(fhat.overlap, abs.cont = 0.5)
  
  # Compute percentage of territory overlapped
  PercentOverlap95a <- (Overlap95 / Vol95a) * 100
  PercentOverlap95b <- (Overlap95 / Vol95b) * 100
  
  # Compute voxel size
  vol.cell <- prod(sapply(fhata$eval.points, diff)[1,])
  
  # Calculate 3D versions of VI and UDOI
  VI3D <- sum(pmin(fhata$estimate, fhatb$estimate) * vol.cell)
  UDOI3D <- Overlap95 * sum(fhata$estimate * fhatb$estimate * vol.cell)
  
  # Return results as a list
  return(list(PercentOverlap95a = PercentOverlap95a, PercentOverlap95b = PercentOverlap95b, VI3D = VI3D, UDOI3D = UDOI3D, Vol95a = Vol95a, Vol95b = Vol95b, Vol50a = Vol50a, Vol50b = Vol50b))
}


