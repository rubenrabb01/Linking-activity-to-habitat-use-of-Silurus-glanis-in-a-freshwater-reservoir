################################################################################
######## Speed up computation time by running iterations in parallel ###########
########       Function Script                                      ###########
################################################################################
#CODE FROM:https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-starting-values.pdf


# Number of tries with different starting values
#niter <- 25

# Create list of starting values
#allPar0 <- lapply(as.list(1:niter), function(x) {
  # Step length mean
#  stepMean0 <- runif(2,
#                     min = c(5,25),
#                     max = c(25,50))
  # Step length standard deviation
#  stepSD0 <- runif(2,
 #                  min = c(5,25),
  #                 max = c(25,50))
  #zero mass step
#  zeromass0 <-  runif(2,
 #                     min = c(0.001,0.01),
  #                    max = c(0.01,0.05))
  # Turning angle mean
#  angleMean0 <- c(0, pi)
  # Turning angle concentration
 # angleCon0 <- runif(2,
  #                   min = c(0.5, 2),
   #                  max = c(1, 5))
  # Return vectors of starting values
#  stepPar0 <- c(stepMean0, stepSD0, zeromass0)
 # anglePar0 <- c(angleMean0, angleCon0)
  #return(list(step = stepPar0, angle = anglePar0))
#})


###########################################################################
##########################  3 state function     ########################## 
###########################################################################

# Number of tries with different starting values
niter <- 5

# Create list of starting values for a 3-state model
allPar0_3s <- lapply(as.list(1:niter), function(x) {
  # Step length mean for 3 states
  stepMean0 <- runif(3,
                     min = c(3, 25, 55),
                     max = c(50, 75, 100))
  # Step length standard deviation for 3 states
  stepSD0 <- runif(3,
                   min = c(3, 25, 55),
                   max = c(50, 75, 100))
 
  # Turning angle mean for 3 states
  angleMean0 <- c(0, pi, pi/2)  # Mean turning angles for the 3 states
  
  # Turning angle concentration for 3 states
  angleCon0 <- runif(3,
                     min = c(1, 2, 4),
                     max = c(3, 5, 6))
  
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)  # Concatenate all step parameters
  anglePar0 <- c(angleMean0, angleCon0)         # Concatenate all angle parameters
  
  return(list(step = stepPar0, angle = anglePar0))  # Return both step and angle parameters as a list
})

