# toy model script
# implementing model for:
# t = 20 years
# n x m = 25 x 25 pixels = 125 pixels = 3750 m2

# states
# 1 = Undisturbed
# 2 = Degraded
# 3 = Deforested
# 4 = Regrown


############################################################
#### Dependencies and directory setup ####
#
#
library(terra)
library(glue)

setwd("C:/Users/nburg/Documents/c119b/c119b/toymodel")


############################################################
#### Set parameters ####
#
#
tf <- 20 # total number of timesteps
n <- 25 # map height
m <- 25 # map width


#####################################################################
#### Intialize variables ####
#
#


r0 <- rast(nrows = n, ncols = m) #10x10 empty raster
values(r0) <- sample(1, ncell(r0), replace = TRUE) #fill raster w all undisturbed forest
r0[1] <- 3 #start with one deforested cell in the top left corner

levels(r0) <- data.frame( # this MUST be called AFTER setting all the raster values. raster metadata gets stripped when setting new values
  ID = 1:4,
  State = c("Undisturbed", "Degraded", "Deforested", "Regrown")
)
# i think these probabilities can be thought of as "natural" transitions between states without human/climate impact? might be wrong
P <- matrix( #row = current state, col = next state
  c( # U  G     F     R (to)
    0.99, 0.01, 0,    0,   # U (from)
    0,    0.95, 0.01, 0.04,   # G
    0,    0.01, 0.99, 0,   # F         # 0 prob of going from F -> R bc it needs to go back to G first? thoughts?
    0,    0.09, 0.01, 0.95    # R 
  ),
  nrow = 4,
  byrow = TRUE
)

w <- matrix(1, 3, 3) #create a 3x3 neighborhood weight matrix
w[2, 2] <- 0



##############################################################
#### Define funcs ####
#
#

neighbor_frac <- function(r, state) { #computes fraction of neighboring cells in given state
  focal(r == state, w, sum, na.rm = TRUE)  / 8 #terra::focal(), compute moving window values for each cell. window must have odd dims. 
} # w = window/weight matrix, sum = function. 

modify_P <- function(P, u_frac, g_frac, f_frac, r_frac) { #injects spatial component . why choice of 0.3? arbitrary?
  P2 <- P
  #P2[3, 3] <- P2[3, 3] + 0.3 * f_frac #3,3 is prob of going from F to F. increases depending on f_frac
  P2[1, 2] <- P2[1, 2] + 0.3 * g_frac #1,2 is prob of going from U to G. increases depending on g_frac
  P2[1,3]  <- P2[1, 3] + 0.5 * f_frac #1,3 is prob of going from U to F. increases depending on f_frac
  P2[3, 4] <- P2[3, 4] + 0.1 * u_frac #3,4 is prob of going from F to R. increases depending on u_frac
  P2[2, 4] <- P2[2, 4] + 0.1 * u_frac #2,1 is prob of going from G to R. increases depending on u_frac
  P2[3, ] <- P2[3, ] / sum(P2[3, ]) # rescale all transitions from F so they sum to 1.
  P2[2, ] <- P2[2, ] / sum(P2[2, ]) # rescale all transitions from G so they sum to 1.
  P2
}

step_forest <- function(r) {
  mf <- neighbor_frac(r, 3) #sum neighbor states that are deforested (state 3)
  mg <- neighbor_frac(r, 2) # ditto for G
  mu <- neighbor_frac(r, 1) #sum neighbor states that are undisturbed (state 1)
  mr <- 1 - mf - mg- mu
  r_new <- r 
  
  for (i in 1:ncell(r)) { #for each cell
    s <- values(r)[i] #get state
    Pi <- modify_P(P, values(mu)[i], values(mg)[i], values(mf)[i], values(mr)[i])[s, ]  #get fraction of F neighbors, modify P matrix, extract row s

    values(r_new)[i] <- sample(1:4, 1, prob = Pi) # sample new state based on new P
  }
  
  r_new
}





##################################################################
#### Setup graph labels ####
#
#

forest_cols <- c(
  "#1b7837",  # 1 = Undisturbed (dark green)
  "#90d127",  # 2 = Degraded (light green)
  "#d73027",  # 3 = Deforested (red)
  "#26fd55"   # 4 = Regrown (teal)
)



######################################################################
#### Main loop ####
#
#

r <- r0

cat("Plotting initial state ...")  # print progress
png("map_t_series/toy_map_t0.png", width=800, height=800)  # open device
dev.off()
png("map_t_series/toy_map_live.png", width=800, height=800)  # open device
terra::plot(
  r0,
  col = forest_cols,
  type = "classes",
  all_levels = TRUE, 
  main = "Timestep 0",
  axes = FALSE
)
dev.off()     

for (t in 1:tf) {
  cat("Simulating timestep", t, "...\n")
  r <- step_forest(r)
  
  # Ensure levels persist
  levels(r) <- levels(r0)
  
  # Save the timestamped file
  png(glue("map_t_series/toy_map_t{t}.png"), width=800, height=800)
  terra::plot(r, col = forest_cols, type = "classes", all_levels = TRUE, 
              main = glue("Timestep {t}"), axes = FALSE)
  dev.off()
  
  # Save the "live" (latest) file separately
  png("map_t_series/toy_map_live.png", width=800, height=800)
  terra::plot(r, col = forest_cols, type = "classes", all_levels = TRUE, 
              main = glue("Timestep {t}"), axes = FALSE)
  dev.off()
}













