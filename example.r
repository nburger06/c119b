#example.r

library(terra)

# states
# 1 = Undisturbed
# 2 = Degraded
# 3 = Deforested
# 4 = Regrown

r0 <- rast(nrows = 16, ncols = 16) #10x10 empty raster
values(r0) <- sample(1, ncell(r0), replace = TRUE) #fill raster w random integers representing forest states

P <- matrix( #row = current state, col = next state
  c( # U G F R (to)
    0.8, 0.1, 0.1, 0,   # U (from)
    0, 0.9, 0.25, 0.05,   # G
    0, 0, 0.95, 0.05,   # F
    0, 0.1, 0.1, 0.8     # R
  ),
  nrow = 4,
  byrow = TRUE
)

w <- matrix(1, 3, 3) #create a 3x3 neighborhood weight matrix
w[2, 2] <- 0

neighbor_frac <- function(r, state) {
  focal(r == state, w, sum, na.rm = TRUE) / 8 #terra::focal(), compute moving window values for each cell. window must have odd dims. 
} # w = window/weight matrix, sum = function. 

modify_P <- function(P, mature_frac) {
  P2 <- P
  P2[3, 3] <- P2[3, 3] + 0.3 * mature_frac
  P2[3, ] <- P2[3, ] / sum(P2[3, ])
  P2
}

step_forest <- function(r) {
  mf <- neighbor_frac(r, 3)
  r_new <- r
  
  for (i in 1:ncell(r)) {
    s <- values(r)[i]
    Pi <- modify_P(P, values(mf)[i])[s, ]
    values(r_new)[i] <- sample(1:4, 1, prob = Pi)
  }
  
  r_new
}

r <- r0
levels(r) <- data.frame(
  ID = 1:4,
  State = c("Undisturbed", "Degraded", "Deforested", "Regrown")
)
forest_cols <- c(
  "#1b7837",  # 1 = Undisturbed (dark green)
  "#90d127",  # 2 = Degraded (light green)
  "#d73027",  # 3 = Deforested (red)
  "#26fd55"   # 4 = Regrown (teal)
)

png("forest_map.png", width=800, height=800)  # open device
terra::plot(
    r,
    col = forest_cols,
    axes = FALSE,
    main = paste("Timestep", t)
  )  
dev.off()     
for (t in 1:10) {
  cat("Simulating timestep", t, "...\n")  # print progress
  r <- step_forest(r)
  png("forest_map.png", width=800, height=800)  # open device
  terra::plot(
    r,
    col = forest_cols,
    axes = FALSE,
    main = paste("Timestep", t)
  )  
  dev.off()     
}

# save as PNG
     

















