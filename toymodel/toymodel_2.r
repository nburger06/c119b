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
install.packages("tidyterra")
library(tidyterra)
library(ggplot2)


setwd()


############################################################
## Global script settings and constants ##
#
#
tf <- 60  # total number of timesteps
n <- 10 # map height. set to 10 to run fast while editing. 
m <- 10 # map width

h0 <- 2 # human land use impact scenario. 1,2, or 3. 1 is low impact, 2 is business as usual, 3 is high impact.
Kappa <- 1.5 # climate change impact scenario, defined as the temperature by 2100.

# initial natural transition states
alpha_0 <- 0
beta_0 <- 0
gamma_0 <- 0.01
psi_0 <- 0.04

## Calculate parameters based on scenarios.
kappa <- (Kappa / 74) # this is the climate change constant. calculated as:
# 1.5 degrees C warming between 2026 and 2100, / 74 years to estimate per-year
# linear change. this corresponds to the Intermediate (SSP2-4.5) scenario, total of 2.7 degrees warming
# since pre-industrial levels. 
alpha_c <- 0.05 # climate change deforestation impact constant
beta_c <- 0.1 # climate change degradation impact constant
gamma_c <- 0.005 # climate change secondary deforestation impact constant

if (h0 == 1) { #low impact rates, based on 1990-2000 reported rates in Ernst et al 2013, DRC
  ar_h <- 0.0015 # human deforestation rate
  br_h <- 0.0007 # human degradation rate
  gr_h <- 0.0015 # human secondary rate
  pr_h <- 0.0020 # human regrowth rate
} else if (h0 == 2) { #business-as-usual (BAU) rates based on Ernst et al 2013, DRC
  ar_h <- 0.0032 # human deforestation rate
  br_h <- 0.0016 # human degradation rate
  gr_h <- 0.0032 # human secondary rate
  pr_h <- 0.0010 # human regrowth rate
} else if (h0 == 3) { #high impact rates, based on doubling the BAU rates
  ar_h <- 0.0064 # human deforestation rate
  br_h <- 0.0032 # human degradation rate
  gr_h <- 0.0064 # human secondary rate
  pr_h <- 0.0010 # human regrowth rate
} else {
  stop("Invalid h0 value. Must be 1, 2, or 3.")
}
# gonna consolidate regrowth into one param for now



#####################################################################
#### Intialize variables ####
#
#



# initialize forest raster and set initial states
r0 <- rast(nrows = n, ncols = m) #10x10 empty raster
crs(r0) <- "local" # explicitly set non-geographic CRS to avoid bleeding between edges (otherwise terra assumes the edges wrap around and form a "globe")
values(r0) <- sample(1, ncell(r0), replace = TRUE) #fill raster w all undisturbed forest
r0[1] <- 3 #start with one deforested cell in the top left corner
r0[2] <- 3 #start with one deforested cell in the top left corner
r0[3] <- 3 #start with one deforested cell in the top left corner



levels(r0) <- data.frame( # this MUST be called AFTER setting all the raster values. raster metadata gets stripped when setting new values
  ID = 1:4,
  State = c("Undisturbed", "Degraded", "Deforested", "Regrown")
)

# thinking these probabilities can be thought of as "natural" transitions between states without human/climate impact? might be wrong
ug <- beta_0 
uf <- alpha_0 
ur <- 0
uu <- 1 - ug - uf - ur
gu <- 0
gf <- gamma_0
gr <- psi_0
gg <- 1 - gu - gf - gr
fu <- 0
fg <- psi_0
fr <- 0
ff <- 1 - fg - fu - fr
ru <- 0
rg <- beta_0
rf <- alpha_0
rr <- 1 - ru - rg - rf

P <- matrix( #row = current state, col = next state
  c( # U  G     F     R (to)
    uu,     ug,      uf,     ur,   # U (from)
    gu,     gg,      gf,     gr,   # G
    fu,     fg,      ff,     fr,   # F         # 0 prob of going from F -> R, change this later if we want F -> R possibility
    ru,     rg,      rf,     rr    # R 
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
neighbor_frac <- function(r, state) {
  # Count matching neighbors
  n_match <- focal(r == state, w, sum, na.rm = TRUE, pad = TRUE, fillvalue = NA)
  # Count total existing neighbors (not NA)
  n_total <- focal(!is.na(r), w, sum, na.rm = TRUE, pad = TRUE)
  
  return(n_match / n_total)
}

modify_P <- function(P, u_frac, g_frac, f_frac, r_frac, alpha, beta, gamma, psi) { #update p matrix based on spatial fracs and time-dependent params
  P2 <- P
  #P2[3, 3] <- P2[3, 3] + 0.3 * f_frac #3,3 is prob of going from F to F. increases depending on f_frac
  P2[1,3]  <- P2[1, 3] + alpha * (f_frac + g_frac) #1,3 is prob of going from U to F. increases depending on f_frac
  P2[1, 2] <- P2[1, 2] + beta * (g_frac + f_frac) #1,2 is prob of going from U to G. dependence on f = edge effects
  P2[2, 3] <- P2[2, 3] + gamma * f_frac #2,3 is prob of going from G to F. increases depending on f_frac
  P2[3, 2] <- P2[3, 2] + psi * (u_frac + r_frac + g_frac) #3,2 is prob of going from F to G. increases depending on u and r and g frac
  P2[2, 4] <- P2[2, 4] + psi * (u_frac + r_frac) #2,4 is prob of going from G to R. increases depending on u and r frac
  P2[3, ] <- P2[3, ] / sum(P2[3, ]) # rescale all transitions from F so they sum to 1
  P2[2, ] <- P2[2, ] / sum(P2[2, ]) # rescale all transitions from G so they sum to 1
  P2[1, ] <- P2[1, ] / sum(P2[1, ]) # rescale all transitions from U so they sum to 1
  P2
}

step_forest <- function(r, alpha, beta, gamma, psi) {
  mf <- neighbor_frac(r, 3) #sum neighbor states that are deforested (state 3)
  mg <- neighbor_frac(r, 2) # ditto for G
  mu <- neighbor_frac(r, 1) #sum neighbor states that are undisturbed (state 1)
  mr <- 1 - mf - mg- mu
  r_new <- r 
  
  for (i in 1:ncell(r)) { #for each cell
    s <- values(r)[i] #get state
    Pi <- modify_P(P, values(mu)[i], values(mg)[i], values(mf)[i], values(mr)[i], alpha, beta, gamma, psi)[s, ]  #get fraction of F neighbors, modify P matrix, extract row s

    values(r_new)[i] <- sample(1:4, 1, prob = Pi) # sample new state based on new P
  }
  
  r_new
}

land_use <- function(tf, intensity) {
  # human pop at time t * intensity
  impact = h * alpha_h * beta_h * gamma_h 
  return (impact)
}



##################################################################
#### Setup graph labels ####
#
#

forest_cols <- c(
  "#6A9C52",  # 1 = Undisturbed (dark green)
  "#D96900",  # 2 = Degraded (light green)
  "#AA2A0E",  # 3 = Deforested (red)
  "#5DA5D3"   # 4 = Regrown (teal)
)



######################################################################
#### Main loop ####
#
#

r <- r0

#store values for plotting
alpha_values <- c()
beta_values <- c()
gamma_values <- c()
psi_values <- c()


for (t in 0:tf) {
  cat("Simulating timestep", t, "...\n")

  alpha_h <- 1 - exp(-ar_h * t)
  beta_h <- 1 - exp(-br_h * t)
  gamma_h <- 1 - exp(-gr_h * t)
  psi_h <- 1 - exp(-pr_h * t)
  
  alpha <- alpha_h + alpha_c * kappa * t
  beta <- beta_h + beta_c * kappa * t
  gamma <- gamma_h + gamma_c * kappa * t
  psi <- psi_h

  alpha_values <- c(alpha_values, alpha_h)    
  beta_values <- c(beta_values, beta)
  gamma_values <- c(gamma_values, gamma)
  psi_values <- c(psi_values, psi)

  if (t == 0) {
    r <- r0
  } else {
    r <- step_forest(r, alpha, beta, gamma, psi)
  }
  levels(r) <- levels(r0)

  r_df <- as.data.frame(r, xy = TRUE)
  colnames(r_df) <- c("x", "y", "state")
  r_df$state <- factor(r_df$state, levels = c("Undisturbed", "Degraded", "Deforested", "Regrown"))
  all_states <- c("Undisturbed", "Degraded", "Deforested", "Regrown")

  p <- ggplot(r_df, aes(x = x, y = y, fill = state)) +
  geom_raster() +
  scale_fill_manual(
    values = setNames(forest_cols, all_states),
    breaks = all_states,   # <-- force all legend entries
    limits = all_states,   # <-- force all legend entries
    drop = FALSE
  )  + labs(title = glue("Timestep {t}"), fill = "State") +
  theme_void() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  legend.key.size = unit(1.5, "cm")
  )
  ggsave(glue("map_t_series/toy_map_t{t}.png"), plot = p, width = 8, height = 8)
  ggsave("toy_map_live.png", plot = p, width = 8, height = 8)
}
# later change width = 1800, height = 1200
t_vector <- 0:tf

# combine all params into one dataframe for easy plotting
param_df <- data.frame(
  t = rep(t_vector, 4),
  value = c(alpha_values, beta_values, gamma_values, psi_values),
  parameter = rep(c("Alpha (U→F)", "Beta (U→G)", "Gamma (G→F)", "Psi (F→G)"), each = length(t_vector))
)

# shared theme
forest_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


# Alpha (U->F)
ggsave("param_figs/alpha_values.png", width = 8, height = 6,
  ggplot(param_df[param_df$parameter == "Alpha (U→F)", ], aes(x = t, y = value)) +
    geom_line(color = "#AA2A0E") + geom_point(color = "#AA2A0E") +
    labs(title = "Undisturbed to Deforested (Alpha)", x = "Timestep", y = "Alpha Value") +
    forest_theme
)

# Beta (U->G)
ggsave("param_figs/beta_values.png", width = 8, height = 6,
  ggplot(param_df[param_df$parameter == "Beta (U→G)", ], aes(x = t, y = value)) +
    geom_line(color = "#D96900") + geom_point(color = "#D96900") +
    labs(title = "Undisturbed to Degraded (Beta)", x = "Timestep", y = "Beta Value") +
    forest_theme
)

# Gamma (G->F)
ggsave("param_figs/gamma_values.png", width = 8, height = 6,
  ggplot(param_df[param_df$parameter == "Gamma (G→F)", ], aes(x = t, y = value)) +
    geom_line(color = "#6A9C52") + geom_point(color = "#6A9C52") +
    labs(title = "Degraded to Deforested (Gamma)", x = "Timestep", y = "Gamma Value") +
    forest_theme
)

# Psi (F->G)
ggsave("param_figs/psi_values.png", width = 8, height = 6,
  ggplot(param_df[param_df$parameter == "Psi (F→G)", ], aes(x = t, y = value)) +
    geom_line(color = "#5DA5D3") + geom_point(color = "#5DA5D3") +
    labs(title = "Deforested to Regrowth (Psi)", x = "Timestep", y = "Psi Value") +
    forest_theme
)

# bonus: all params on one plot
ggsave("param_figs/all_params.png", width = 8, height = 6,
  ggplot(param_df, aes(x = t, y = value, color = parameter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c(
      "Alpha (U→F)" = "#AA2A0E",
      "Beta (U→G)"  = "#D96900",
      "Gamma (G→F)" = "#6A9C52",
      "Psi (F→G)"   = "#5DA5D3"
    )) +
    labs(title = "All Transition Parameters Over Time", x = "Timestep", y = "Value", color = "Parameter") +
    forest_theme
)

trans_df <- data.frame(
  t = rep(0:tf, 4),
  value = c(uf, ug, gf, fr),
  transition = rep(c("mean P(U → F)", "mean P(U → G)", "mean P(G → F)", "mean P(F → R)"), each = length(0:tf))
)
trans_df$transition <- factor(trans_df$transition, levels = c("mean P(U → F)", "mean P(U → G)", "mean P(G → F)", "mean P(F → R)"))

ggsave("param_figs/average_transition_probs.png", width = 8, height = 6,
  ggplot(trans_df, aes(x = t, y = value, color = transition)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c(
      "mean P(U → F)" = "#AA2A0E",
      "mean P(U → G)" = "#D96900",
      "mean P(G → F)" = "#6A9C52",
      "mean P(F → R)" = "#5DA5D3"
    )) +
    labs(title = "Landscape Avg Transition Probabilities", x = "Timestep", y = "Probability", color = "Transition") +
    forest_theme
)












