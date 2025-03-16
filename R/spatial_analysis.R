library(spatstat)
library(stlnpp)
library(sf)
library(RColorBrewer)
library(reshape2)
library(foreach)
library(doMC)

# Define some general parameters
moy <- c("january", "february", "march", "april", "may", "june", "july",
         "august", "september", "october", "november", "december")
gr_palette <- grey.colors(length(moy), start = 0.9, end = 0.1)
col_palette <- brewer.pal(length(moy), "Paired")

# Read in spatial data (note, geometry type is points for trees and polygon for
# site boundary)
trees <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Baseline")
roads <- read_sf(dsn = "data/gis/", layer = "Averley_Road_Centerlines")
# boundary <- read_sf(dsn = "data/gis/", layer = "Averley_Boundary")

# Union all of the road segments
roads <- st_union(roads)

# Project geometry to a coordinate system that uses meters if in
# decimal-based system
trees <- st_transform(trees, crs = 28355) # e.g. GDA94 / MGA zone 55
roads <- st_transform(roads, crs = 28355) # e.g. GDA94 / MGA zone 55
# boundary <- st_transform(boundary, crs = 28355) # e.g. GDA94 / MGA zone 55

# Change case to all lower in spatial data attribute names
colnames(trees) <- tolower(colnames(trees))

# Extract coordinates from points
trees_coords <- st_coordinates(trees)

# Combine coordinates with tree attribute data
trees <- cbind(trees_coords, trees)

# Assign IDs to records in data
trees$id <- seq_len(nrow(trees))

# Record number of datapoints
n_data <- nrow(trees)

# Determine which months have flowering trees
m_id <- vector()
for(i in 1:12) {
  if(!all(is.na(st_drop_geometry(trees[, moy[i]])))) m_id <- c(m_id, i)
}

# Determine months with no flowering trees
missing_id <- setdiff(1:12, m_id)

#### Descriptive/Visual Analysis ####
# Examine the distribution flowering trees by month
png("figs/month_flowering_baseline.png", width = 900, height = 1100, res = 100)
plot(trees[, moy], max.plot = 12, pal = "black", pch = 20, cex = 0.5)
dev.off()

# Calculate total month that each tree is in flower - note that the geometry needs
# to be dropped before passing to the rowSums function
idx_cols <- which(colnames(trees) %in% moy)
trees$total_months <- rowSums(st_drop_geometry(trees[, moy]), na.rm = TRUE)

# Calculate the overall distribution of annual proportion flowering
png("figs/dist_flowering_baseline.png", width = 900, height = 900, res = 100)
par(mar = c(5, 5, 1.5, 1.5))
h <- hist(trees$total_months,
          main = "",
          xlab = "Total Months Flowering",
          probability = TRUE,
          xaxt = 'n',
          breaks = 12,
          xlim = c(0, 13))
axis(side = 1, at = seq(0.5, 12.5), labels = seq(0, 12))
lines(density(trees$total_months))
dev.off()

# Examine the spatial arrangement of annual proportion flowering
png("figs/total_months_flowering_baseline.png", width = 1100, height = 900, res = 100)
plot(trees["total_months"],
     pal = c("white", gr_palette),
     pch = 20,
     cex = 0.7,
     main = "",
     breaks = 0:13)
dev.off()

# Record how many trees are non-flowering
n_non_flower <- sum(trees$total_months == 0)

# Calculate the overall Shannon Diversity Index
sp_count <- aggregate(id ~ species, data = st_drop_geometry(trees), FUN = length)
sp_props <- sp_count$id / sum(sp_count$id)
sdi <- -sum(sp_props * log(sp_props))
paste0("Shannon Diversity Index: ", sdi)

# Calculate the overall Shannon Equitability Index
sei <- sdi / log(length(unique(sp_count$species)))
paste0("Shannon Diversity Index: ", sei)

#### Spatial Statistical Analysis ####
# Create planar segment pattern from road geometry & convert to linear network
roads_psp <- as.psp(roads)
L <- as.linnet(roads_psp)

# Create planar point pattern from tree locations
trees_ppp <- as.ppp(trees$geometry)
marks(trees_ppp) <- st_drop_geometry(trees)[, moy]

# Create a hypothetical point pattern assuming complete spatial randomness (CSR)
linear_network_CSR <- rpoislpp(intensity(lpp(as.ppp(trees$geometry), L)), L)

png("figs/lpp_CSR.png", width = 1100, height = 900, res = 100)
plot(linear_network_CSR, pch = 20, cex = 0.7, main = "")
dev.off()

# Evaluate the distribution using a modified K-function
expected_K <- envelope(linear_network_CSR, linearK, nsim = 20) # takes a while

png("figs/K_CSR.png", width = 900, height = 900, res = 100)
plot(expected_K, main = "Expected K Trend", legend = FALSE)
dev.off()

# Create a linear point pattern from existing data
linear_network <- lpp(trees_ppp, L)

png("figs/lpp_baseline.png", width = 1100, height = 900, res = 100)
plot(linear_network, legend = FALSE, pch = 20, cex = 0.7, main = "")
dev.off()

# Create point patterns on a linear network for each month (in parallel by default,
# please note that this may use up to 16GB of memory - use sequential processing
# if memory is a limiting factor)
registerDoMC(cores = 6)
lpps <- foreach(i = 1:12) %dopar% {
  #uncomment below and ignore previous two lines for sequential execution
  #lpps <- foreach(i = m_id) %do% {
  if (i %in% missing_id) { NULL
  }else{
    ppp <- as.ppp(trees[which(trees[[moy[i]]] == 1), "geometry"])
    linear_network <- lpp(ppp, L)
    envelope(linear_network, linearK, nsim = 20)
  }
}
names(lpps) <- moy

png("figs/K_baseline.png", width = 900, height = 1100, res = 100)
par(mfrow = c(4, 3))
for(i in 1:12) {
  par(mar = c(3, 1.5, 3, 1.5))
  if (i %in% missing_id) {plot.new()
  }else{
    plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "")
  }
}
dev.off()

# png("figs/K_baseline.png", width = 900, height = 1100, res = 100)
# par(mfrow = c(4, 3))
# for(i in 1:length(lpps)) {
#   if(i %in% c(1, 4, 7)) {
#     par(mar = c(1.5, 3, 3, 1.5))
#     plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "" , ylab = "K")
#   }
#   if(i %in% c(11, 12)){
#     par(mar = c(3, 1.5, 3, 1.5))
#     plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "")
#   }
#   if(i %in% c(10)){
#     par(mar = c(3, 3, 3, 1.5))
#     plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "K")
#   }
#   if(i %in% c(2, 3, 5, 6, 8, 9)){
#     par(mar = c(1.5, 1.5, 3, 1.5))
#     plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "" , ylab = "")
#   }
# }
# dev.off()


#### Proposed layouts ####
# Previous steps were for the baseline case, now test effects of changing layouts
# based on three scenarios

# Read in script for function that analyses a proposed GIS layout
source("R/spatial_analysis_functions.R")

# Read in and re-project GIS data for realistic scenario and run spatial analysis
baseline <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Baseline")
baseline <- st_transform(baseline, crs = 28355) # e.g. GDA94 / MGA zone 55
baseline_k <- run_spatial_analysis(point_pattern = baseline, roads = roads, scenario = "baseline", return_lpps = TRUE)
baseline_n_species <- length(unique(baseline$species))
baseline_pct_nonflowering <- sum(apply(st_drop_geometry(baseline)[, 3:14], 1, function(x) all(is.na(x)))) / nrow(baseline)

# Read in and re-project GIS data for improved scenario and run spatial analysis
improved <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Improved")
improved <- st_transform(improved, crs = 28355) # e.g. GDA94 / MGA zone 55
improved_k <- run_spatial_analysis(point_pattern = improved, roads = roads, scenario = "improved", return_lpps = TRUE)
improved_n_species <- length(unique(improved$species))
improved_pct_nonflowering <- sum(apply(st_drop_geometry(improved)[, 3:14], 1, function(x) all(is.na(x)))) / nrow(improved)

# Read in and re-project GIS data for gold-standard scenario and run spatial analysis
gold_standard <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Gold_Standard")
gold_standard <- st_transform(gold_standard, crs = 28355) # e.g. GDA94 / MGA zone 55
gold_standard_k <- run_spatial_analysis(point_pattern = gold_standard, roads = roads, scenario = "gold_standard", return_lpps = TRUE)
gold_standard_n_species <- length(unique(gold_standard$species))
gold_standard_pct_nonflowering <- sum(apply(st_drop_geometry(gold_standard)[, 3:14], 1, function(x) all(is.na(x)))) / nrow(gold_standard)

png(paste0("figs/K_all_scenarios.png"), width = 1800, height = 2200, res = 200)
plot_k(list(baseline_k, improved_k, gold_standard_k))
dev.off()



# #### Spatiotemporal Statistical Analysis ####
# # Reshape data to provide a datum for each month that each tree is flowering (e.g.
# # some number between X trees and X trees multiplied by 12 - based on varying 
# # number of flowering months per species)
# trees_long <- melt(st_drop_geometry(trees), id = c("X", "Y", "id", "species"))
# trees_long <- na.omit(trees_long[trees_long$value == 1, ])
# 
# # Record how many trees are flowering (note, unique values must be counted due to
# # the long format of the data)
# n_flower <- length(unique(trees_long$id))
# 
# # Verify that the data was reshaped correctly
# (n_non_flower + n_flower == n_data)
# 
# # Examine example species to verify flowering in expected months - e.g. "Angophora
# # costa" should only flower OCT-DEC
# trees_long[trees_long$species == "Angophora costa", ]
# 
# # Update the "value" column to code the respective month
# for (m in m_id) {
#   idx <- which(trees_long$variable == moy[m])
#   trees_long$value[idx] <- m
# }
# 
# colnames(trees_long)[ncol(trees_long)] <- "t"
# 
# # Create planar segment pattern from road geometry & convert to linear network
# roads_psp <- as.psp(roads)
# L <- as.linnet(roads_psp)
# 
# # Create planar point pattern from tree locations
# trees_ppp <- as.ppp(trees_long[ , c("X", "Y")], W = boundary)
# 
# linear_network_st <- stlpp(trees_ppp, L, T = trees_long$t)
# 
# mod <- STLK(linear_network_st)
# 
# save(mod, file = "output/STLK")
# 
# plot(mod)
# 
# # Extract coordinates of vertices from the boundary polygon
# boundary_coords <- st_coordinates(boundary)[ , 1:2]
# 
# # Create a point process for analysis and re-scale spatial units to kilometers -
# # use all month of year data as marks
# data_ppp <- as.ppp(st_drop_geometry(data[ , c("X", "Y")]), W = boundary)
# marks(data_ppp) <- st_drop_geometry(data[ , moy])
# data_ppp <- rescale(data_ppp, 1000, "km")
# 
# plot(data_ppp, s.region = boundary, pch = 20, use.marks = TRUE, cols = "black", maxsize = 0.001, legend = FALSE, main = "")
# 
# # Calculate density
# k <- density(data_ppp)
# plot(k, main=NULL, las=1)
# 
# # Test for clustering or dispersion based on point pair correlation
# g <- pcf(data_ppp)
# plot(g, main = NULL, las = 1)
# 
# 
# 
# # Create a point process for analysis and re-scale spatial units to kilometers -
# # duplicate points based on flowering months of year 
# data_long_ppp <- as.ppp(data_long[ , c("X", "Y", "value")], W = boundary)
# data_long_ppp <- rescale(data_long_ppp, 1000, "km")
# 
# plot(data_long_ppp, s.region = boundary, pch = 20, use.marks = TRUE, cols = "black", maxsize = 0.001, legend = FALSE, main = "")
# 
# # Test for clustering or dispersion based on point pair correlation
# g <- pcf(data_long_ppp)
# plot(g, main = NULL, las = 1)


##### Scenarios #####
# Baseline - Boulevard in non-permitted sections gets more flowering trees: Changing median trees on boulevard outside permit 1 and regular street trees on the boulevard east-west to east of creek - Yellow gum, On the creek - Black Iron Bark
# Improved - Entire Boulevard gets more flowering trees: Realistic but also across permit 1 - median and roadside
# Gold Standard - Mix of flowering trees throughout development: Replacing Acers (2), Japanese Zelkova, and Elms (2) with ??? - create map with non-flowering...