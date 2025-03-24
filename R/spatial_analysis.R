library(spatstat)
library(stlnpp)
library(sf)
library(RColorBrewer)
library(reshape2)
library(foreach)
library(doMC)

# Read in spatial data for linear network (roads)
roads <- read_sf(dsn = "data/gis/", layer = "Averley_Road_Centerlines")

# Union all of the road segments
roads <- st_union(roads)

# Project geometry to a coordinate system that uses meters if in
# decimal-based system
roads <- st_transform(roads, crs = 28355) # e.g. GDA94 / MGA zone 55

#### Test effects of changing layouts based on three scenarios ####

# Read in script for function that analyses a proposed GIS layout
source("R/spatial_analysis_functions.R")

# Read in and re-project GIS point data (trees) for baseline scenario and run spatial analysis
baseline <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Baseline")
baseline <- st_transform(baseline, crs = 28355) # e.g. GDA94 / MGA zone 55
baseline_k <- run_spatial_analysis(point_pattern = baseline,
                                   linear_network = roads,
                                   scenario = "baseline",
                                   return_lpps = TRUE)

# Read in and re-project GIS point data (trees) for improved scenario and run spatial analysis
improved <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Improved")
improved <- st_transform(improved, crs = 28355) # e.g. GDA94 / MGA zone 55
improved_k <- run_spatial_analysis(point_pattern = improved,
                                   linear_network = roads,
                                   scenario = "improved",
                                   return_lpps = TRUE)

# Read in and re-project GIS point data (trees) for gold-standard scenario and run spatial analysis
gold_standard <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Gold_Standard")
gold_standard <- st_transform(gold_standard, crs = 28355) # e.g. GDA94 / MGA zone 55
gold_standard_k <- run_spatial_analysis(point_pattern = gold_standard,
                                        linear_network = roads,
                                        scenario = "gold_standard",
                                        return_lpps = TRUE)

# Generate an expected point pattern of complete spatial randomness (CSR)
roads_psp <- as.psp(roads)
L <- as.linnet(roads_psp)
linear_network_CSR <- rpoislpp(intensity(lpp(as.ppp(baseline), L)), L)
expected_K <- envelope(linear_network_CSR, linearK, nsim = 20)

# Plot to compare all scenarios
png(paste0("figs/K_all_scenarios.png"), width = 2700, height = 3300, res = 300)
plot_k(lpps_list = list(baseline_k, improved_k, gold_standard_k), lpps_csr = expected_K)
dev.off()
