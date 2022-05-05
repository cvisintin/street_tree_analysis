library(spatstat)
library(stlnpp)
library(sf)
library(RColorBrewer)
library(reshape2)

# Define some general parameters
moy <- c("january", "february", "march", "april", "may", "june", "july",
         "august", "september", "october", "november", "december")
m_id <- 1:12
gr_palette <- grey.colors(length(moy), start = 0.8, end = 0.1)
col_palette <- brewer.pal(length(moy), "Paired")

# Read in spatial data (note, geometry type is points for trees and polygon for
# site boundary)
trees <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout")
boundary <- read_sf(dsn = "data/gis/", layer = "Averley_Boundary")
roads <- read_sf(dsn = "data/gis/", layer = "Averley_Road_Centerlines")

# Union all of the road segments
roads <- st_union(roads)

# Optional step: project geometry to a coordinate system that uses meters if in
# decimal-based system
trees <- st_transform(trees, crs = 8807) # e.g. WGS 84 / Pseudo-Mercator
boundary <- st_transform(boundary, crs = 8807) # e.g. WGS 84 / Pseudo-Mercator
roads <- st_transform(roads, crs = 8807) # e.g. WGS 84 / Pseudo-Mercator

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

# Assign colour palettes to points


#### Descriptive/Visual Analysis ####
# Examine the distribution flowering trees by month
plot(trees[, moy], max.plot = 12, pal = "black", pch = 20, cex = 0.5)

# Calculate proportion of year (total months / 12) that each tree is in flower -
# note that the geometry needs to be dropped before passing to the rowSums function
idx_cols <- which(colnames(trees) %in% moy)
trees$prop_year <- rowSums(st_drop_geometry(trees[, moy]), na.rm = TRUE) / 12

# What is the overall distribution of annual proportion flowering
h <- hist(trees$prop_year, main = "", xlab = "Annual Proportion Flowering", probability = TRUE)
lines(density(trees$prop_year, main = "", xlab = "Annual Proportion Flowering"))

# Examine the spatial arrangement of annual proportion flowering
plot(trees["prop_year"], pal = gr_palette[1:length(unique(trees$prop_year)) + 1], pch = 20, cex = 0.5, main = "")

# Record how many trees are non-flowering
n_non_flower <- sum(trees$prop_year == 0)

#### Spatial Statistical Analysis ####
# Create planar segment pattern from road geometry & convert to linear network
roads_psp <- as.psp(roads)
L <- as.linnet(roads_psp)

# Create planar point pattern from tree locations
trees_ppp <- as.ppp(trees$geometry)
marks(trees_ppp) <- st_drop_geometry(trees)[, moy]

# Create point patterns on a linear network for each month


for(i in m_id) {
  ppp <- as.ppp(trees[moy[i] == 1, "geometry"])
  linear_network <- lpp(ppp, L)
  
}

# summary(linear_network)
# plot(linear_network, legend = FALSE, main = "", pch = 20, maxsize = 1.0, col = "lightgray")


# # Move all tree locations to be coincident with the road centerlines - this is a
# # requirement of the linear pattern analysis (note, this may take a while depending
# # on the complexity of the road network and number of trees)
# n <- nrow(trees)
# trees_centerlines <- do.call(c,
#                              lapply(seq(n), function(i) {
#                                nrst = st_nearest_points(st_geometry(trees)[i], roads)
#                                nrst_len = st_length(nrst)
#                                nrst_mn = which.min(nrst_len)
#                                if (as.vector(nrst_len[nrst_mn]) > 25) return(st_geometry(trees)[i])
#                                return(st_cast(nrst[nrst_mn], "POINT")[2])
#                              }))


#### Spatiotemporal Statistical Analysis ####
# Reshape data to provide a datum for each month that each tree is flowering (e.g.
# some number between X trees and X trees multiplied by 12 - based on varying 
# number of flowering months per species)
trees_long <- melt(st_drop_geometry(trees), id = c("X", "Y", "id", "species"))
trees_long <- na.omit(trees_long[trees_long$value == 1, ])

# Record how many trees are flowering (note, unique values must be counted due to
# the long format of the data)
n_flower <- length(unique(trees_long$id))

# Verify that the data was reshaped correctly
(n_non_flower + n_flower == n_data)

# Examine example species to verify flowering in expected months - e.g. "Angophora
# costa" should only flower OCT-DEC
trees_long[trees_long$species == "Angophora costa", ]

# Update the "value" column to code the respective month
for (m in m_id) {
  idx <- which(trees_long$variable == moy[m])
  trees_long$value[idx] <- m
}

colnames(trees_long)[ncol(trees_long)] <- "t"

# Create planar segment pattern from road geometry & convert to linear network
roads_psp <- as.psp(roads)
L <- as.linnet(roads_psp)

# Create planar point pattern from tree locations
trees_ppp <- as.ppp(trees_long[ , c("X", "Y")], W = boundary)

linear_network_st <- stlpp(trees_ppp, L, T = trees_long$t)

STLK(linear_network_st)






# Extract coordinates of vertices from the boundary polygon
boundary_coords <- st_coordinates(boundary)[ , 1:2]

# Create a point process for analysis and re-scale spatial units to kilometers -
# use all month of year data as marks
data_ppp <- as.ppp(st_drop_geometry(data[ , c("X", "Y")]), W = boundary)
marks(data_ppp) <- st_drop_geometry(data[ , moy])
data_ppp <- rescale(data_ppp, 1000, "km")

plot(data_ppp, s.region = boundary, pch = 20, use.marks = TRUE, cols = "black", maxsize = 0.001, legend = FALSE, main = "")

# Calculate density
k <- density(data_ppp)
plot(k, main=NULL, las=1)

# Test for clustering or dispersion based on point pair correlation
g <- pcf(data_ppp)
plot(g, main = NULL, las = 1)



# Create a point process for analysis and re-scale spatial units to kilometers -
# duplicate points based on flowering months of year 
data_long_ppp <- as.ppp(data_long[ , c("X", "Y", "value")], W = boundary)
data_long_ppp <- rescale(data_long_ppp, 1000, "km")

plot(data_long_ppp, s.region = boundary, pch = 20, use.marks = TRUE, cols = "black", maxsize = 0.001, legend = FALSE, main = "")

# Test for clustering or dispersion based on point pair correlation
g <- pcf(data_long_ppp)
plot(g, main = NULL, las = 1)
