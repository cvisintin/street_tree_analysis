library(spatstat)
library(spatstat.Knet)
library(sf)
library(foreach)
library(future)
library(future.apply)
plan(multisession)

moy <- c("january", "february", "march", "april", "may", "june", "july",
         "august", "september", "october", "november", "december")

trees <- read_sf(dsn = "data/gis/", layer = "Averley_Tree_Layout-Baseline")
colnames(trees) <- tolower(colnames(trees))
for(i in moy){
  trees[ , i] <- as.factor(trees[ , i, drop = TRUE])
}
trees$id <- seq_len(nrow(trees))
trees <- st_transform(trees, crs = 28355) # e.g. GDA94 / MGA zone 55

roads <- read_sf(dsn = "data/gis/", layer = "Averley_Road_Centerlines")
roads <- st_union(roads)
roads <- st_transform(roads, crs = 28355) # e.g. GDA94 / MGA zone 55

pts <- st_crop(trees, xmin = 371820, ymin = 5786320, xmax = 371970, ymax = 5786470)
lns <- st_crop(roads, xmin = 371820, ymin = 5786320, xmax = 371970, ymax = 5786470)

pts <- trees
lns <- roads

m_id <- vector()
for(i in 1:12) {
  if(!all(is.na(st_drop_geometry(pts[, moy[i]])))) m_id <- c(m_id, i)
}
missing_id <- setdiff(1:12, m_id)

pts_on_lns <- st_cast(st_nearest_points(pts, lns),
                      "POINT")[seq(2, nrow(pts)*2, 2)]
pts_ppp <- as.ppp(pts_on_lns)

plot(lns, col = "grey")
points(pts$geometry, cex = 0.5, pch = 20)
points(pts_on_lns, cex = 0.5, pch = 20, col = "darkred")

pts_ppp <- as.ppp(pts)
marks(pts_ppp) <- st_drop_geometry(pts)[, moy]

lns_psp <- as.psp(lns)
L <- as.linnet(lns_psp)

linear_network_CSR <- rpoislpp(intensity(lpp(pts_ppp, L)), L)

expected_K <- envelope(linear_network_CSR, linearK, nsim = 5)
plot(expected_K, main = "Expected K Trend", legend = FALSE)

linear_network <- lpp(pts_ppp, L)
measured_k <- envelope(linear_network, linearK, nsim = 5)
plot(measured_k)

lpps <- foreach(i = 1:12) %dopar% {
  #uncomment below and ignore previous two lines for sequential execution
  #lpps <- foreach(i = m_id) %do% {
  if (i %in% missing_id) { NULL
  }else{
    ppp <- as.ppp(pts[which(pts[[moy[i]]] == 1), "geometry"])
    linear_network <- lpp(ppp, L)
    envelope(linear_network, linearK, nsim = 5)
  }
}
names(lpps) <- moy

par(mfrow = c(4, 3))
for(i in 1:12) {
  par(mar = c(3, 1.5, 3, 1.5))
  if (i %in% missing_id) {plot.new()
  }else{
    plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "")
  }
}
