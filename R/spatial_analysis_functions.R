run_spatial_analysis <- function(point_pattern, roads, scenario, return_lpps = FALSE) {
  
  moy <- c("january", "february", "march", "april", "may", "june", "july",
           "august", "september", "october", "november", "december")
  gr_palette <- grey.colors(length(moy), start = 0.9, end = 0.1)
  col_palette <- brewer.pal(length(moy), "Paired")
  
  colnames(point_pattern) <- tolower(colnames(point_pattern))
  point_pattern_coords <- st_coordinates(point_pattern)
  point_pattern <- cbind(point_pattern_coords, point_pattern)
  point_pattern$id <- seq_len(nrow(point_pattern))
  n_data <- nrow(point_pattern)
  
  m_id <- vector()
  for(i in 1:12) {
    if(!all(is.na(st_drop_geometry(point_pattern[, moy[i]])))) m_id <- c(m_id, i)
  }
  
  png(paste0("figs/month_flowering_", scenario,".png"), width = 900, height = 1100, res = 100)
  plot(point_pattern[, moy], max.plot = 12, pal = "black", pch = 20, cex = 0.5)
  dev.off()
  
  idx_cols <- which(colnames(point_pattern) %in% moy)
  point_pattern$total_months <- rowSums(st_drop_geometry(point_pattern[, moy]), na.rm = TRUE)
  
  png(paste0("figs/dist_flowering_", scenario,".png"), width = 900, height = 900, res = 100)
  h <- hist(point_pattern$total_months, main = "", xlab = "Total Months Flowering", probability = TRUE)
  lines(density(point_pattern$total_months))
  dev.off()
  
  png(paste0("figs/total_months_flowering_", scenario,".png"), width = 1100, height = 900, res = 100)
  plot(point_pattern["total_months"],
       pal = c("white", gr_palette),
       pch = 20,
       cex = 0.7,
       main = "",
       breaks = 0:13)
  dev.off()
  
  roads_psp <- as.psp(roads)
  L <- as.linnet(roads_psp)
  point_pattern_ppp <- as.ppp(point_pattern$geometry)
  marks(point_pattern_ppp) <- st_drop_geometry(point_pattern)[, moy]
  linear_network <- lpp(point_pattern_ppp, L)
  
  png(paste0("figs/lpp_", scenario,".png"), width = 1100, height = 900, res = 100)
  plot(linear_network, legend = FALSE, pch = 20, cex = 0.7, main = "")
  dev.off()
  
  registerDoMC(cores = 6)
  lpps <- foreach(i = m_id) %dopar% {
    #uncomment below and ignore previous two lines for sequential execution
    #lpps <- foreach(i = m_id) %do% {
    ppp <- as.ppp(point_pattern[which(point_pattern[[moy[i]]] == 1), "geometry"])
    linear_network <- lpp(ppp, L)
    envelope(linear_network, linearK, nsim = 20)
  }
  names(lpps) <- moy[m_id]
  
  png(paste0("figs/K_", scenario,".png"), width = 900, height = 1100, res = 100)
  par(mfrow = c(4, 3))
  for(i in 1:length(lpps)) {
    par(mar = c(3, 1.5, 3, 1.5))
    plot(lpps[[i]], main = moy[m_id][i], legend = FALSE, xlab = "T" , ylab = "")
  }
  dev.off()
  
  if (return_lpps == TRUE) return(lpps)
  
}
