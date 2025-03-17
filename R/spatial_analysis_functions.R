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
  
  missing_id <- setdiff(1:12, m_id)
  
  png(paste0("figs/month_flowering_", scenario,".png"), width = 900, height = 1100, res = 100)
  plot(point_pattern[, moy], max.plot = 12, pal = "black", pch = 20, cex = 0.5)
  dev.off()
  
  idx_cols <- which(colnames(point_pattern) %in% moy)
  point_pattern$total_months <- rowSums(st_drop_geometry(point_pattern[, moy]), na.rm = TRUE)
  
  png(paste0("figs/dist_flowering_", scenario,".png"), width = 900, height = 900, res = 100)
  par(mar = c(5, 5, 1.5, 1.5))
  h <- hist(point_pattern$total_months,
            main = "",
            xlab = "Total Months Flowering",
            ylab = "Proportion of Trees",
            probability = TRUE,
            xaxt = 'n',
            breaks = 12,
            xlim = c(0, 13))
  axis(side = 1, at = seq(0.5, 12.5), labels = seq(0, 12))
  lines(density(trees$total_months))
  dev.off()
  
  png(paste0("figs/total_months_flowering_", scenario,".png"), width = 1100, height = 900, res = 100)
  plot(point_pattern["total_months"],
       pal = c("white", gr_palette),
       pch = 20,
       cex = 0.7,
       main = "",
       breaks = 0:13)
  dev.off()
  
  sp_count <- aggregate(id ~ species, data = st_drop_geometry(point_pattern), FUN = length)
  sp_props <- sp_count$id / sum(sp_count$id)
  sdi <- -sum(sp_props * log(sp_props))
  print(paste0("Shannon Diversity Index: ", sdi))
  
  # Calculate the overall Shannon Equitability Index (He = H / ln(S))
  sei <- sdi / log(length(unique(sp_count$species)))
  print(paste0("Shannon Equitability Index: ", sei))
  
  roads_psp <- as.psp(st_cast(roads, "MULTILINESTRING"))
  L <- as.linnet(roads_psp)
  point_pattern_ppp <- as.ppp(point_pattern$geometry)
  marks(point_pattern_ppp) <- st_drop_geometry(point_pattern)[, moy]
  linear_network <- lpp(point_pattern_ppp, L)
  
  png(paste0("figs/lpp_", scenario,".png"), width = 1100, height = 900, res = 100)
  plot(linear_network, legend = FALSE, pch = 20, lwd = 0.5, size = 1.25, cols = "palegreen3", cex = 0.7, main = "")
  dev.off()
  
  registerDoMC(cores = 6)
  lpps <- foreach(i = 1:12) %dopar% {
    #uncomment below and ignore previous two lines for sequential execution
    #lpps <- foreach(i = m_id) %do% {
    if (i %in% missing_id) { NULL
    }else{
      ppp <- as.ppp(point_pattern[which(point_pattern[[moy[i]]] == 1), "geometry"])
      linear_network <- lpp(ppp, L)
      envelope(linear_network, linearK, nsim = 20)
    }
  }
  names(lpps) <- moy
  
  png(paste0("figs/K_", scenario,".png"), width = 900, height = 1100, res = 100)
  par(mfrow = c(4, 3))
  for(i in 1:12) {
    par(mar = c(3, 1.5, 3, 1.5))
    if (i %in% missing_id) {plot.new()
    }else{
      plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "")
    }
  }
  dev.off()
  
  if (return_lpps == TRUE) return(lpps)
  
}

plot_k <- function(lpps_list) {
  par(mfrow = c(4, 3))
  ymax <- sapply(1:length(lpps_list), function(x) sapply(1:length(lpps_list[[1]]), function(y) max(lpps_list[[x]][[y]]$obs)))
  for(i in 1:12) {
    par(mar = c(0.5, 0.5, 2.5, 0.5))
    month <- paste0(toupper(substr(moy[i], 1, 1)), substr(moy[i], 2, nchar(moy[i])))
    plot(lpps_list[[1]][[i]]$r, lpps_list[[1]][[i]]$theo,
         col = "white", type = 'l', main = month, xlab = "",
         ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, max(ymax[i, ])))
    polygon(x = c(lpps_list[[1]][[i]]$r, rev(lpps_list[[1]][[i]]$r)),
            y = c(lpps_list[[1]][[i]]$lo, 
                  rev(lpps_list[[1]][[i]]$hi)),
            col =  "lightgray", border = NA)
    lines(lpps_list[[1]][[i]]$r, lpps_list[[1]][[i]]$theo,
         col = "darkgray", lty = 3)
    lines(lpps_list[[1]][[i]]$r, lpps_list[[1]][[i]]$obs, lty = 1)
    lines(lpps_list[[1]][[i]]$r, lpps_list[[2]][[i]]$obs, lty = 2)
    lines(lpps_list[[1]][[i]]$r, lpps_list[[3]][[i]]$obs, lty = 6)
  }
}
