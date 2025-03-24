run_spatial_analysis <- function(point_pattern, linear_network, scenario, return_lpps = FALSE) {
  
  moy <- c("january", "february", "march", "april", "may", "june", "july",
           "august", "september", "october", "november", "december")
  gr_palette <- grey.colors(length(moy), start = 0.9, end = 0.1)
  col_palette <- brewer.pal(length(moy), "Paired")
  
  colnames(point_pattern) <- tolower(colnames(point_pattern))
  point_pattern$id <- seq_len(nrow(point_pattern))
  
  m_id <- vector()
  for(i in 1:12) {
    if(!all(is.na(st_drop_geometry(point_pattern[, moy[i]])))) m_id <- c(m_id, i)
  }
  
  missing_id <- setdiff(1:12, m_id)
  
  png(paste0("figs/month_flowering_", scenario,".png"), width = 2700, height = 3300, res = 300)
  plot(point_pattern[, moy], max.plot = 12, pal = "black", pch = 20, cex = 0.4)
  dev.off()
  
  idx_cols <- which(colnames(point_pattern) %in% moy)
  point_pattern$total_months <- rowSums(st_drop_geometry(point_pattern[, moy]), na.rm = TRUE)

  png(paste0("figs/total_months_flowering_", scenario,".png"), width = 1100, height = 900, res = 100)
  plot(point_pattern["total_months"],
       pal = c("white", gr_palette),
       pch = 20,
       cex = 0.7,
       main = "",
       breaks = 0:13)
  dev.off()

  linear_network_psp <- as.psp(linear_network)
  L <- as.linnet(linear_network_psp)

  registerDoMC(cores = 6)
  lpps <- foreach(i = 1:12) %dopar% {
    # uncomment below and ignore (comment) previous two lines for sequential execution
    # lpps <- foreach(i = m_id) %do% {
    if (i %in% missing_id) { NULL
    }else{
      ppp <- as.ppp(point_pattern[which(point_pattern[[moy[i]]] == 1), "geometry"])
      linear_network <- lpp(ppp, L)
      envelope(linear_network, linearK, nsim = 20)
    }
  }
  names(lpps) <- moy
  
  png(paste0("figs/K_", scenario,".png"), width = 2700, height = 3300, res = 300)
  par(mfrow = c(4, 3))
  for(i in 1:12) {
    par(mar = c(3, 1.5, 3, 1.5))
    if (i %in% missing_id) { plot.new()
    }else{
      plot(lpps[[i]], main = moy[i], legend = FALSE, xlab = "T" , ylab = "")
    }
  }
  dev.off()
  
  if (return_lpps == TRUE) return(lpps)
  
}


plot_k <- function(lpps_list, lpps_csr) {
  moy <- c("january", "february", "march", "april", "may", "june", "july",
           "august", "september", "october", "november", "december")
  ymax <- ceiling(
    max(
      sapply(1:length(lpps_list),
             function(x) sapply(1:length(lpps_list[[1]]),
                                function(y) max(lpps_list[[x]][[y]]$obs)))))
  xmax <- ceiling(
    max(
      sapply(1:length(lpps_list),
             function(x) sapply(1:length(lpps_list[[1]]),
                                function(y) max(lpps_list[[x]][[y]]$r)))))
  
  par(mfrow = c(4, 3))
  
  for(i in 1:12) {
    par(mar = c(0.5, 0.5, 2.5, 0.5))
    month <- paste0(toupper(substr(moy[i], 1, 1)), substr(moy[i], 2, nchar(moy[i])))
    plot(NULL,
         col = "white", type = 'l', main = month, xlab = "",
         ylab = "", xaxt = "n", yaxt = "n", xlim = c(0, xmax),
         ylim = c(0, ymax))
    # lines(lpps_list[[1]][[i]]$r, lpps_list[[1]][[i]]$theo,
    #      col = "white", type = 'l', main = month, xlab = "",
    #      ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, max(ymax[i, ])))
    polygon(x = c(lpps_csr$r, rev(lpps_csr$r)),
            y = c(lpps_csr$lo, rev(lpps_csr$hi)),
            col =  "lightgray", border = NA)
    lines(lpps_csr$r, lpps_csr$theo, col = "darkgray", lty = 3)
    lines(lpps_list[[1]][[i]]$r, lpps_list[[1]][[i]]$obs, lty = 1)
    lines(lpps_list[[2]][[i]]$r, lpps_list[[2]][[i]]$obs, lty = 2)
    lines(lpps_list[[3]][[i]]$r, lpps_list[[3]][[i]]$obs, lty = 6)
  }
}
