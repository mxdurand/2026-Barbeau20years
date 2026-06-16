library(data.table)

df <- fread("data/B20Y_BM_all.csv")

# To long format
SWC_cols <- names(df)[grep("SWC", names(df))]
df <- df[rowSums(!is.na(df[, ..SWC_cols])) > 0] # Removing rows for which all SWC is NA (2005-2011 and more)
df <- df[, c("TIMESTAMP", SWC_cols), with = FALSE]
df <- melt(df, id.vars = "TIMESTAMP", variable.name = "sensor", value.name = "swc")
rm(SWC_cols)

# Date formatting
df[, y := substring(TIMESTAMP, 1 , 4)]
df[, m := substring(TIMESTAMP, 6 , 7)]
df[, d := substring(TIMESTAMP, 9 , 10)]
df[, h := substring(TIMESTAMP, 12 , 16)]
date_lookup <- unique(df[, .(y, m, d)])
date_lookup[, doy := yday(as.Date(paste(y, m, d, sep = "-")))]
df <- merge(df, date_lookup, by = c("y", "m", "d"))
rm(date_lookup)

# Sensor information in column name
df[, sensor := as.character(sensor)]
df[, h_pos  := as.integer(sub("SWC_([0-9]+)_([0-9]+)_1", "\\1", sensor))]
df[, v_dep  := as.integer(sub("SWC_([0-9]+)_([0-9]+)_1", "\\2", sensor))]

# Depth information
d0 <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45)
depth_map <- data.table(v_dep = 1:15, depth = d0)
df <- merge(df, depth_map, by = "v_dep")
rm(depth_map, d0)


# Standardize by position
pos_means <- df[, .(pos_mean = mean(swc, na.rm = TRUE)), by = .(h_pos, depth)]
df <- merge(df, pos_means, by = c("h_pos", "depth"))
df[, swc_anom := swc - pos_mean]
rm(pos_means)

# Remove top layer (often missing)
df <- df[depth != 0.05]


# Getting mean of positions
dfST <- df[, .(
  swc_mean  = mean(swc_anom, na.rm = TRUE), # position-corrected mean
  swc_abs   = mean(swc, na.rm = TRUE),      # raw absolute value if needed
  n_pos     = sum(!is.na(swc))              # now this is number of positions (max 6)
), by = .(y, m, d, doy, h, depth)]


# Add layer definition
dfST[depth <= 0.25, layer := "shallow"]
dfST[depth > 0.25 & depth <= 0.75, layer := "mid"]
dfST[depth > 0.75, layer := "deep"]
dfSTf <- copy(dfST)
dfSTf[, layer := "full"]
dfST <- rbind(dfST, dfSTf)
rm(dfSTf)

# Layer means with coverage tracking
dfSL <- dfST[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),
  swc_abs   = mean(swc_abs,  na.rm = TRUE),
  n_depths  = sum(!is.na(swc_mean)),  # how many depths contributed
  n_pos_mean = mean(n_pos, na.rm = TRUE)   # average positions per timestep
), by = .(y, m, d, doy, h, layer)]


# Daily averages
dfSJ <- dfST[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),  
  swc_abs   = mean(swc_abs, na.rm = TRUE),        
  n_hours    = sum(!is.na(swc_mean))        # valid half-hours in the day
), by = .(y, m, d, doy, layer)]
dfSJ[, y := as.numeric(y)]

min_hours <- 24   # at least 24 out of 48 half-hours otherwise day removed
dfSJ[n_hours < min_hours, swc_mean := NA]
dfSJ[n_hours < min_hours, swc_abs  := NA]
rm(min_hours)

# Calculate REW
rew_params <- dfST[, .(
  FC = quantile(swc_abs, 0.95, na.rm = TRUE),
  WP = quantile(swc_abs, 0.05, na.rm = TRUE)
), by = layer]

dfSJ <- merge(dfSJ, rew_params, by = "layer")
dfSJ[, REW := (swc_abs - WP) / (FC - WP)]
dfSJ[REW < 0, REW := 0] ; dfSJ[REW > 1, REW := 1]
rm(rew_params)

# Subset for the growing season
dfGS <- fread("data/df_inGS.csv")
dfGS[, m := ifelse(m < 10, paste0("0", as.character(m)), as.character(m))]
dfGS[, d := ifelse(d < 10, paste0("0", as.character(d)), as.character(d))]
dfSJ <- merge(dfSJ, dfGS, by = c("y", "m", "d", "doy"), all.x = T)
rm(dfGS)

# Yearly average
dfSY <- dfSJ[in_gs == TRUE & !is.na(swc_mean), .(
  swc_mean = mean(swc_mean, na.rm = TRUE),
  swc_abs  = mean(swc_abs,  na.rm = TRUE),
  rew_mean  = mean(REW,  na.rm = TRUE),
  n_days   = sum(!is.na(swc_mean))
), by = .(y, layer)]


# Computing number of stress days & drought intensity
stress_threshold <- 0.3

stress_days <- dfSJ[in_gs == TRUE, .(
  n_valid    = sum(!is.na(REW)),
  pct_stress = round(mean(REW < stress_threshold, na.rm = TRUE) * 100, 1)
), by = .(y, layer)]

dfSJ[, deficit := pmax(stress_threshold - REW, 0)]
drought_intensity <- dfSJ[in_gs == TRUE, .(
  n_valid      = sum(!is.na(deficit)),
  mean_deficit = mean(deficit, na.rm = TRUE)
), by = .(y, layer)]


# Plots for yearly trends -------------------------------------------------------------------------

# mSS <- lm(swc_mean ~ y, dfSY[layer == "shallow"]) ; summary(mSS)
# mSM <- lm(swc_mean ~ y, dfSY[layer == "mid"]) ; summary(mSM)
# mSD <- lm(swc_mean ~ y, dfSY[layer == "deep"]) ; summary(mSD)
# mSF <- lm(swc_mean ~ y, dfSY[layer == "full"]) ; summary(mSF)
# 
# mSSa <- lm(swc_abs ~ y, dfSY[layer == "shallow"]) ; summary(mSSa)
# mSMa <- lm(swc_abs ~ y, dfSY[layer == "mid"]) ; summary(mSMa)
# mSDa <- lm(swc_abs ~ y, dfSY[layer == "deep"]) ; summary(mSDa)
# mSFa <- lm(swc_abs ~ y, dfSY[layer == "full"]) ; summary(mSFa)
# 
# mRS <- lm(rew_mean ~ y, dfSY[layer == "shallow"]) ; summary(mRS)
# mRM <- lm(rew_mean ~ y, dfSY[layer == "mid"]) ; summary(mRM)
# mRD <- lm(rew_mean ~ y, dfSY[layer == "deep"]) ; summary(mRD)
# mRF <- lm(rew_mean ~ y, dfSY[layer == "full"]) ; summary(mRF)
# 
# mPS <- lm(pct_stress ~ y, stress_days[layer == "shallow"]) ; summary(mPS)
# mPM <- lm(pct_stress ~ y, stress_days[layer == "mid"]) ; summary(mPM)
# mPD <- lm(pct_stress ~ y, stress_days[layer == "deep"]) ; summary(mPD)
# mPF <- lm(pct_stress ~ y, stress_days[layer == "full"]) ; summary(mPF)
# 
# mDS <- lm(mean_deficit ~ y, drought_intensity[layer == "shallow"]) ; summary(mDS)
# mDM <- lm(mean_deficit ~ y, drought_intensity[layer == "mid"]) ; summary(mDM)
# mDD <- lm(mean_deficit ~ y, drought_intensity[layer == "deep"]) ; summary(mDD)
# mDF <- lm(mean_deficit ~ y, drought_intensity[layer == "full"]) ; summary(mDF)
# 
# png("figs/B20Y_annualSWC.png", height = 10000, width = 12000, res = 900)
# par(mfrow = c(2,2), bty = "L", mar = c(5,4,1,1), oma = c(0,1,0,0))
# 
# plot(-500, xlim = c(2011,2024), ylim = c(0,0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(swc_abs ~ y, dfSY[layer == "full"], col = "#3B1F2B", type = "l", lwd = 2)
# points(swc_abs ~ y, dfSY[layer == "deep"], col = "#E15554", type = "l", lwd = 2)
# points(swc_abs ~ y, dfSY[layer == "mid"], col = "#F6AE2D", type = "l", lwd = 2)
# points(swc_abs ~ y, dfSY[layer == "shallow"], col = "#2E86AB", type = "l", lwd = 2)
# points(swc_abs ~ y, dfSY[layer == "full"], pch = 21, bg = "#3B1F2B", type = "p", cex = 1.4)
# points(swc_abs ~ y, dfSY[layer == "deep"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
# points(swc_abs ~ y, dfSY[layer == "mid"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
# points(swc_abs ~ y, dfSY[layer == "shallow"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("Soil water content (m"^"3", " m"^"-3", ")", sep = ""))))
# # abline(a = coef(mSFa)[1], b = coef(mSFa)[2], col = "#3B1F2B", lwd = 1, lty = 3)
# legend("bottomleft", bty = "n", legend = c("Full profile         [0.15-1.45]: p = 0.65", "Deep layers       [0.75-1.45]: p = 0.48", "Mid layers         [0.25-0.75]: p = 0.63", "Shallow layers [0.15-0.25]: p = 0.66"), fill = c("#3B1F2B", "#E15554", "#F6AE2D", "#2E86AB"), text.font = 4, cex = 1.2)
# 
# plot(-500, xlim = c(2011,2024), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(rew_mean ~ y, dfSY[layer == "full"], col = "#3B1F2B", type = "l", lwd = 2)
# points(rew_mean ~ y, dfSY[layer == "deep"], col = "#E15554", type = "l", lwd = 2)
# points(rew_mean ~ y, dfSY[layer == "mid"], col = "#F6AE2D", type = "l", lwd = 2)
# points(rew_mean ~ y, dfSY[layer == "shallow"], col = "#2E86AB", type = "l", lwd = 2)
# points(rew_mean ~ y, dfSY[layer == "full"], pch = 21, bg = "#3B1F2B", type = "p", cex = 1.4)
# points(rew_mean ~ y, dfSY[layer == "deep"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
# points(rew_mean ~ y, dfSY[layer == "mid"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
# points(rew_mean ~ y, dfSY[layer == "shallow"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("Relative extractable water (%)", sep = ""))))
# legend("bottomleft", bty = "n", legend = c("Full profile         [0.15-1.45]: p = 0.65", "Deep layers       [0.75-1.45]: p = 0.48", "Mid layers         [0.25-0.75]: p = 0.64", "Shallow layers [0.15-0.25]: p = 0.70"), fill = c("#3B1F2B", "#E15554", "#F6AE2D", "#2E86AB"), text.font = 4, cex = 1.2)
# 
# plot(-500, xlim = c(2011,2024), ylim = c(0,100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(pct_stress ~ y, stress_days[layer == "full"], col = "#3B1F2B", type = "l", lwd = 2)
# points(pct_stress ~ y, stress_days[layer == "deep"], col = "#E15554", type = "l", lwd = 2)
# points(pct_stress ~ y, stress_days[layer == "mid"], col = "#F6AE2D", type = "l", lwd = 2)
# points(pct_stress ~ y, stress_days[layer == "shallow"], col = "#2E86AB", type = "l", lwd = 2)
# points(pct_stress ~ y, stress_days[layer == "full"], pch = 21, bg = "#3B1F2B", type = "p", cex = 1.4)
# points(pct_stress ~ y, stress_days[layer == "deep"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
# points(pct_stress ~ y, stress_days[layer == "mid"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
# points(pct_stress ~ y, stress_days[layer == "shallow"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("% days below threshold (REW < 0.3)", sep = ""))))
# legend("topleft", bty = "n", legend = c("Full profile         [0.15-1.45]: p = 0.59", "Deep layers       [0.75-1.45]: p = 0.69", "Mid layers         [0.25-0.75]: p = 0.24", "Shallow layers [0.15-0.25]: p = 0.90"), fill = c("#3B1F2B", "#E15554", "#F6AE2D", "#2E86AB"), text.font = 4, cex = 1.2)
# 
# 
# plot(-500, xlim = c(2011,2024), ylim = c(0,0.2), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(mean_deficit ~ y, drought_intensity[layer == "full"], col = "#3B1F2B", type = "l", lwd = 2)
# points(mean_deficit ~ y, drought_intensity[layer == "deep"], col = "#E15554", type = "l", lwd = 2)
# points(mean_deficit ~ y, drought_intensity[layer == "mid"], col = "#F6AE2D", type = "l", lwd = 2)
# points(mean_deficit ~ y, drought_intensity[layer == "shallow"], col = "#2E86AB", type = "l", lwd = 2)
# points(mean_deficit ~ y, drought_intensity[layer == "full"], pch = 21, bg = "#3B1F2B", type = "p", cex = 1.4)
# points(mean_deficit ~ y, drought_intensity[layer == "deep"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
# points(mean_deficit ~ y, drought_intensity[layer == "mid"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
# points(mean_deficit ~ y, drought_intensity[layer == "shallow"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# abline(a = coef(mDM)[1], b = coef(mDM)[2], col = "#F6AE2D", lwd = 1, lty = 3)
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("Mean deficit below threshold (REW < 0.3)", sep = ""))))
# legend("topleft", bty = "n", legend = c("Full profile         [0.15-1.45]: p = 0.38", "Deep layers       [0.75-1.45]: p = 0.54", "Mid layers         [0.25-0.75]: p = 0.06 .", "Shallow layers [0.15-0.25]: p = 0.60"), fill = c("#3B1F2B", "#E15554", "#F6AE2D", "#2E86AB"), text.font = 4, cex = 1.2)
# 
# dev.off()

# Partition -------------------------------------------------------------------
