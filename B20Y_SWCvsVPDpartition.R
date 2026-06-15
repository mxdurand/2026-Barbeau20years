library(data.table)

df <- fread("data/B20Y_BM_all.csv")

# To long format
SWC_cols <- names(df)[grep("SWC", names(df))]
df <- df[, c("TIMESTAMP", SWC_cols), with = FALSE]
df <- melt(df, id.vars = "TIMESTAMP", variable.name = "sensor", value.name = "swc")
rm(SWC_cols)

df[, y := substring(TIMESTAMP, 1 , 4)]
df[, m := substring(TIMESTAMP, 6 , 7)]
df[, d := substring(TIMESTAMP, 9 , 10)]
df[, h := substring(TIMESTAMP, 12 , 16)]
date_lookup <- unique(df[, .(y, m, d)])
date_lookup[, doy := yday(as.Date(paste(y, m, d, sep = "-")))]
df <- merge(df, date_lookup, by = c("y", "m", "d"))
rm(date_lookup)

df[, sensor := as.character(sensor)]
df[, h_pos  := as.integer(sub("SWC_([0-9]+)_([0-9]+)_1", "\\1", sensor))]
df[, v_dep  := as.integer(sub("SWC_([0-9]+)_([0-9]+)_1", "\\2", sensor))]

d0 <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45)
depth_map <- data.table(v_dep = 1:15, depth = d0)
df <- merge(df, depth_map, by = "v_dep")
rm(depth_map)

# Standardize by position
pos_means <- df[, .(pos_mean = mean(swc, na.rm = TRUE)), by = .(h_pos, depth)]
df <- merge(df, pos_means, by = c("h_pos", "depth"))
df[, swc_anom := swc - pos_mean]
rm(pos_means)

# Remove missing year
df <- df[depth != 0.05]

# Tracking how many positions contribute
dfST <- df[, .(
  swc_mean  = mean(swc_anom, na.rm = TRUE), # position-corrected mean
  swc_abs   = mean(swc, na.rm = TRUE),      # raw absolute value if needed
  n_pos     = sum(!is.na(swc))              # now this is number of positions (max 6)
), by = .(y, m, d, doy, h, depth)]

# Add layer definition
dfST[depth <= 0.25, layer := "shallow"]
dfST[depth > 0.25 & depth <= 0.75, layer := "mid"]
dfST[depth > 0.75, layer := "deep"]

# every hour by full profile
dfSTD <- dfST[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),  
  swc_abs   = mean(swc_abs, na.rm = TRUE),        
  n_depths  = sum(!is.na(swc_mean))  # how many depths contributed
), by = .(y, m, d, doy, h)]


# Layer means with coverage tracking
dfSL <- dfST[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),
  swc_abs   = mean(swc_abs,  na.rm = TRUE),
  n_depths  = sum(!is.na(swc_mean)),  # how many depths contributed
  n_pos_mean = mean(n_pos, na.rm = TRUE)   # average positions per timestep
), by = .(y, m, d, doy, h, layer)]


# Daily data
dfSJ <- dfST[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),  
  swc_abs   = mean(swc_abs, na.rm = TRUE),        
  n_hours    = sum(!is.na(swc_mean))        # valid half-hours in the day
), by = .(y, m, d, doy, layer)]

min_hours <- 24   # at least 24 out of 48 half-hours
dfSJ[n_hours < min_hours, swc_mean := NA]
dfSJ[n_hours < min_hours, swc_abs  := NA]

dfSD <- dfSJ[, .(
  swc_mean  = mean(swc_mean, na.rm = TRUE),  
  swc_abs   = mean(swc_abs, na.rm = TRUE),        
  n_depths  = sum(!is.na(swc_mean))  # how many depths contributed
), by = .(y, m, d, doy)]

rew_params <- dfST[, .(
  FC = quantile(swc_abs, 0.95, na.rm = TRUE),
  WP = quantile(swc_abs, 0.05, na.rm = TRUE)
), by = layer]

dfSJ <- merge(dfSJ, rew_params, by = "layer")
dfSD[, FC := quantile(dfSTD$swc_abs, 0.95, na.rm = TRUE)]
dfSD[, WP := quantile(dfSTD$swc_abs, 0.05, na.rm = TRUE)]

dfSJ[, REW := (swc_abs - WP) / (FC - WP)]
dfSD[, REW := (swc_abs - WP) / (FC - WP)]
dfSJ[REW < 0, REW := 0] ; dfSJ[REW > 1, REW := 1]
dfSD[REW < 0, REW := 0] ; dfSD[REW > 1, REW := 1]

# Get growing season data
dfGS <- fread("data/df_inGS.csv")
dfGS[, y := as.character(y)]
dfGS[, m := as.character(m)]
dfGS[, d := as.character(d)]
dfSJ <- merge(dfSJ, dfGS, by = c("y", "m", "d", "doy"), all.x = T)
dfSD <- merge(dfSD, dfGS, by = c("y", "m", "d", "doy"), all.x = T)

dfSLY <- dfSJ[in_gs == TRUE & !is.na(swc_mean), .(
  swc_mean = mean(swc_mean, na.rm = TRUE),
  swc_abs  = mean(swc_abs,  na.rm = TRUE),
  rew_mean  = mean(REW,  na.rm = TRUE),
  n_days   = sum(!is.na(swc_mean))
), by = .(y, layer)]

dfSY <- dfSJ[in_gs == TRUE & !is.na(swc_mean), .(
  swc_mean = mean(swc_mean, na.rm = TRUE),
  swc_abs  = mean(swc_abs,  na.rm = TRUE),
  rew_mean  = mean(REW,  na.rm = TRUE),
  n_days   = sum(!is.na(swc_mean))
), by = .(y)]


dfSLY[, y := as.numeric(y)]
dfSY[, y := as.numeric(y)]

mSS <- lm(swc_mean ~ y, dfSLY[layer == "shallow"]) ; summary(mSS)
mSM <- lm(swc_mean ~ y, dfSLY[layer == "mid"]) ; summary(mSM)
mSD <- lm(swc_mean ~ y, dfSLY[layer == "deep"]) ; summary(mSD)
mSA <- lm(swc_mean ~ y, dfSY) ; summary(mSA)
mSSa <- lm(swc_abs ~ y, dfSLY[layer == "shallow"]) ; summary(mSSa)
mSMa <- lm(swc_abs ~ y, dfSLY[layer == "mid"]) ; summary(mSMa)
mSDa <- lm(swc_abs ~ y, dfSLY[layer == "deep"]) ; summary(mSDa)
mSAa <- lm(swc_abs ~ y, dfSY) ; summary(mSAa)
mRS <- lm(rew_mean ~ y, dfSLY[layer == "shallow"]) ; summary(mRS)
mRM <- lm(rew_mean ~ y, dfSLY[layer == "mid"]) ; summary(mRM)
mRD <- lm(rew_mean ~ y, dfSLY[layer == "deep"]) ; summary(mRD)
mRA <- lm(rew_mean ~ y, dfSY) ; summary(mRA)

plot(swc_mean ~ y, dfSLY[layer == "shallow"]) ; abline(a = coef(mSS)[1], b = coef(mSS)[2], col = "red")
plot(swc_mean ~ y, dfSLY[layer == "mid"]) ; abline(a = coef(mSM)[1], b = coef(mSM)[2], col = "red")
plot(swc_mean ~ y, dfSLY[layer == "deep"]) ; abline(a = coef(mSD)[1], b = coef(mSD)[2], col = "red")
plot(swc_mean ~ y, dfSY) ; abline(a = coef(mSA)[1], b = coef(mSA)[2], col = "red")

plot(swc_abs ~ y, dfSLY[layer == "shallow"]) ; abline(a = coef(mSSa)[1], b = coef(mSSa)[2], col = "red")
plot(swc_abs ~ y, dfSLY[layer == "mid"]) ; abline(a = coef(mSMa)[1], b = coef(mSMa)[2], col = "red")
plot(swc_abs ~ y, dfSLY[layer == "deep"]) ; abline(a = coef(mSDa)[1], b = coef(mSDa)[2], col = "red")
plot(swc_abs ~ y, dfSY) ; abline(a = coef(mSAa)[1], b = coef(mSAa)[2], col = "red")

plot(rew_mean ~ y, dfSLY[layer == "shallow"]) ; abline(a = coef(mRS)[1], b = coef(mRS)[2], col = "red")
plot(rew_mean ~ y, dfSLY[layer == "mid"]) ; abline(a = coef(mRM)[1], b = coef(mRM)[2], col = "red")
plot(rew_mean ~ y, dfSLY[layer == "deep"]) ; abline(a = coef(mRD)[1], b = coef(mRD)[2], col = "red")
plot(rew_mean ~ y, dfSY) ; abline(a = coef(mRA)[1], b = coef(mRA)[2], col = "red")









#       
#       
# # -------------------------------------------------------------------------
# 
# # Define layers
# d1 <- c(0.05, 0.15, 0.25)                             # V = 1,2,3
# d2 <- c(0.35, 0.45, 0.55, 0.65, 0.75)                 # V = 4,5,6,7,8
# d3 <- c(0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45)     # V = 9-15
# d0 <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45)
# pos <- 1:6  # which horizontal positions to use
# 
# # Get layer columns
# get_swc_cols <- function(h_pos, depth_vec, all_depths = d0) {
#   v_idx <- which(all_depths %in% depth_vec)
#   as.vector(outer(h_pos, v_idx, function(h, v) paste0("SWC_", h, "_", v, "_1")))
# }
# 
# cols0 <- get_swc_cols(pos, d0)  # full profile
# cols1 <- get_swc_cols(pos, d1)  # Shallow
# cols2 <- get_swc_cols(pos, d2)  # Mid
# cols3 <- get_swc_cols(pos, d3)  # Deep
# 
# # Verify columns exist in df
# cols1 <- intersect(cols1, names(df))
# cols2 <- intersect(cols2, names(df))
# cols3 <- intersect(cols3, names(df))
# cols0 <- intersect(cols0, names(df))
# 
# # cat("Shallow:", length(cols1), "columns\n")
# # cat("Mid:",     length(cols2),     "columns\n")
# # cat("Deep:",    length(cols3),    "columns\n")
# # cat("Full profile:", length(cols0), "columns\n")
# 
# # Compute layer means
# df[, SWC1 := rowMeans(.SD, na.rm = TRUE), .SDcols = cols1]
# df[, SWC2 := rowMeans(.SD, na.rm = TRUE), .SDcols = cols2]
# df[, SWC3 := rowMeans(.SD, na.rm = TRUE), .SDcols = cols3]
# df[, SWC0 := rowMeans(.SD, na.rm = TRUE), .SDcols = cols0]
# 
# # -------------------------------------------------------------------------
# 
# # Compare position 6 vs others where both are available
# df[!is.na(SWC_6_1_1), .(
#   mean_pos1_4 = mean(rowMeans(.SD[, get_swc_cols(1:4, d1) |> 
#                                     intersect(names(df)), with = FALSE], na.rm = TRUE)),
#   mean_pos6   = mean(rowMeans(.SD[, get_swc_cols(6, d1) |> 
#                                     intersect(names(df)), with = FALSE], na.rm = TRUE))
# )]
# 
# # Simpler version
# pos1_4_mean <- df[, rowMeans(.SD, na.rm = TRUE), 
#                   .SDcols = intersect(get_swc_cols(1:4, d0), names(df))]
# pos6_mean   <- df[, rowMeans(.SD, na.rm = TRUE), 
#                   .SDcols = intersect(get_swc_cols(6, d0), names(df))]
# 
# # Are they systematically different?
# mean(pos6_mean - pos1_4_mean, na.rm = TRUE)
