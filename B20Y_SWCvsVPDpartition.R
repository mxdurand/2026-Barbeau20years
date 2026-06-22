library(data.table)
library(ppcor)
library(rdacca.hp)
library(mgcv)

source("B20Y_PETEC.R")
source("B20Y_SWCannualTrends.R")

dfSJ[, m := as.integer(m)] ; dfSJ[, d := as.integer(d)]
dfPJ <- merge(dfEJ, dfSJ, by = c("y", "m", "d", "doy"))
dfPJ <- dfPJ[in_gs.x == TRUE]

### IMPORTANT Restrincting dataset (otherwise confounding factors drive results)

dfPJ <- dfPJ[layer == "full"]   # Taking the full SWC profile
dfPJ <- dfPJ[rn > 10]           # Only bright days (remove VPD dependence on radiation + some pheno effect)
dfPJ <- dfPJ[m %in% c(6,7,8)]   # Keep only month where canopy fully developped (wet and dry days have pheno dependence)

# Create bins
dfPJ[, REW_tercile := cut(REW, breaks = quantile(REW, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          labels = c("dry", "mod", "wet"), include.lowest = TRUE)]
dfPJ[, vpd_bin := cut(vpd, breaks = quantile(vpd, probs = seq(0, 1, 0.1), na.rm = TRUE),
                      include.lowest = TRUE)]

# Subset
# dfPJ <- dfPJ[y < 2018]   # Subset for early period
dfPJ <- dfPJ[y >= 2018]  # Subset for late period


# Explaining the early vs late --------------------------------------------

dfPJ_full <- dfPJ
dfPJ_early <- dfPJ[y < 2018]
dfPJ_late <- dfPJ[y >= 2018]

# No change in correlative power
cor.test(dfPJ_full[,vpd], dfPJ_full[,REW])
cor.test(dfPJ_early[,vpd], dfPJ_early[,REW])
cor.test(dfPJ_late[,vpd], dfPJ_late[,REW])

# A. More and stronger drought in later period

dfPJ_early[REW < 0.3, .N] # Doubled no days with REW < 0.3
dfPJ_late[REW < 0.3, .N]

dfPJ_early[REW_tercile == "dry", mean(REW)] # Stress intensit a bit stronger in later years
dfPJ_late[REW_tercile == "dry", mean(REW)]

# B. GPP-REW slope in increasing in late period 
fit_early <- lm(gpp ~ vpd + REW, data = dfPJ_early)
fit_late  <- lm(gpp ~ vpd + REW, data = dfPJ_late)

coef(fit_early)["REW"]  # 4.6
coef(fit_late)["REW"]   # 7.3 35% increase

# Is the difference significant?
dfPJ_full[, period := fifelse(y < 2018, "early", "late")]
fit_interaction <- lm(gpp ~ vpd + REW * period, data = dfPJ_full)
summary(fit_interaction)  # REW:periodlate is significant: change in slope significant

fit_full <- lm(gpp ~ vpd + REW * period, data = dfPJ_full)
fit_et  <- lm(e   ~ vpd + REW * period, data = dfPJ_full)
fit_ter <- lm(ter ~ vpd + REW * period, data = dfPJ_full)
summary(fit_full)
summary(fit_et)
summary(fit_ter)

dfPJ_full[, .(
  rew_mean = mean(REW, na.rm = TRUE),
  rew_p10  = quantile(REW, 0.10, na.rm = TRUE),
  rew_p25  = quantile(REW, 0.25, na.rm = TRUE)
), by = .(period, REW_tercile)][order(REW_tercile, period)]

### Summary:
# A. [Minor] More extreme drought days: 84% more days with REW <0.3 but stress not much stronger
# B. [Major] Steeper flux-REW sensitivity: Slope almst doubled for all fluxes, highly signif. so same SWD produces larger flux reduction
# C. [unexpexted]: All fluxes are lowed in late period: thinning????

# Binned response curve ---------------------------------------------------

dfPJ[, REW_tercile := cut(REW, breaks = quantile(REW, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          labels = c("dry", "mod", "wet"), include.lowest = TRUE)]
dfPJ[, vpd_bin := cut(vpd, breaks = quantile(vpd, probs = seq(0, 1, 0.1), na.rm = TRUE),
                      include.lowest = TRUE)]
# Binned response curves
dfB <- dfPJ[, .(
  gpp_mean = mean(gpp, na.rm = TRUE),
  gpp_se   = sd(gpp, na.rm = TRUE) / sqrt(.N),
  gpp_sd   = sd(gpp, na.rm = TRUE),
  et_mean = mean(e, na.rm = TRUE),
  et_se   = sd(e, na.rm = TRUE) / sqrt(.N),
  et_sd   = sd(e, na.rm = TRUE),
  ter_mean = mean(ter, na.rm = TRUE),
  ter_se   = sd(ter, na.rm = TRUE) / sqrt(.N),
  ter_sd   = sd(ter, na.rm = TRUE),
  nee_mean = mean(nee, na.rm = TRUE),
  nee_se   = sd(nee, na.rm = TRUE) / sqrt(.N),
  nee_sd   = sd(nee, na.rm = TRUE),
  vpd_mean = mean(vpd, na.rm = TRUE),
  vpd_se   = sd(vpd, na.rm = TRUE) / sqrt(.N),
  vpd_sd   = sd(vpd, na.rm = TRUE),
  n        = .N
), by = .(REW_tercile, vpd_bin)]
setorder(dfB, REW_tercile, vpd_mean)


# Variance partitioning ---------------------------------------------------

# GPP
mf_gpp <- lm(gpp ~ vpd + REW, data = dfPJ)     # Full model
mv_gpp <- lm(gpp ~ vpd, data = dfPJ)           # VPD only
mr_gpp <- lm(gpp ~ REW, data = dfPJ)           # REW only
rf_gpp <- summary(mf_gpp)$r.squared            # R2 of full model
rv_gpp <- summary(mv_gpp)$r.squared            # R2 of VPD only model
rr_gpp <- summary(mr_gpp)$r.squared            # R2 of REW only model
uv_gpp <- rf_gpp - rr_gpp                      # Contribution of VPD
ur_gpp <- rf_gpp - rv_gpp                      # Contribution of REW
sh_gpp <- rv_gpp + rr_gpp - rf_gpp             # Contribution shared

# ETR
mf_etr <- lm(e ~ vpd + REW, data = dfPJ)       # Full model
mv_etr <- lm(e ~ vpd, data = dfPJ)             # VPD only
mr_etr <- lm(e ~ REW, data = dfPJ)             # REW only
rf_etr <- summary(mf_etr)$r.squared            # R2 of full model
rv_etr <- summary(mv_etr)$r.squared            # R2 of VPD only model
rr_etr <- summary(mr_etr)$r.squared            # R2 of REW only model
uv_etr <- rf_etr - rr_etr                      # Contribution of VPD
ur_etr <- rf_etr - rv_etr                      # Contribution of REW
sh_etr <- rv_etr + rr_etr - rf_etr             # Contribution shared

# TER
mf_ter <- lm(ter ~ vpd + REW, data = dfPJ)     # Full model
mv_ter <- lm(ter ~ vpd, data = dfPJ)           # VPD only
mr_ter <- lm(ter ~ REW, data = dfPJ)           # REW only
rf_ter <- summary(mf_ter)$r.squared            # R2 of full model
rv_ter <- summary(mv_ter)$r.squared            # R2 of VPD only model
rr_ter <- summary(mr_ter)$r.squared            # R2 of REW only model
uv_ter <- rf_ter - rr_ter                      # Contribution of VPD
ur_ter <- rf_ter - rv_ter                      # Contribution of REW
sh_ter <- rv_ter + rr_ter - rf_ter             # Contribution shared

# NEE
mf_nee <- lm(nee ~ vpd + REW, data = dfPJ)     # Full model
mv_nee <- lm(nee ~ vpd, data = dfPJ)           # VPD only
mr_nee <- lm(nee ~ REW, data = dfPJ)           # REW only
rf_nee <- summary(mf_nee)$r.squared            # R2 of full model
rv_nee <- summary(mv_nee)$r.squared            # R2 of VPD only model
rr_nee <- summary(mr_nee)$r.squared            # R2 of REW only model
uv_nee <- rf_nee - rr_nee                      # Contribution of VPD
ur_nee <- rf_nee - rv_nee                      # Contribution of REW
sh_nee <- rv_nee + rr_nee - rf_nee             # Contribution shared

dtc <- data.table("nee" = 100*round(c(rf_nee, uv_nee, ur_nee, sh_nee), 3),
                  "gpp" = 100*round(c(rf_gpp, uv_gpp, ur_gpp, sh_gpp), 3),
                  "etr" = 100*round(c(rf_etr, uv_etr, ur_etr, sh_etr), 3),
                  "ter" = 100*round(c(rf_ter, uv_ter, ur_ter, sh_ter), 3))
row.names(dtc) <- c("tot", "vpd", "rew", "sha") 

fwrite(dtc, "data/output/PartitionSWC_VPD_RN10_JJA_late.csv", row.names = TRUE)

# Plots ---------------------------------------------------------------------

# Calculate slopes of the linear parts
get_slopes <- function(target_var, data) {
  data[, {
    m <- lm(get(target_var) ~ vpd)
    .(target = target_var,
      slope  = coef(m)[2],
      int    = coef(m)[1],
      se     = summary(m)$coefficients["vpd", "Std. Error"],
      p      = summary(m)$coefficients["vpd", "Pr(>|t|)"])
  }, by = REW_tercile]
}
targets <- c("gpp", "e", "ter", "nee")
slopes <- rbindlist(lapply(targets, get_slopes, data = dfPJ))
x_vpd <- seq(0.25,1.9,0.01)

# Significance
getSignif <- function(p){
  if(p < 0.001){
    s <- "p < 0.001 ***"
  } else if (p < 0.01){
    s <- "p < 0.01 **"
  } else if (p < 0.05){
    s <- "p < 0.05 *"
  } else {
    s <- paste0(round(p,1), " ns") 
  }
  return(s)
}

### GPP
png("figs/B20Y_PartitionSWC_VPD_RN10_JJA_late.png", height = 5000, width = 15000, res = 1000)
par(mfrow = c(1,3), bty = "L", mar = c(5,4,1,1), oma = c(0,1,0,0))

plot(-500, xlim = c(0,2), ylim = c(0,18), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
if(slopes[REW_tercile == "wet" & target == "gpp", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "gpp",slope]+slopes[REW_tercile == "wet" & target == "gpp",int], col = "#2E86AB", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "gpp",slope]+slopes[REW_tercile == "wet" & target == "gpp",int], col = "#2E86AB", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "mod" & target == "gpp", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "gpp",slope]+slopes[REW_tercile == "mod" & target == "gpp",int], col = "#F6AE2D", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "gpp",slope]+slopes[REW_tercile == "mod" & target == "gpp",int], col = "#F6AE2D", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "dry" & target == "gpp", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "gpp",slope]+slopes[REW_tercile == "dry" & target == "gpp",int], col = "#E15554", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "gpp",slope]+slopes[REW_tercile == "dry" & target == "gpp",int], col = "#E15554", lty = 3, lwd = 2)
}
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "wet"], col = "#2E86AB", type = "l", lwd = 1, lty = 3)
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "mod"], col = "#F6AE2D", type = "l", lwd = 1, lty = 3)
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "dry"], col = "#E15554", type = "l", lwd = 1, lty = 3)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean],  x1 = dfB[REW_tercile == "wet", vpd_mean], y0 = dfB[REW_tercile == "wet", gpp_mean] + dfB[REW_tercile == "wet", gpp_se], y1 = dfB[REW_tercile == "wet", gpp_mean] - dfB[REW_tercile == "wet", gpp_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean],  x1 = dfB[REW_tercile == "mod", vpd_mean], y0 = dfB[REW_tercile == "mod", gpp_mean] + dfB[REW_tercile == "mod", gpp_se], y1 = dfB[REW_tercile == "mod", gpp_mean] - dfB[REW_tercile == "mod", gpp_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean],  x1 = dfB[REW_tercile == "dry", vpd_mean], y0 = dfB[REW_tercile == "dry", gpp_mean] + dfB[REW_tercile == "dry", gpp_se], y1 = dfB[REW_tercile == "dry", gpp_mean] - dfB[REW_tercile == "dry", gpp_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean] + dfB[REW_tercile == "wet", vpd_se],  x1 = dfB[REW_tercile == "wet", vpd_mean] - dfB[REW_tercile == "wet", vpd_se], y0 = dfB[REW_tercile == "wet", gpp_mean], y1 = dfB[REW_tercile == "wet", gpp_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean] + dfB[REW_tercile == "mod", vpd_se],  x1 = dfB[REW_tercile == "mod", vpd_mean] - dfB[REW_tercile == "mod", vpd_se], y0 = dfB[REW_tercile == "mod", gpp_mean], y1 = dfB[REW_tercile == "mod", gpp_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean] + dfB[REW_tercile == "dry", vpd_se],  x1 = dfB[REW_tercile == "dry", vpd_mean] - dfB[REW_tercile == "dry", vpd_se], y0 = dfB[REW_tercile == "dry", gpp_mean], y1 = dfB[REW_tercile == "dry", gpp_mean], code = 3, angle = 90, length = 0.03)
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "wet"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "mod"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
points(gpp_mean ~ vpd_mean, dfB[REW_tercile == "dry"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "VPD (kPa)")
mtext(side = 2, font = 4, line = 2.3, cex = 1.4, text = expression(bolditalic(paste("GPP (gC m"^"-2", "day"^"-1", ")", sep = ""))))
legend("bottomright", bty = "n", title = "Linear part only:", legend = c(paste0("Wet: ", getSignif(slopes[REW_tercile == "wet" & target == "gpp", p])), paste0("Mid: ", getSignif(slopes[REW_tercile == "mod" & target == "gpp", p])), paste0("Dry: ", getSignif(slopes[REW_tercile == "dry" & target == "gpp", p]))), fill = c("#2E86AB", "#F6AE2D", "#E15554"), text.font = 4, cex = 1.2)
legend(x = -0.1, y = 17.9, bty = "n", legend = c("Variance decomposition (full data)"), text.font = 4, cex = 1)
legend(x = -0.1, y = 17.3, bty = "n", legend = c(paste("Total", dtc[1,gpp], "%"), paste("VPD", dtc[2,gpp], "%"), paste("REW", dtc[3,gpp], "%"), paste("Shared", dtc[4,gpp], "%")), text.font = 1, cex = 1)
title(main = "Days split into equal-sized REW/VPD bins")


### ETR
plot(-500, xlim = c(0,2), ylim = c(0,6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "L")
if(slopes[REW_tercile == "wet" & target == "e", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "e",slope]+slopes[REW_tercile == "wet" & target == "e",int], col = "#2E86AB", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "e",slope]+slopes[REW_tercile == "wet" & target == "e",int], col = "#2E86AB", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "mod" & target == "e", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "e",slope]+slopes[REW_tercile == "mod" & target == "e",int], col = "#F6AE2D", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "e",slope]+slopes[REW_tercile == "mod" & target == "e",int], col = "#F6AE2D", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "dry" & target == "e", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "e",slope]+slopes[REW_tercile == "dry" & target == "e",int], col = "#E15554", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "e",slope]+slopes[REW_tercile == "dry" & target == "e",int], col = "#E15554", lty = 3, lwd = 2)
}
points(et_mean ~ vpd_mean, dfB[REW_tercile == "wet"], col = "#2E86AB", type = "l", lwd = 1, lty = 3)
points(et_mean ~ vpd_mean, dfB[REW_tercile == "mod"], col = "#F6AE2D", type = "l", lwd = 1, lty = 3)
points(et_mean ~ vpd_mean, dfB[REW_tercile == "dry"], col = "#E15554", type = "l", lwd = 1, lty = 3)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean],  x1 = dfB[REW_tercile == "wet", vpd_mean], y0 = dfB[REW_tercile == "wet", et_mean] + dfB[REW_tercile == "wet", et_se], y1 = dfB[REW_tercile == "wet", et_mean] - dfB[REW_tercile == "wet", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean],  x1 = dfB[REW_tercile == "mod", vpd_mean], y0 = dfB[REW_tercile == "mod", et_mean] + dfB[REW_tercile == "mod", et_se], y1 = dfB[REW_tercile == "mod", et_mean] - dfB[REW_tercile == "mod", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean],  x1 = dfB[REW_tercile == "dry", vpd_mean], y0 = dfB[REW_tercile == "dry", et_mean] + dfB[REW_tercile == "dry", et_se], y1 = dfB[REW_tercile == "dry", et_mean] - dfB[REW_tercile == "dry", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean] + dfB[REW_tercile == "wet", vpd_se],  x1 = dfB[REW_tercile == "wet", vpd_mean] - dfB[REW_tercile == "wet", vpd_se], y0 = dfB[REW_tercile == "wet", et_mean], y1 = dfB[REW_tercile == "wet", et_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean] + dfB[REW_tercile == "mod", vpd_se],  x1 = dfB[REW_tercile == "mod", vpd_mean] - dfB[REW_tercile == "mod", vpd_se], y0 = dfB[REW_tercile == "mod", et_mean], y1 = dfB[REW_tercile == "mod", et_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean] + dfB[REW_tercile == "dry", vpd_se],  x1 = dfB[REW_tercile == "dry", vpd_mean] - dfB[REW_tercile == "dry", vpd_se], y0 = dfB[REW_tercile == "dry", et_mean], y1 = dfB[REW_tercile == "dry", et_mean], code = 3, angle = 90, length = 0.03)
points(et_mean ~ vpd_mean, dfB[REW_tercile == "wet"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
points(et_mean ~ vpd_mean, dfB[REW_tercile == "mod"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
points(et_mean ~ vpd_mean, dfB[REW_tercile == "dry"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "VPD (kPa)")
mtext(side = 2, font = 4, line = 2.3, cex = 1.4, text = expression(bolditalic(paste("ET (mm day"^"-1", ")", sep = ""))))
legend("bottomright", bty = "n", title = "Linear part only:", legend = c(paste0("Wet: ", getSignif(slopes[REW_tercile == "wet" & target == "e", p])), paste0("Mid: ", getSignif(slopes[REW_tercile == "mod" & target == "e", p])), paste0("Dry: ", getSignif(slopes[REW_tercile == "dry" & target == "e", p]))), fill = c("#2E86AB", "#F6AE2D", "#E15554"), text.font = 4, cex = 1.2)
legend(x = -0.1, y = 6, bty = "n", legend = c("Variance decomposition (full data)"), text.font = 4, cex = 1)
legend(x = -0.1, y = 5.5, bty = "n", legend = c(paste("Total", dtc[1,etr], "%"), paste("VPD", dtc[2,etr], "%"), paste("REW", dtc[3,etr], "%"), paste("Shared", dtc[4,etr], "%")), text.font = 1, cex = 1)
title(main = "Days split into equal-sized REW/VPD bins")


### TER
plot(-500, xlim = c(0,2), ylim = c(0,10), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "L")
if(slopes[REW_tercile == "wet" & target == "ter", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "ter",slope]+slopes[REW_tercile == "wet" & target == "ter",int], col = "#2E86AB", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "ter",slope]+slopes[REW_tercile == "wet" & target == "ter",int], col = "#2E86AB", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "mod" & target == "ter", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "ter",slope]+slopes[REW_tercile == "mod" & target == "ter",int], col = "#F6AE2D", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "ter",slope]+slopes[REW_tercile == "mod" & target == "ter",int], col = "#F6AE2D", lty = 3, lwd = 2)
}
if(slopes[REW_tercile == "dry" & target == "ter", p] < 0.05){
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "ter",slope]+slopes[REW_tercile == "dry" & target == "ter",int], col = "#E15554", lty = 1, lwd = 2)
} else {
  lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "ter",slope]+slopes[REW_tercile == "dry" & target == "ter",int], col = "#E15554", lty = 3, lwd = 2)
}
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "wet"], col = "#2E86AB", type = "l", lwd = 1, lty = 3)
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "mod"], col = "#F6AE2D", type = "l", lwd = 1, lty = 3)
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "dry"], col = "#E15554", type = "l", lwd = 1, lty = 3)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean],  x1 = dfB[REW_tercile == "wet", vpd_mean], y0 = dfB[REW_tercile == "wet", ter_mean] + dfB[REW_tercile == "wet", et_se], y1 = dfB[REW_tercile == "wet", ter_mean] - dfB[REW_tercile == "wet", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean],  x1 = dfB[REW_tercile == "mod", vpd_mean], y0 = dfB[REW_tercile == "mod", ter_mean] + dfB[REW_tercile == "mod", et_se], y1 = dfB[REW_tercile == "mod", ter_mean] - dfB[REW_tercile == "mod", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean],  x1 = dfB[REW_tercile == "dry", vpd_mean], y0 = dfB[REW_tercile == "dry", ter_mean] + dfB[REW_tercile == "dry", et_se], y1 = dfB[REW_tercile == "dry", ter_mean] - dfB[REW_tercile == "dry", et_se], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "wet", vpd_mean] + dfB[REW_tercile == "wet", vpd_se],  x1 = dfB[REW_tercile == "wet", vpd_mean] - dfB[REW_tercile == "wet", vpd_se], y0 = dfB[REW_tercile == "wet", ter_mean], y1 = dfB[REW_tercile == "wet", ter_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "mod", vpd_mean] + dfB[REW_tercile == "mod", vpd_se],  x1 = dfB[REW_tercile == "mod", vpd_mean] - dfB[REW_tercile == "mod", vpd_se], y0 = dfB[REW_tercile == "mod", ter_mean], y1 = dfB[REW_tercile == "mod", ter_mean], code = 3, angle = 90, length = 0.03)
arrows(x0 = dfB[REW_tercile == "dry", vpd_mean] + dfB[REW_tercile == "dry", vpd_se],  x1 = dfB[REW_tercile == "dry", vpd_mean] - dfB[REW_tercile == "dry", vpd_se], y0 = dfB[REW_tercile == "dry", ter_mean], y1 = dfB[REW_tercile == "dry", ter_mean], code = 3, angle = 90, length = 0.03)
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "wet"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "mod"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
points(ter_mean ~ vpd_mean, dfB[REW_tercile == "dry"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "VPD (kPa)")
mtext(side = 2, font = 4, line = 2.3, cex = 1.4, text = expression(bolditalic(paste("TER (gC m"^"-2", "day"^"-1", ")", sep = ""))))
legend("bottomright", bty = "n", title = "Linear part only:", legend = c(paste0("Wet: ", getSignif(slopes[REW_tercile == "wet" & target == "ter", p])), paste0("Mid: ", getSignif(slopes[REW_tercile == "mod" & target == "ter", p])), paste0("Dry: ", getSignif(slopes[REW_tercile == "dry" & target == "ter", p]))), fill = c("#2E86AB", "#F6AE2D", "#E15554"), text.font = 4, cex = 1.2)
legend(x = -0.1, y = 10, bty = "n", legend = c("Variance decomposition (full data)"), text.font = 4, cex = 1)
legend(x = -0.1, y = 9.5, bty = "n", legend = c(paste("Total", dtc[1,ter], "%"), paste("VPD", dtc[2,ter], "%"), paste("REW", dtc[3,ter], "%"), paste("Shared", dtc[4,ter], "%")), text.font = 1, cex = 1)
title(main = "Days split into equal-sized REW/VPD bins")

dev.off()

### NEE
# x_vpd <- seq(0.06,1,0.01)
# plot(-500, xlim = c(0,1.8), ylim = c(0,-9), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "L")
# lines(x_vpd, x_vpd*slopes[REW_tercile == "wet" & target == "nee",slope]+slopes[REW_tercile == "wet" & target == "nee",int], col = "#2E86AB", lty = 1, lwd = 2)
# lines(x_vpd, x_vpd*slopes[REW_tercile == "mod" & target == "nee",slope]+slopes[REW_tercile == "mod" & target == "nee",int], col = "#F6AE2D", lty = 1, lwd = 2)
# lines(x_vpd, x_vpd*slopes[REW_tercile == "dry" & target == "nee",slope]+slopes[REW_tercile == "dry" & target == "nee",int], col = "#E15554", lty = 1, lwd = 2)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "wet"], col = "#2E86AB", type = "l", lwd = 1, lty = 3)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "mod"], col = "#F6AE2D", type = "l", lwd = 1, lty = 3)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "dry"], col = "#E15554", type = "l", lwd = 1, lty = 3)
# arrows(x0 = dfB[REW_tercile == "wet", vpd_mean],  x1 = dfB[REW_tercile == "wet", vpd_mean], y0 = dfB[REW_tercile == "wet", nee_mean] + dfB[REW_tercile == "wet", et_se], y1 = dfB[REW_tercile == "wet", nee_mean] - dfB[REW_tercile == "wet", et_se], code = 3, angle = 90, length = 0.03)
# arrows(x0 = dfB[REW_tercile == "mod", vpd_mean],  x1 = dfB[REW_tercile == "mod", vpd_mean], y0 = dfB[REW_tercile == "mod", nee_mean] + dfB[REW_tercile == "mod", et_se], y1 = dfB[REW_tercile == "mod", nee_mean] - dfB[REW_tercile == "mod", et_se], code = 3, angle = 90, length = 0.03)
# arrows(x0 = dfB[REW_tercile == "dry", vpd_mean],  x1 = dfB[REW_tercile == "dry", vpd_mean], y0 = dfB[REW_tercile == "dry", nee_mean] + dfB[REW_tercile == "dry", et_se], y1 = dfB[REW_tercile == "dry", nee_mean] - dfB[REW_tercile == "dry", et_se], code = 3, angle = 90, length = 0.03)
# arrows(x0 = dfB[REW_tercile == "wet", vpd_mean] + dfB[REW_tercile == "wet", vpd_se],  x1 = dfB[REW_tercile == "wet", vpd_mean] - dfB[REW_tercile == "wet", vpd_se], y0 = dfB[REW_tercile == "wet", nee_mean], y1 = dfB[REW_tercile == "wet", nee_mean], code = 3, angle = 90, length = 0.03)
# arrows(x0 = dfB[REW_tercile == "mod", vpd_mean] + dfB[REW_tercile == "mod", vpd_se],  x1 = dfB[REW_tercile == "mod", vpd_mean] - dfB[REW_tercile == "mod", vpd_se], y0 = dfB[REW_tercile == "mod", nee_mean], y1 = dfB[REW_tercile == "mod", nee_mean], code = 3, angle = 90, length = 0.03)
# arrows(x0 = dfB[REW_tercile == "dry", vpd_mean] + dfB[REW_tercile == "dry", vpd_se],  x1 = dfB[REW_tercile == "dry", vpd_mean] - dfB[REW_tercile == "dry", vpd_se], y0 = dfB[REW_tercile == "dry", nee_mean], y1 = dfB[REW_tercile == "dry", nee_mean], code = 3, angle = 90, length = 0.03)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "wet"], pch = 21, bg = "#2E86AB", type = "p", cex = 1.4)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "mod"], pch = 21, bg = "#F6AE2D", type = "p", cex = 1.4)
# points(nee_mean ~ vpd_mean, dfB[REW_tercile == "dry"], pch = 21, bg = "#E15554", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "VPD (kPa)")
# mtext(side = 2, font = 4, line = 2.3, cex = 1.4, text = expression(bolditalic(paste("NEE (gC m"^"-2", "day"^"-1", ")", sep = ""))))
# legend("bottomright", bty = "n", legend = c("Wet: p < 0.001 ***", "Mid: p < 0.001 ***", "Dry: p < 0.001 ***"), fill = c("#2E86AB", "#F6AE2D", "#E15554"), text.font = 4, cex = 1.2)
# legend("topleft", bty = "n", legend = c("Variance decomposition (full data)"), text.font = 4, cex = 1)
# title(main = "Days split into equal-sized REW/VPD bins")
# 
# ### Variance partitioning (linear only)
# 
# m_full <- lm(nee ~ vpd + REW, data = dfPJ)
# m_vpd  <- lm(nee ~ vpd, data = dfPJ)
# m_rew  <- lm(nee ~ REW, data = dfPJ)
# 
# r2_full <- summary(m_full)$r.squared
# r2_vpd  <- summary(m_vpd)$r.squared
# r2_rew  <- summary(m_rew)$r.squared
# 
# unique_vpd <- r2_full - r2_rew
# unique_rew <- r2_full - r2_vpd
# shared     <- r2_vpd + r2_rew - r2_full
# 
# # Table for variance partitioning 
# tbl <- rbind(
#   c("Term",   "Linear"),
#   c("VPD",    paste0(round(unique_vpd,3)*100, "%")),
#   c("REW",    paste0(round(unique_rew,3)*100, "%")),
#   c("Shared", paste0(round(shared,3)*100, "%")),
#   c("Total",  paste0(round(r2_full,3)*100, "%"))
# )
# 
# x0 <- 0.15 ; y0 <- -8.5 ; line_h <- 0.35
# xrange <- diff(par("usr")[1:2])
# col_gap <- xrange * 0.06
# for (i in seq_len(nrow(tbl))) {
#   text(x = x0,             y = y0 + (i-1)*line_h, labels = tbl[i,1], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
#   text(x = x0 + col_gap,   y = y0 + (i-1)*line_h, labels = tbl[i,2], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
# }
# 


# -------------------------------------------------------------------------

### Other methods
# Partial correlations
#pc <- pcor(dfPJ[, .(gpp, vpd, REW)])
#pc$estimate
#pc$p.value

# Multiple regression
#fit_int <- lm(gpp ~ vpd + REW, data = dfPJ)
# summary(fit_int)

#fit_int_lowvpd <- lm(gpp ~ vpd * REW, data = dfPJ[vpd < 1])
# summary(fit_int_lowvpd)

#fit_gam <- gam(gpp ~ s(vpd, REW), data = dfPJ)
# summary(fit_gam)
# plot(fit_gam, scheme = 2)  # gives a nice 2D contour/perspective of the interaction surface
