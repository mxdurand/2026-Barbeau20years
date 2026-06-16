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
# dfPJ <- dfPJ[y >= 2018]   # Check some subset years


# VPD / SWC are largely independent (important for partial correlation)
# cor(dfPJ[, .(vpd, REW)]) # -0.088


# Binned response curves
dfPJ[, REW_tercile := cut(REW, breaks = quantile(REW, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          labels = c("dry", "mod", "wet"), include.lowest = TRUE)]
dfPJ[, vpd_bin := cut(vpd, breaks = quantile(vpd, probs = seq(0, 1, 0.1), na.rm = TRUE),
                      include.lowest = TRUE)]

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

# GPP ---------------------------------------------------------------------

### Variance partitioning

# Linear
m_full <- lm(gpp ~ vpd + REW, data = dfPJ)
m_vpd  <- lm(gpp ~ vpd, data = dfPJ)
m_rew  <- lm(gpp ~ REW, data = dfPJ)

r2_full <- summary(m_full)$r.squared
r2_vpd  <- summary(m_vpd)$r.squared
r2_rew  <- summary(m_rew)$r.squared

unique_vpd <- r2_full - r2_rew
unique_rew <- r2_full - r2_vpd
shared     <- r2_vpd + r2_rew - r2_full

# Non linear (GAM)
m_full2 <- gam(gpp ~ s(vpd) + s(REW) + ti(vpd, REW), data = dfPJ)
m_vpd2 <- gam(gpp ~ s(vpd), data = dfPJ) # VPD only
m_rew2 <- gam(gpp ~ s(REW), data = dfPJ) # REW only
r2_full2 <- summary(m_full2)$r.sq
r2_vpd2 <- summary(m_vpd2)$r.sq
r2_rew2 <- summary(m_rew2)$r.sq

unique_vpd2 <- r2_full2 - r2_rew2   # what VPD adds beyond REW alone... 
unique_rew2 <- r2_full2 - r2_vpd2   # what REW adds beyond VPD alone
shared2     <- r2_vpd2 + r2_rew2 - r2_full2

# Table for variance partitioning 
# tbl <- rbind(
#   c("Term",   "Linear", "GAM"),
#   c("VPD",    paste0(round(unique_vpd,3)*100, "%"),  paste0(round(unique_vpd2,3)*100, "%")),
#   c("REW",    paste0(round(unique_rew,3)*100, "%"),  paste0(round(unique_rew2,3)*100, "%")),
#   c("Shared", paste0(round(shared,3)*100, "%"),  paste0(round(shared2,3)*100, "%")),
#   c("Total",  paste0(round(r2_full,3)*100, "%"),  paste0(round(r2_full2,3)*100, "%"))
# )
tbl <- rbind(
  c("Term",   "Linear"),
  c("VPD",    paste0(round(unique_vpd,3)*100, "%")),
  c("REW",    paste0(round(unique_rew,3)*100, "%")),
  c("Shared", paste0(round(shared,3)*100, "%")),
  c("Total",  paste0(round(r2_full,3)*100, "%"))
)


### Other methods
# Partial correlations
pc <- pcor(dfPJ[, .(gpp, vpd, REW)])
pc$estimate
pc$p.value

# Multiple regression
fit_int <- lm(gpp ~ vpd + REW, data = dfPJ)
# summary(fit_int)

fit_int_lowvpd <- lm(gpp ~ vpd * REW, data = dfPJ[vpd < 1])
# summary(fit_int_lowvpd)

fit_gam <- gam(gpp ~ s(vpd, REW), data = dfPJ)
# summary(fit_gam)
# plot(fit_gam, scheme = 2)  # gives a nice 2D contour/perspective of the interaction surface


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

# Plot
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
title(main = "Days split into equal-sized REW/VPD bins")

x0 <- 0.1 ; y0 <- 16.7 ; line_h <- 0.5
xrange <- diff(par("usr")[1:2])
col_gap <- xrange * 0.06
for (i in seq_len(nrow(tbl))) {
  text(x = x0,             y = y0 - (i-1)*line_h, labels = tbl[i,1], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
  text(x = x0 + col_gap,   y = y0 - (i-1)*line_h, labels = tbl[i,2], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
}


# ET ----------------------------------------------------------------------

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
title(main = "Days split into equal-sized REW/VPD bins")

### Variance partitioning (linear only)

m_full <- lm(e ~ vpd + REW, data = dfPJ)
m_vpd  <- lm(e ~ vpd, data = dfPJ)
m_rew  <- lm(e ~ REW, data = dfPJ)

r2_full <- summary(m_full)$r.squared
r2_vpd  <- summary(m_vpd)$r.squared
r2_rew  <- summary(m_rew)$r.squared

unique_vpd <- r2_full - r2_rew
unique_rew <- r2_full - r2_vpd
shared     <- r2_vpd + r2_rew - r2_full

# Table for variance partitioning 
tbl <- rbind(
  c("Term",   "Linear"),
  c("VPD",    paste0(round(unique_vpd,3)*100, "%")),
  c("REW",    paste0(round(unique_rew,3)*100, "%")),
  c("Shared", paste0(round(shared,3)*100, "%")),
  c("Total",  paste0(round(r2_full,3)*100, "%"))
)

x0 <- 0.15 ; y0 <- 5.5 ; line_h <- 0.2
xrange <- diff(par("usr")[1:2])
col_gap <- xrange * 0.06
for (i in seq_len(nrow(tbl))) {
  text(x = x0,             y = y0 - (i-1)*line_h, labels = tbl[i,1], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
  text(x = x0 + col_gap,   y = y0 - (i-1)*line_h, labels = tbl[i,2], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
}


# TER ---------------------------------------------------------------------

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
title(main = "Days split into equal-sized REW/VPD bins")

### Variance partitioning (linear only)

m_full <- lm(ter ~ vpd + REW, data = dfPJ)
m_vpd  <- lm(ter ~ vpd, data = dfPJ)
m_rew  <- lm(ter ~ REW, data = dfPJ)

r2_full <- summary(m_full)$r.squared
r2_vpd  <- summary(m_vpd)$r.squared
r2_rew  <- summary(m_rew)$r.squared

unique_vpd <- r2_full - r2_rew
unique_rew <- r2_full - r2_vpd
shared     <- r2_vpd + r2_rew - r2_full

# Table for variance partitioning 
tbl <- rbind(
  c("Term",   "Linear"),
  c("VPD",    paste0(round(unique_vpd,3)*100, "%")),
  c("REW",    paste0(round(unique_rew,3)*100, "%")),
  c("Shared", paste0(round(shared,3)*100, "%")),
  c("Total",  paste0(round(r2_full,3)*100, "%"))
)

x0 <- 0.15 ; y0 <- 9.2 ; line_h <- 0.35
xrange <- diff(par("usr")[1:2])
col_gap <- xrange * 0.06
for (i in seq_len(nrow(tbl))) {
  text(x = x0,             y = y0 - (i-1)*line_h, labels = tbl[i,1], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
  text(x = x0 + col_gap,   y = y0 - (i-1)*line_h, labels = tbl[i,2], adj = 0, cex = 0.65, font = ifelse(i==1, 2, 1))
}

dev.off()

# NEE ---------------------------------------------------------------------

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
