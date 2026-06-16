library(data.table)

dfE <- fread("data/B20Y_EC_all.csv")
dfE[dfE == -9999] <- NA
dfE <- dfE[, c(1:8, 10:19, 26, 30, 36, 38, 39, 43, 46:51, 60)]
colnames(dfE) <- c("y", "m", "d", "h", "doy", "par", "rg", "ta", "rh", "vpd", "press", "p", "ws", "ustar", "fc", "sc", "co2", "nee", "gpp", "ter", "le", "h2o", "e", "he", "rn", "swin", "swout", "lwin", "lwout", "gmoy", "etp")

# Constants
lam     <- 2.45e6      # latent heat of vaporisation (J/kg)
rho     <- 1000        # density of water (kg/m3)
cpa     <- 1013        # specific heat of air (J/kg/K)
eps     <- 0.622       # ratio molecular weight water/dry air

# Conversions
dfE[, press := press / 10] # in kPa
dfE[, del := 4098 * (0.6108 * exp(17.27 * ta / (ta + 237.3))) / (ta + 237.3)^2] # Slope of saturation vapour pressure curve (before conversion to celsius)
dfE[, ta := ta + 273.15]   # in Kelvins
dfE[, gam := (cpa * press) / (eps * lam)]
dfE[, ws := ws * 4.87 / log(67.8 * 37 - 5.42)]

# Soil heat flux 
dfE[, g  := gmoy * 1800 / 1e6]     # in MJ/m2 per half hour [1800 seconds in 30 min]
dfE[, rn := rn * 1800 / 1e6]       # in MJ/m2 per half hour [1800 seconds in 30 min]
dfE[is.na(g), g := 0.05 * rn]   # rn now in MJ, so this is consistent

# Calculate ETP
dfE[, ETP_mm := (0.408 * del * (rn - g) + gam * (37 / ta) * ws * vpd) / (del + gam * (1 + 0.34 * ws))] # in mm/half-hour
dfE[, ETP := ETP_mm / 1800 * lam]

dfE[ETP < 0, ETP := 0]
dfE[ETP_mm < 0, ETP_mm := 0]
dfE[etp < 0, etp := 0]

# plot(ETP ~ etp, dfE, xlim = c(0,1200), ylim = c(0,1200), pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.1))
# abline(0, 1, col = "red")
# summary(lm(ETP ~ etp, dfE))

# ET / PET ratio
dfE[, etpet := e / ETP_mm]

dfE[ETP_mm <= 0.001, etpet := NA]   # avoid division by near-zero
dfE[rn <= 0, etpet := NA]           # nighttime meaningless
dfE[etpet < 0, etpet := NA]        # negative ET is not meaningful here

# Aggregate to daily
dfEJ <- dfE[, .(
  etpet   = mean(etpet, na.rm = TRUE),
  gpp     = mean(gpp, na.rm = TRUE),
  e       = sum(e, na.rm = TRUE),
  etp     = sum(ETP_mm, na.rm = TRUE),
  ta      = mean(ta - 273.15, na.rm = TRUE),
  vpd     = mean(vpd, na.rm = TRUE),
  rn      = sum(rn, na.rm = TRUE)
), by = .(y, m, d, doy)]

dfEJ[, etpet_sum := e / etp] # Daily sum more accurate perhaps

# Define growing season by GPP
GPP_max <- dfEJ[, .(gppmax = max(gpp, na.rm = TRUE)), by = y]
dfEJ <- merge(dfEJ, GPP_max, by = "y")
dfEJ[, in_gs := gpp > 0.1 * gppmax]

dfGS <- dfEJ[, .(y, m, d, doy, in_gs)]
# fwrite(dfGS, "data/df_inGS.csv")
rm(dfGS, GPP_max)

dfEY <- dfEJ[in_gs == TRUE, .(
  e       = sum(e, na.rm = TRUE),
  etp     = sum(etp, na.rm = TRUE),
  ta      = mean(ta, na.rm = TRUE),
  vpd     = mean(vpd, na.rm = TRUE),
  rn      = sum(rn, na.rm = TRUE)
), by = y]
dfEY[, etpet := e / etp]

# plot(e ~ y, dfEY, cex = 1.4, bg = "gray", pch = 21)
# plot(etp ~ y, dfEY, cex = 1.4, bg = "gray", pch = 21)
# plot(etpet ~ y, dfEY, cex = 1.4, bg = "gray", pch = 21)

mE     <- lm(e ~ y, dfEY)
mETP   <- lm(etp ~ y, dfEY)
mETPET <- lm(etpet ~ y, dfEY)
summary(mE)
summary(mETP)
summary(mETPET)
rm(mE, mETP, mETPET)

# Cleanup
rm(cpa, eps, lam, rho)
# -------------------------------------------------------------------------

# png("figs/B20Y_annualETPET.png", height = 6000, width = 15000, res = 1000)
# par(mfrow = c(1,3), bty = "L", mar = c(5,4,1,1), oma = c(0,1,0,0))
# 
# plot(-500, xlim = c(2004,2026), ylim = c(450,650), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(e ~ y, dfY, pch = 21, bg = "dodgerblue3", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("Actual evapotranspiration [ET] (mm year"^"-1", ")", sep = ""))))
# abline(a = coef(mE)[1], b = coef(mE)[2], col = "dodgerblue1", lwd = 3, lty = 3)
# legend("topleft", bty = "n", legend = c("p = 0.60 | R2 = 0"), text.font = 4, cex = 1.4)
# 
# plot(-500, xlim = c(2004,2026), ylim = c(650,1050), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(etp ~ y, dfY, pch = 21, bg = "dodgerblue3", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("Potential evapotranspiration [PET] (mm year"^"-1", ")", sep = ""))))
# abline(a = coef(mETP)[1], b = coef(mETP)[2], col = "dodgerblue1", lwd = 3, lty = 1)
# legend("topleft", bty = "n", legend = c("p = 0.0009 | R2 = 0.41"), text.font = 4, cex = 1.4)
# 
# plot(-500, xlim = c(2004,2026), ylim = c(0.5,0.85), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# points(etpet ~ y, dfY, pch = 21, bg = "dodgerblue3", type = "p", cex = 1.4)
# axis(side = 1, font = 4)
# axis(side = 2, font = 4, las = 2)
# mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
# mtext(side = 2, font = 4, line = 2.8, cex = 1.4, text = expression(bolditalic(paste("ET / PET (ratio actual to potential)", sep = ""))))
# abline(a = coef(mETPET)[1], b = coef(mETPET)[2], col = "dodgerblue1", lwd = 3, lty = 1)
# legend("topleft", bty = "n", legend = c("p = 0.004 | R2 = 0.32"), text.font = 4, cex = 1.4)
# 
# dev.off()

