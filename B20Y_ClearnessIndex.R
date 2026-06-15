library(progress)
library(data.table)
dt <- fread("data/ICOSFLUXNET_allsites_allyears.csv")
dt <- dt[site == "FR-Fon"]

# Deal with duplicate data
dtu <- dt[, if (.N == 1) .SD, by = .(timestamp, site)]
dtd <- dt[, if (.N > 1) .SD, by = .(timestamp, site)]
dtd <- dtd[file == "icos"]
dt <- rbind(dtu, dtd)
setorder(dt, site, timestamp)
rm(dtu,dtd)

dst_offset <- function(date) { # Clear sky data uses local time but I fed it UTC so I need to add +1 or +2 hours
  # DST in Europe: last Sunday of March to last Sunday of October
  ifelse(format(date, "%m") %in% c("04","05","06","07","08","09") |
           (format(date, "%m") == "03" & as.integer(format(date, "%d")) >= 25) |
           (format(date, "%m") == "10" & as.integer(format(date, "%d")) < 25),
         2L, 1L)
}
dt[, date := as.Date(paste(y, m, d, sep = "-"))]
dt[, dst := dst_offset(date)]


### Find the clear PAR (done once) ------------------------------------------

# Import source
# source("C:/Users/duran/Dropbox/CNRS/R/active/repos/2025-NikyOvercast/NOV_MAESPA_Utils.R")
# 
# dfC <- clearpar(localtime = "00:00:00", lat = 48.476357, lon = -2.780096, tz = 0, date = "2005-01-01", quantum = TRUE)[0,]
# dfC[1:368160,] <- NA
# pb <- progress_bar$new(format = "(:spin) [:bar] :percent | ETA: :eta", total = nrow(dt))
# 
# for(i in 1:nrow(dt))
# {
#   dfi <- dt[i]
#   idate <- paste(dfi[1,y], ifelse(dfi[1,m] < 10, paste0("0", dfi[1,m]), dfi[1,m]), ifelse(dfi[1,d] < 10, paste0("0", dfi[1,d]), dfi[1,d]), sep = "-")
#   itime <- paste(ifelse(dfi[1,h] < 10, paste0("0", dfi[1,h]), dfi[1,h]), ifelse(dfi[1,min] < 10, paste0("0", dfi[1,min]), dfi[1,min]), "00", sep = ":")
#   itz <- dfi[1,dst]
# 
#   dfC[i,] <- clearpar(localtime = itime, lat = 48.476357, lon = -2.780096, tz = itz, date = idate, quantum = FALSE)
#   pb$tick()
# }
# fwrite(dfC, "Barbeau20years_clearsky.csv")


# Correct local time ------------------------------------------------------
dst_offset <- function(date) { # Clear sky data uses local time but I fed it UTC so I need to add +1 or +2 hours
  # DST in Europe: last Sunday of March to last Sunday of October
  ifelse(format(date, "%m") %in% c("04","05","06","07","08","09") |
           (format(date, "%m") == "03" & as.integer(format(date, "%d")) >= 25) |
           (format(date, "%m") == "10" & as.integer(format(date, "%d")) < 25),
         7200L, 3600L)
}

dfC <- fread("Barbeau20years_clearsky.csv")

dfC$h <- as.numeric(substring(dfC$localtime, 1, 2))
dfC$h <- ifelse(dfC$h < 10, paste0("0", dfC$h), dfC$h)
dfC$min <- as.numeric(substring(dfC$localtime, 4, 5))
dfC$min <- ifelse(dfC$min < 10, paste0("0", dfC$min), dfC$min)
dfC$posix <- as.POSIXct(paste0(dfC$date, " ", dfC$h, ":", dfC$min, ":00"), format = "%Y-%m-%d %H:%M:%S")
# dfC$posix <- dfC$posix + dst_offset(dfC$date)
dfC$posix <- dfC$posix + 3600

dfC$y <- as.numeric(substring(dfC$posix, 1, 4))
dfC$m <- as.numeric(substring(dfC$posix, 6, 7))
dfC$d <- as.numeric(substring(dfC$posix, 9, 10))
dfC$h <- as.numeric(substring(dfC$posix, 12, 13))
dfC$h <- ifelse(dfC$h < 10, paste0("0", dfC$h), dfC$h)
dfC$min <- as.numeric(substring(dfC$posix, 15, 16))
dfC$min <- ifelse(dfC$min < 10, paste0("0", dfC$min), dfC$min)
dfC$date <- as.Date(paste(dfC$y, dfC$m, dfC$d, sep = "-"))
dfC$timestamp <- paste(gsub("-", "", dfC$date), dfC$h, dfC$min, sep = "")

dfC$RG <- dfC$RG / 2.1 # Also need to correct RG units

# Merge data --------------------------------------------------------------
dt[,timestamp := as.character(timestamp)]
dt[,h := NULL] ; dt[,min := NULL] ; dt[,y := NULL] ; dt[,m := NULL] ; dt[,d := NULL] ; dt[,date := NULL] ;   # Remove doublon columns
dtX <- merge(dt, dfC, by = "timestamp")
setDT(dtX)
dtX[, time := as.numeric(h) + (as.numeric(min)/60)]

# Scaling by year ---------------------------------------------------------

# Get the closest days between Mod and Mes
# Those are clear days and calculate difference
# Apply it as a yearly scaling factor

plotDay <- function(day = sample(unique(dtX$date), 1), data = dtX)
{
  dfi <- data[date == day]
  plot(-500, xlim = c(0,24), ylim = c(0,1000), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  points(RG ~ time, dfi, pch = 21, bg = "gray", type = "o")
  points(rtot ~ time, dfi, pch = 21, bg = "red", type = "o")
  axis(side = 1, font = 4)
  axis(side = 2, font = 4, las = 2)
  mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Time (h.)")
  mtext(side = 2, font = 4, line = 3, cex = 1.4, text = "PPFD")
  title(main = paste(as.character(day)))
  legend("topleft", bty = "n", legend = c("Clear sky", "Measured"), fill = c("gray", "red"), cex = 1, text.font = 2)
}

# Calculate error
dtX[, error := abs(RG - rtot)]
dtE <- aggregate(dtX[,.(error, RG, rtot)], by = list("date" = dtX$date), sum, na.rm = T)
setDT(dtE)
dtE[, y := substring(date, 1, 4)]
dtE[, percError := error/RG]

# Look at the nine closest per year
iyear <- 2005
for(iyear in unique(dtE$y))
{
  dty <- dtE[y == iyear]
  dty <- dty[order(percError)]
  dates <- head(dty$date, n = 9)
  
  par(mfrow = c(3,3))
  for(i in seq_along(dates))
  {
    idate <- dates[i]
    plotDay(day = idate, data = dtX)
  }
}

# Clearness index ---------------------------------------------------------
dtE[, kt := rtot / RG]
dtE[, rtot_mean := rtot / 48]
# dtE[kt > 1, kt := 1]
# dtE[, snowsuspect := month(date) %in% c(10, 11, 12, 1, 2, 3, 4) & kt > 1.05]
# dtE <- dtE[snowsuspect == FALSE]

dtZ <- dtE[, .(frac_clear = mean(kt > 0.65)), by = year(date)]

dtY <- aggregate(dtE[,.(kt, rtot_mean)], by = list("y" = dtE$y), mean, na.rm = T)
setDT(dtY)
dtY[, y := as.numeric(y)]
dtY <- dtY[!y == 2024]

dtY25 <- aggregate(dtE[,.(kt, rtot_mean)], by = list("y" = dtE$y), quantile, probs = 0.25, na.rm = T)
dtY50 <- aggregate(dtE[,.(kt, rtot_mean)], by = list("y" = dtE$y), quantile, probs = 0.50, na.rm = T)
dtY75 <- aggregate(dtE[,.(kt, rtot_mean)], by = list("y" = dtE$y), quantile, probs = 0.75, na.rm = T)
setDT(dtY25) ; setDT(dtY50) ; setDT(dtY75)
dtY25[, y := as.numeric(y)]
dtY50[, y := as.numeric(y)]
dtY75[, y := as.numeric(y)]
dtY25 <- dtY25[!y == 2024]
dtY50 <- dtY50[!y == 2024]
dtY75 <- dtY75[!y == 2024]

mK25 <- lm(kt ~ y, dtY25)
mK50 <- lm(kt ~ y, dtY50)
mK75 <- lm(kt ~ y, dtY75)
summary(mK25)
summary(mK50)
summary(mK75)

mK <- lm(kt ~ y, dtY)
mR <- lm(rtot_mean ~ y, dtY)
mF <- lm(frac_clear ~ year, dtZ)
mRK <- lm(rtot_mean ~ kt, dtY)

summary(mR)
summary(mK)
summary(mF)
summary(mRK)

plot(kt ~ rtot_mean, dtY)

par(mfrow = c(1,1), bty = "L")
plot(-500, xlim = c(2005,2025), ylim = c(0.6,0.9), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(kt ~ y, dtY, pch = 21, bg = "gray", type = "p", cex = 1.2)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
mtext(side = 2, font = 4, line = 3, cex = 1.4, text = "Clearness Index")
abline(a = coef(mK)[1], b = coef(mK)[2], col = "red")
legend("topleft", bty = "n", legend = c("With 2024:  +2.5% clearness in 20 years (ns)", "Without 2024: +5.3% clearness in 20 years (*)"), text.font = 4)
#abline(a = coef(mK50)[1], b = coef(mK50)[2], col = "orangered")

plot(-500, xlim = c(2005,2025), ylim = c(120,160), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(rtot_mean ~ y, dtY, pch = 21, bg = "gray", type = "p", cex = 1.2)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
mtext(side = 2, font = 4, line = 3, cex = 1.4, text = "Global radiation")
abline(a = coef(mR)[1], b = coef(mR)[2], col = "red")


plot(-500, xlim = c(2005,2025), ylim = c(0.45,0.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(frac_clear ~ year, dtZ, pch = 21, bg = "gray", type = "p", cex = 1.2)
axis(side = 1, font = 4)
axis(side = 2, font = 4, las = 2)
mtext(side = 1, font = 4, line = 2.5, cex = 1.4, text = "Years")
mtext(side = 2, font = 4, line = 3, cex = 1.4, text = "% clear days")
abline(a = coef(mF)[1], b = coef(mF)[2], col = "red")

predict(mR)
predict(mK)

library(quantreg)
fit_q <- rq(kt ~ y, tau = c(0.25, 0.5, 0.75, 0.9), data = dtY)
summary(fit_q)

# old ---------------------------------------------------------------------


