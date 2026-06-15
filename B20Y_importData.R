library(data.table)
library(readxl)

# EC data
path <- "data/EC/"
files <- list.files(path, full.names = TRUE)[grep("LI7500", list.files(path, full.names = TRUE))]
files <- files[grep("hh", files)]

df0 <- data.table()
for(i in 1:length(files))
{
  ifile <- files[i]
  df <- read_excel(ifile, sheet = "data")
  df0 <- rbind(df0, df)
}
fwrite(df0, "data/B20Y_EC_all.csv")

# BM data
path <- "data/BM/"
files <- list.files(path, full.names = TRUE)[grep(".csv", list.files(path, full.names = TRUE))]

i = 1
df0 <- data.table()
for(i in 1:length(files))
{
  ifile <- files[i]
  print(ifile)
  header <- names(fread(ifile, nrows = 0))
  df <- fread(ifile, skip = 2, header = FALSE, col.names = header)
  df <- df[!TIMESTAMP == ""]
  df[df == -9999] <- NA
  df0 <- rbind(df0, df, fill = T)
}
fwrite(df0, "data/B20Y_BM_all.csv")
