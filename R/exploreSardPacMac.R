# Code to explore P. mackerel/P. sardine incidental catch product
# Created 5/16/2025, Robert Wildermuth

library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(cubelyr)

# modified from code in ProcessData.Rmd in recrmntDFA
datPath <- "C:/Users/r.wildermuth/Documents/CEFI/SDMClimateForecasts/SDMoutput/"

# sardine spawning season
cpsFiles <- expand.grid(1:12, 1998:2023, "SDMs.nc") %>%
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

cpsFiles <- paste0(datPath, "HistoricalCPSGAMsBRTsERDDAPROMSDomain/", cpsFiles)
cpsSDMs <- read_stars(cpsFiles, proxy = FALSE, quiet = TRUE)
# Time fix from Barb
wrongTimes <- as.Date(st_get_dimension_values(cpsSDMs, "time"))
wrongTimes[c(1, length(wrongTimes))] # Can see is wrong baseline
rightTimes <- wrongTimes + 693579 # 693579 is the number of days since 0000-00-00
cpsSDMs <- st_set_dimensions(cpsSDMs, 3, values = rightTimes, names = "time")
st_get_dimension_values(cpsSDMs, "time")[c(1, length(wrongTimes))] # Check times look ok now

st_crs(cpsSDMs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
cpsSDMs <- st_transform(cpsSDMs, "+proj=longlat +ellps=WGS84 +datum=WGS84")
#st_get_dimension_values(cpsSDMs, "time")

sardMackSDMs <- cpsSDMs %>% select(sardGAM, chubGAM, jackGAM)


# Use CalCurrent Atlantis extent to limit sample frame
ccAtl <- st_read("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/emocc_whole_domain.shp")
ccAtl <- st_transform(ccAtl, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# # restrict to area south of 40deg N
# newAtlbbox <- st_bbox(ccAtl)
# newAtlbbox$ymax <- 40.0
# newAtlbbox <- st_bbox(unlist(newAtlbbox))
# ccAtl <- st_crop(ccAtl, newAtlbbox)

# crop to Atlantis extent
sardMackSDMs <- st_crop(x = sardMackSDMs,
                    y = ccAtl)

# show first 5 days of SDM
sardSlice <- slice(sardMackSDMs, index = 1:5, along = "time")
ggplot() + geom_stars(data = sardSlice, color = NA) + facet_wrap(~time)

# discrimination threshold from Barb Muhling:
# sardHabThresh <- 0.45 # this is for no-SSB, shape-constrained sardine GAM
# values for ERDDAP GAMs from Barb, 4/16/2025
sardThrGAM <- 0.43
pmackThrGAM <- 0.18 
jackThrGAM <- 0.27  

# Try masking low probability cells
# sardSpawn <- sardSDMs
# sardSpawn[sardSDMs < sardHabThresh] <- NA
# Don't allow spawning habitat above 40deg north
#sardSpawn2 <- sardSpawn %>% filter(y < 40)

# For now assume occurrences are independent
sardMackSDMs <- sardMackSDMs %>% mutate(coocrSarChb = sardGAM * chubGAM,
                                        coocrSarJck = sardGAM * jackGAM)


# Map prep
pac.coast <- borders("world", colour="brown", fill="brown", xlim = c(-140, -100), ylim = c(20, 60))
mycols <- RColorBrewer::brewer.pal(9, "Greens")#colors()[c(473,562,71,610,655,653,621,34)]
mypalette <- colorRampPalette(mycols)(255)

# sardWk <- aggregate(sardMackSDMs, by = "weeks", FUN = sum, na.rm = TRUE)
# # reorder
# sardWk <- aperm(sardWk, c(2,3,1))
# 
# ggplot() +
#   # geom_stars(data = sardWk, aes(fill = coocrSarChb), color = NA) +
#   geom_stars(data = sardWk, aes(fill = coocrSarJck), color = NA) +
#   # scale_fill_gradientn(colours = mypalette, limits = c(0, max(sardWk$coocrSarChb,
#   scale_fill_gradientn(colours = mypalette, limits = c(0, max(sardWk$coocrSarJck,
#                                                          na.rm = TRUE)), 
#                        na.value = "transparent") +
#   guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
#   pac.coast + 
#   geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
#   coord_sf(xlim = c(-130, -113), ylim = c(27, 50)) +
#   facet_wrap(~time)

sardMackMon <- aggregate(sardMackSDMs, by = "months", FUN = sum, na.rm = TRUE)
st_get_dimension_values(sardMackMon, "time")

# reorder and then make 0s NAs to turn them grey
sardMackMon <- aperm(sardMackMon, c(2,3,1))
sardMackMon[sardMackMon == 0] <- NA

ggplot() +
  geom_stars(data = sardMackMon, aes(fill = coocrSarChb), color = NA) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  scale_fill_gradientn(colours = mycols[3:9],
                       limits = c(0, max(sardMackMon$coocrSarChb,
                                         na.rm = TRUE)),
                       na.value = "transparent") +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  coord_sf(xlim = c(-130, -115), ylim = c(28, 50)) +
  facet_wrap(~time, nrow = 3) +
  theme_minimal()

# plot multiple attributes with facet_grid (from: https://rpubs.com/michaeldorman/646276)
coocrDims <- st_redimension(sardMackMon %>% 
                              select(sardGAM, chubGAM, coocrSarChb)) %>% 
                slice(index = 11, along = "time")
names(coocrDims) <- "value"

ggplot() +
  geom_stars(data = coocrDims, color = NA) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  scale_fill_gradientn(colours = mycols,
                       limits = c(0, max(coocrDims$value,
                                         na.rm = TRUE)),
                       na.value = "transparent") +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  coord_sf(xlim = c(-130, -115), ylim = c(28, 50)) +
  facet_wrap(~new_dim) +
  theme_minimal()

test1 <- sardMackMon %>% select(coocrSarChb) %>%
            aggregate(by = ccAtl, FUN = sum, na.rm = TRUE)
# 
test2 <- data.frame(time = st_get_dimension_values(test1, "time"),
                      coocrSarChb = t(test1$coocrSarChb))

test2 %>% 
  ggplot() +
  geom_line(aes(x = time, y = coocrSarChb))
