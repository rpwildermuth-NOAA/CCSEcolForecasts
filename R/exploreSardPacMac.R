# Code to explore P. mackerel/P. sardine incidental catch product
# Created 5/16/2025, Robert Wildermuth

library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(cubelyr)
library(nngeo)
library(pROC)
library(dismo)


# modified from code in ProcessData.Rmd in recrmntDFA
datPath <- "C:/Users/r.wildermuth/Documents/CEFI/SDMClimateForecasts/"

# sardine spawning season
cpsFiles <- expand.grid(1:12, 1998:2023, "SDMs.nc") %>%
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

cpsFiles <- paste0(datPath, "SDMoutput/HistoricalCPSGAMsBRTsERDDAPROMSDomain/", cpsFiles)
cpsSDMs <- read_stars(cpsFiles, proxy = FALSE, quiet = TRUE, along = "time")
# Time fix from Barb
wrongTimes <- as.Date(st_get_dimension_values(cpsSDMs, "time"))
wrongTimes[c(1, length(wrongTimes))] # Can see is wrong baseline
# rightTimes <- wrongTimes + 693579 # 693579 is the number of days since 0000-00-00
rightTimes <- wrongTimes + 10226 # 35794 is the number of days since 1970-01-01
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
  geom_stars(data = sardMackMon[,,,1:6], aes(fill = coocrSarChb), color = NA) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  scale_fill_gradientn(colours = mycols[3:9],
                       limits = c(0, max(sardMackMon$coocrSarChb,
                                         na.rm = TRUE)),
                       na.value = "transparent") +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  coord_sf(xlim = c(-130, -115), ylim = c(28, 50)) +
  facet_wrap(~time, nrow = 2) +
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

# look at co-occurrence observations from survey --------------------------

trawlDat <- read_csv(file = paste0(datPath, "CPS_Trawl_LifeHistory_Haulcatch_dl20250804.csv"))

# find hauls with sardine and PacMac together
sardHauls <- trawlDat %>% filter(scientific_name == "Sardinops sagax")
pacmacHauls <- trawlDat %>% filter(scientific_name == "Scomber japonicus")
coocrObs <- inner_join(x = sardHauls, y = pacmacHauls, 
                       by = c("cruise", "ship", "haul", "collection", 
                              "start_latitude", "start_longitude", "stop_latitude",
                              "stop_longitude", "equilibrium_time", "haulback_time"))
coocrObs <- coocrObs %>% select("cruise", "ship", "haul", "collection", 
                                "start_latitude", "start_longitude", "stop_latitude",
                                "stop_longitude", "equilibrium_time", "haulback_time") %>%
              left_join(y = trawlDat, by = c("cruise", "ship", "haul", "collection", 
                                             "start_latitude", "start_longitude", 
                                             "stop_latitude", "stop_longitude", 
                                             "equilibrium_time", "haulback_time"))

# get presence/absence of sardine in pacmac hauls
presabsSardinMac <- left_join(x = pacmacHauls, y = sardHauls,  
                               by = c("cruise", "ship", "haul", "collection", 
                                      "start_latitude", "start_longitude", "stop_latitude",
                                      "stop_longitude", "equilibrium_time", "haulback_time"))
presabsSardinMac <- presabsSardinMac %>% select("cruise", "ship", "haul", "collection", 
                                                "start_latitude", "start_longitude", "stop_latitude",
                                                "stop_longitude", "equilibrium_time", "haulback_time",
                                                scientific_name.x, presence_only.x,
                                                scientific_name.y, presence_only.y) %>%
                      mutate(pa = as.numeric(!is.na(scientific_name.y)))
presabsSardinMac %>% summarize(macAndSard = sum(pa),
                               macNoSard = sum(pa == 0))
317/nrow(presabsSardinMac)
# ~70% of hauls with P. mackerel also have P. sardine

presabsMacinSard <- left_join(x = sardHauls, y = pacmacHauls,  
                              by = c("cruise", "ship", "haul", "collection", 
                                     "start_latitude", "start_longitude", "stop_latitude",
                                     "stop_longitude", "equilibrium_time", "haulback_time"))
presabsMacinSard <- presabsMacinSard %>% select("cruise", "ship", "haul", "collection", 
                                                "start_latitude", "start_longitude", "stop_latitude",
                                                "stop_longitude", "equilibrium_time", "haulback_time",
                                                scientific_name.x, presence_only.x,
                                                scientific_name.y, presence_only.y) %>%
                      mutate(pa = as.numeric(!is.na(scientific_name.y)))
presabsMacinSard %>% summarize(sardAndMac = sum(pa),
                               sardNoMac = sum(pa == 0))
317/nrow(presabsMacinSard)

# see if %sard by weight can be calculated
# max incidental take of sardine in other directed CPS fisheries is 20% of landed weight when sardine are overfished
# see sec4.5.1 Rebuilding Plan for Pacific Sardine, https://www.pcouncil.org/documents/2023/06/coastal-pelagic-species-fishery-management-plan.pdf/
# assume "landed weight" is weight of CPS finfish in each haul - avoids dealing with presence_only = Y
coocrObs <- coocrObs %>% mutate(remaining_weight = case_when(is.na(remaining_weight) ~ 0,
                                                             TRUE ~ remaining_weight), 
                                totWeight = subsample_weight + remaining_weight)
landedPct <- coocrObs %>% filter(scientific_name %in% c("Sardinops sagax", 
                                                        "Scomber japonicus", 
                                                        "Trachurus symmetricus", 
                                                        "Engraulis mordax")) %>%
                group_by(cruise, haul, collection) %>%
                summarize(landedWeight = sum(totWeight))
landedPct <- landedPct %>% filter(!is.na(landedWeight)) %>% # remove when subsample weight not available
                left_join(coocrObs %>% filter(scientific_name %in% c("Sardinops sagax", 
                                                                     "Scomber japonicus", 
                                                                     "Trachurus symmetricus", 
                                                                     "Engraulis mordax")),
                          by = c("cruise", "haul", "collection")) %>%
                mutate(pctLandWt = totWeight/landedWeight*100) 
landedPct %>%
  select(haul, collection, scientific_name, totWeight, landedWeight, pctLandWt)

# number of co-occurring hauls
coocrObs %>% filter(scientific_name == "Sardinops sagax") %>% nrow()
# number of co-occuring hauls with sardine weight above limit
landedPct %>% filter(scientific_name == "Sardinops sagax",
                     pctLandWt >= 20) %>% nrow()
155/317
# ~49% of hauls with both P. mackerel and P. sardine are above the percentage landing limit

# look at spread of co-occurrences in space and time
coocrObs %>% ggplot(aes(x = stop_longitude, y = stop_latitude)) + 
  geom_point()
coocrObs %>% mutate(Mo = month(equilibrium_time),
                    Yr = year(equilibrium_time)) %>%
  distinct(cruise, haul, collection, stop_latitude, stop_longitude, Mo, Yr) %>%
  summary()
coocrObs %>% mutate(Mo = month(equilibrium_time),
                    Yr = year(equilibrium_time)) %>%
  distinct(cruise, haul, collection, stop_latitude, stop_longitude, Mo, Yr) %>%
  select(Mo, Yr) %>% table()

landedPct %>% filter(scientific_name == "Sardinops sagax",
                     pctLandWt >= 20) %>% 
  ggplot(aes(x = stop_longitude, y = stop_latitude)) + 
  geom_point()
landedPct %>% filter(scientific_name == "Sardinops sagax",
                     pctLandWt >= 20) %>%
  mutate(Mo = month(equilibrium_time),
         Yr = year(equilibrium_time)) %>%
  summary()

# Figure out how to overlap predictions in sardMackSDMs with presence/absence in presabsSardinMac
presabsSardinMac <- presabsSardinMac %>% mutate(time = date(equilibrium_time))

paSardinMac <- presabsSardinMac %>% st_as_sf(coords=c("stop_longitude","stop_latitude"),
                                             crs="+proj=longlat +ellps=WGS84 +datum=WGS84",
                                             remove=FALSE)
# probably need to use extraction function Stephanie and Owen figured out
# Try Owen's point:nearest neighbor method
sfSardMackSDMs <- sardMackSDMs %>% select(coocrSarChb) %>% st_as_sf(as_points = TRUE)
# only need dates with pres/abs data available
obsDates <- as.character(unique(presabsSardinMac$time))
obsDates <- obsDates[-grep("2024", obsDates)]
sfSardMackSDMs <- sfSardMackSDMs %>% select(all_of(obsDates))
# spatial join using nearest neighbors
# loop over observation days
paPreds <- tibble(time = '', 
                  pa = 0,
                  pred = 0)[0,]
for(dd in 1:length(obsDates)){
  trawl_nn <- st_join(paSardinMac, sfSardMackSDMs[, dd], st_nn, k=1)
  trawl_nn <- trawl_nn %>% select(time, pa, all_of(obsDates[dd])) %>% 
                filter(time == obsDates[dd]) %>%
                rename(pred = obsDates[dd])
  paPreds <- rbind(paPreds, trawl_nn)
}

# calculate point-wise AUC
ptAUC <- auc(paPreds$pa, paPreds$pred, direction = "<", quiet = TRUE)
# RW: not great performance :(

# do like Nerea does it
pres <- paPreds %>% filter(pa == 1) %>% pull(pred)
abs <- paPreds %>% filter(pa == 0) %>% pull(pred)
evalPtAUC <- evaluate(p=pres, a=abs)
evalPtAUC
