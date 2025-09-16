rm(list = ls())

library(ncdf4)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# Load trained BRT model
brt_model <- readRDS("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/native/SWOR.res1.tc3.lr03.single_native.rds")
out_dir="C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\4_predict/native/"
# Set your .nc path and variables
nc_files <- list.files("C:/Users/nereo/Dropbox (Personal)/Nerea/NOAA/PROJECTS & COLLABORATIONS/PROJECTS/Forecast/Swordfish_forecast/6_dynamic_native/", 
                       pattern = "\\.nc$", full.names = TRUE)
var_names <- c("deptho", "sos", "ssh", "ssu", "ssv", "tos")

# Define raster template
r_template <- rast(extent = c(-180, 180, -90, 90), resolution = 0.25, crs = "EPSG:4326")

# Get number of time steps from first variable
test_nc <- nc_open(nc_files[1])
nt <- dim(ncvar_get(test_nc, var_names[1]))[3]


nc <- nc_open(nc_files[2])
time_raw <- ncvar_get(nc, "time")  # may be in hours/days since a start date
time_units <- ncatt_get(nc, "time", "units")$value
nc_close(nc)

# Convert time to Date
# Example: "days since 1950-01-01"
origin_str <- sub(".*since ", "", time_units)
dates <- as.Date(time_raw, origin = origin_str)


nc_close(test_nc)

# Load coastline
world <- ne_countries(scale = "medium", returnclass = "sf")

# Custom palette
pal <- colorRampPalette(c("#9b59b6", "#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))
ncolors <- 100
custom_colors <- pal(ncolors)

n_days <- nlayers(dates)
dim(dates)

# Loop over each day
for (i in 1:dates) {  # fixed loop to go over length of dates
  cat("Processing day", i, "\n")
  
  raster_list <- list()
  
  date_str <- format(dates[i], "%Y-%m-%d")
  
  
  for (i in seq_along(nc_files)) {
    nc <- nc_open(nc_files[i])
    varname <- var_names[i]
    var_dims <- nc$var[[varname]]$dim
    dim_names <- sapply(var_dims, function(x) x$name)
    
    # Identify lon/lat
    if ("xh" %in% dim_names) lon <- ncvar_get(nc, "xh")
    if ("xq" %in% dim_names) lon <- ncvar_get(nc, "xq")
    if ("yh" %in% dim_names) lat <- ncvar_get(nc, "yh")
    if ("yq" %in% dim_names) lat <- ncvar_get(nc, "yq")
    
    lon_mod <- ifelse(lon > 180, lon - 360, lon)
    
    # Extract time slice
    if (length(var_dims) == 2) {
      var_array <- ncvar_get(nc, varname)
    } else if (length(var_dims) == 3) {
      var_array <- ncvar_get(nc, varname, start = c(1, 1, t), count = c(-1, -1, 1))
    } else {
      stop(paste("Unsupported dimensions in", varname))
    }
    
    nc_close(nc)
    
    lon_exp <- rep(lon_mod, times = length(lat))
    lat_exp <- rep(lat, each = length(lon_mod))
    var_vec <- as.vector(var_array)
    
    points_sf <- st_as_sf(data.frame(lon = lon_exp, lat = lat_exp, var = var_vec),
                          coords = c("lon", "lat"), crs = 4326)
    
    points_vect <- vect(points_sf)
    r_var <- rasterize(points_vect, r_template, field = "var", fun = mean, background = NA)
    names(r_var) <- varname
    raster_list[[i]] <- r_var
  }
  
  r_stack <- rast(raster_list)
  names(r_stack) <- c("depth", "sal", "ssh", "ssu", "ssv", "sst")
  predict_df <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)
  predict_df_clean <- na.omit(predict_df)
  
  # Predict
  predictions <- predict(brt_model, newdata = predict_df_clean[, -c(1,2)], type = "response")
  
  # Convert to spatial
  pred_df <- data.frame(x = predict_df_clean$x,
                        y = predict_df_clean$y,
                        prediction = predictions)
  pred_sf <- st_as_sf(pred_df, coords = c("x", "y"), crs = 4326)
  pred_vect <- vect(pred_sf)
  r_pred <- rasterize(pred_vect, r_template, field = "prediction")
  
 
  
  r_pred_df <- as.data.frame(r_pred, xy = TRUE, na.rm = TRUE)
  
  
  # Plot with white background
  p <- ggplot() +
    geom_raster(data = r_pred_df, aes(x = x, y = y, fill = last)) +
    geom_sf(data = world, fill = "gray80", color = "black", size = 0.3) +
    scale_fill_gradientn(name = "Prediction", colors = custom_colors, limits = c(0, 1),na.value = "transparent") +
    coord_sf(xlim = c(-130, -117), ylim = c(30, 47), expand = FALSE) +
    labs(title = paste("BRT Prediction - Day", date_str),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      panel.background = element_rect(fill = "white", color = "white"),  # Set plot background to white
      plot.background = element_rect(fill = "white", color = "white"),   
    )
  
  png_file <- file.path(out_dir, paste0("prediction_", date_str, ".png"))
  ggsave(paste0("prediction_day_", date_str, ".png"), plot = p, width = 8, height = 6, dpi = 300)
 
   # Save raster as .grd file
  grd_file <- file.path(out_dir, paste0("prediction_", date_str, ".grd"))
  writeRaster(r_pred, filename = paste0("prediction_day_", date_str, ".grd"), overwrite = TRUE)
}

