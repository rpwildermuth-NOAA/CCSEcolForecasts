####################################################################################################
# SWORDFISH BRT TEMPORAL FORECAST TESTING
#
# For each terminal year:
#   1. Train the BRT using all observations up to the terminal year
#   2. Evaluate the model 1–5 years after the terminal year
#   3. Calculate training and forecast AUC
####################################################################################################

rm(list = ls())

# ================================================================================================
# 1. LIBRARIES
# ================================================================================================

library(lubridate)
library(pROC)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(gbm)
library(dismo)

set.seed(123)

# ================================================================================================
# 2. FILES AND OUTPUT DIRECTORY
# ================================================================================================

data_file <- "C:/Users/r.wildermuth/Documents/CEFI/SDMClimateForecasts/CCSEcolForecasts/data/swor_env_MOM6_gridded_daily.csv"
# data_file <- paste0(
#   "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/",
#   "FINAL/Gridded/1_Data/swor_env_MOM6_gridded_daily.csv"
# )

output_dir <- paste0(
  "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/",
  "FINAL/Gridded/Forecast_testing"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ================================================================================================
# 3. READ AND PREPARE DATA
# ================================================================================================

dataframe <- read.csv(data_file)

# Convert date explicitly
dataframe$date <- as.Date(dataframe$date)

if (any(is.na(dataframe$date))) {
  stop("Some dates could not be converted. Check the format of the date column.")
}

# Extract year
dataframe$year <- lubridate::year(dataframe$date)

# Keep years through 2017
dataframe <- dataframe %>%
  filter(year <= 2017)

cat("Years available:", min(dataframe$year), "-", max(dataframe$year), "\n")
cat("Number of observations:", nrow(dataframe), "\n")

# ================================================================================================
# 4. RESPONSE AND PREDICTORS
# ================================================================================================

response_variable <- "PresAbs"

gbm.x <- c(
  "bbv",
  "ild",
  "sos",
  "ssh",
  "ssu_rotate",
  "ssv_rotate",
  "tos",
  "deptho",
  "tos_sd",
  "ssh_sd",
  "deptho_sd",
  "eke",
  "moon_phase"
)

# Check that every required column exists
required_columns <- c(response_variable, "date", "lat", "lon", "year", gbm.x) # RW: changed to include lat and long, check for functionality with Nerea's code

missing_columns <- setdiff(required_columns, names(dataframe))

if (length(missing_columns) > 0) {
  stop(
    "These required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Remove records missing the response or predictors
dataframe <- dataframe %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(c(response_variable, gbm.x)),
      ~ !is.na(.)
    )
  )
# Make sure response is numeric 0/1
dataframe[[response_variable]] <- as.numeric(
  as.character(dataframe[[response_variable]])
)

if (!all(dataframe[[response_variable]] %in% c(0, 1))) {
  stop("PresAbs must contain only 0 and 1.")
}

# ================================================================================================
# 5. SAFE AUC FUNCTION
# ================================================================================================

calculate_auc <- function(observed, predicted) {
  
  good <- complete.cases(observed, predicted)
  
  observed <- observed[good]
  predicted <- predicted[good]
  
  # AUC requires both 0 and 1
  if (length(unique(observed)) < 2) {
    return(NA_real_)
  }
  
  auc_value <- pROC::auc(
    response = observed,
    predictor = predicted,
    direction = "<",
    quiet = TRUE
  )
  
  return(as.numeric(auc_value))
}

# ================================================================================================
# 6. FUNCTION TO FIT BRT
# ================================================================================================

fit_brt <- function(data,
                    predictor_names,
                    response_name,
                    learning_rate = 0.01) {
  
  # gbm.step can use column numbers more reliably
  predictor_columns <- match(predictor_names, names(data))
  response_column <- match(response_name, names(data))
  
  model <- dismo::gbm.step(
    data = data,                       # IMPORTANT: use training data, not global dataframe
    gbm.x = predictor_columns,
    gbm.y = response_column,
    family = "bernoulli",
    tree.complexity = 3,
    learning.rate = learning_rate,
    bag.fraction = 0.60,
    plot.main = FALSE,
    silent = TRUE
  )
  
  if (is.null(model)) {
    stop("gbm.step failed to produce a model.")
  }
  
  return(model)
}

# ================================================================================================
# 7. BRT PREDICTION FUNCTION
# ================================================================================================

predict_brt <- function(model, newdata) {
  
  best_trees <- model$gbm.call$best.trees
  
  predictions <- predict(
    object = model,
    newdata = newdata,
    n.trees = best_trees,
    type = "response"
  )
  
  return(as.numeric(predictions))
}

# ================================================================================================
# 8. FUNCTION FOR ONE TERMINAL YEAR
# ================================================================================================

build_brt <- function(data,
                      terminal_year,
                      predictor_names,
                      response_name,
                      learning_rate = 0.01,
                      maximum_horizon = 5) {
  
  cat("\nTerminal year:", terminal_year, "\n")
  
  # Training data includes everything through terminal year
  train <- data %>%
    filter(year <= terminal_year)
  
  # Test data contain only the following five calendar years
  test <- data %>%
    filter(
      year > terminal_year,
      year <= terminal_year + maximum_horizon
    )
  
  cat(
    "Training years:",
    min(train$year),
    "-",
    max(train$year),
    "| observations:",
    nrow(train),
    "\n"
  )
  
  if (nrow(train) == 0) {
    warning("No training observations for terminal year ", terminal_year)
    return(NULL)
  }
  
  if (length(unique(train[[response_name]])) < 2) {
    warning("Training data do not contain both presences and absences.")
    return(NULL)
  }
  
  # Fit model
  brt_model <- fit_brt(
    data = train,
    predictor_names = predictor_names,
    response_name = response_name,
    learning_rate = learning_rate
  )
  
  # Create output with consistent columns
  results <- data.frame(
    terminalYr = terminal_year,
    trainAUC = NA_real_,
    `1YrsOut` = NA_real_,
    `2YrsOut` = NA_real_,
    `3YrsOut` = NA_real_,
    `4YrsOut` = NA_real_,
    `5YrsOut` = NA_real_,
    check.names = FALSE
  )
  
  # Training AUC
  train_predictions <- predict_brt(
    model = brt_model,
    newdata = train
  )
  
  results$trainAUC <- calculate_auc(
    observed = train[[response_name]],
    predicted = train_predictions
  )
  
  # Forecast AUC for each horizon
  for (horizon in seq_len(maximum_horizon)) {
    
    forecast_year <- terminal_year + horizon
    
    forecast_data <- test %>%
      filter(year == forecast_year)
    
    output_column <- paste0(horizon, "YrsOut")
    
    if (nrow(forecast_data) == 0) {
      cat("  Horizon", horizon, ": year missing\n")
      results[[output_column]] <- NA_real_
      next
    }
    
    if (length(unique(forecast_data[[response_name]])) < 2) {
      cat("  Horizon", horizon, ": only one response class; AUC unavailable\n")
      results[[output_column]] <- NA_real_
      next
    }
    
    forecast_predictions <- predict_brt(
      model = brt_model,
      newdata = forecast_data
    )
    
    results[[output_column]] <- calculate_auc(
      observed = forecast_data[[response_name]],
      predicted = forecast_predictions
    )
    
    cat(
      "  Horizon",
      horizon,
      "| forecast year:",
      forecast_year,
      "| AUC:",
      round(results[[output_column]], 3),
      "\n"
    )
  }
  
  return(results)
}

# ================================================================================================
# 9. SELECT TERMINAL YEARS
# ================================================================================================

all_years <- sort(unique(dataframe$year))

minimum_training_years <- 7

if (length(all_years) <= minimum_training_years) {
  stop("There are not enough years for temporal forecasting.")
}

# The seventh available year is the first possible terminal year
# This produces seven available training years, including the terminal year
first_terminal_position <- minimum_training_years

# Do not use the final year as a terminal year because no future test year would exist
terminal_years <- all_years[
  first_terminal_position:(length(all_years) - 1)
]

cat(
  "\nTerminal years:",
  paste(terminal_years, collapse = ", "),
  "\n"
)

# ================================================================================================
# 10. RUN ALL FORECAST TESTS
# ================================================================================================
# RW: not reporting same AUC values as Barb's tests using functions below
all_auc_list <- vector("list", length(terminal_years))

for (j in seq_along(terminal_years)) {
  
  all_auc_list[[j]] <- tryCatch(
    {
      build_brt(
        data = dataframe,
        terminal_year = terminal_years[j],
        predictor_names = gbm.x,
        response_name = response_variable,
        learning_rate = 0.01,
        maximum_horizon = 5
      )
    },
    error = function(e) {
      
      warning(
        "Model failed for terminal year ",
        terminal_years[j],
        ": ",
        conditionMessage(e)
      )
      
      return(NULL)
    }
  )
}

allAUCs <- bind_rows(all_auc_list)

if (nrow(allAUCs) == 0) {
  stop("No BRT models completed successfully.")
}

# write.csv(
#   allAUCs,
#   file.path(output_dir, "AUC_by_terminal_year.csv"),
#   row.names = FALSE
# )

# ================================================================================================
# 11. CONVERT TO LONG FORMAT
# ================================================================================================

aucsLong <- allAUCs %>%
  pivot_longer(
    cols = c(
      "trainAUC",
      "1YrsOut",
      "2YrsOut",
      "3YrsOut",
      "4YrsOut",
      "5YrsOut"
    ),
    names_to = "horizon",
    values_to = "auc"
  )

aucsLong$horizon <- factor(
  aucsLong$horizon,
  levels = c(
    "trainAUC",
    "1YrsOut",
    "2YrsOut",
    "3YrsOut",
    "4YrsOut",
    "5YrsOut"
  )
)

# write.csv(
#   aucsLong,
#   file.path(output_dir, "AUC_long_format.csv"),
#   row.names = FALSE
# )

# ================================================================================================
# 12. CALCULATE SUMMARY STATISTICS
# ================================================================================================

aucsAgg <- aucsLong %>%
  group_by(horizon) %>%
  summarise(
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    n = sum(!is.na(auc)),
    se_auc = sd_auc / sqrt(n),
    .groups = "drop"
  )

# Replace NaN values caused by horizons with no valid observations
aucsAgg <- aucsAgg %>%
  mutate(
    mean_auc = ifelse(is.nan(mean_auc), NA_real_, mean_auc),
    sd_auc = ifelse(is.nan(sd_auc), NA_real_, sd_auc),
    se_auc = ifelse(is.nan(se_auc), NA_real_, se_auc)
  )

# write.csv(
#   aucsAgg,
#   file.path(output_dir, "AUC_summary_table.csv"),
#   row.names = FALSE
# )

print(aucsAgg)

# ================================================================================================
# 13. MEAN AUC BAR PLOT
# ================================================================================================

p_bar <- ggplot(
  aucsAgg,
  aes(x = horizon, y = mean_auc)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_auc - sd_auc,
      ymax = mean_auc + sd_auc
    ),
    width = 0.20,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, by = 0.1),
    oob = scales::rescale_none
  ) +
  labs(
    x = "Forecast horizon",
    y = "Mean AUC",
    title = "Swordfish BRT temporal forecast performance",
    subtitle = "Error bars represent ±1 standard deviation"
  ) +
  theme_bw()

print(p_bar)

# ggsave(
#   filename = file.path(output_dir, "AUC_barplot.png"),
#   plot = p_bar,
#   width = 7,
#   height = 5,
#   dpi = 300
# )
# 
# ggsave(
#   filename = file.path(output_dir, "AUC_barplot.pdf"),
#   plot = p_bar,
#   width = 7,
#   height = 5
# )

# ================================================================================================
# 14. AUC TRAJECTORY FOR EACH TERMINAL YEAR
# ================================================================================================

p_lines <- ggplot(
  aucsLong,
  aes(
    x = horizon,
    y = auc,
    group = factor(terminalYr)
  )
) +
  geom_line(alpha = 0.4, na.rm = TRUE) +
  geom_point(alpha = 0.6, na.rm = TRUE) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, by = 0.1)
  ) +
  labs(
    x = "Forecast horizon",
    y = "AUC",
    title = "Temporal forecast performance by terminal year",
    group = "Terminal year"
  ) +
  theme_bw()

print(p_lines)

# ggsave(
#   filename = file.path(output_dir, "AUC_by_terminal_year.png"),
#   plot = p_lines,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

cat("\nForecast testing completed.\n")
cat("Results saved in:\n", output_dir, "\n")

# RW checking SDMforecastr package functions ------------------------------

# Load the helper functions you will need
source("./R/buildSDM.R")
source("./R/scoreSDM.R")
source("./R/testSkillSDM.R")
source("./R/peelSDM.R")

# Define some parameters for building the SDM/s
# Note that for now the user supplies the SDM-specific parameters. Optimizing parameters for an SDM (esp. a BRT) 
# when you don't have good information on dataset size, number of variables, prevalence rate etc. is difficult
# This could be improved in future iterations of this code!
sdmType <- "brt" # Can currently be "gam" or "brt"
k <- 4 # Number of knots for a GAM. Will be ignored if not building a GAM
tc <- 3 # tree complexity for a BRT. Will be ignored if not building a BRT
lr <- 0.01 # learning rate for a BRT. Will be ignored if not building a BRT
max.trees <- 10000 # max trees for a BRT. Will be ignored if not building a BRT
gbm.x # Names of predictors in the SDM
response_variable  # Target variable for the SDM
aucCutoff <- 10 # If less observations than this cutoff within a month/season etc., AUC will not be calculated

yrsToTrain <- seq(7, 24) # Or minimum_training_years to end of all-years -1
yrsToForecast <- 5
# Optional: do we want to assess persistence skill? That is: can last year's environmental conditions
# predict this years species distributions?
includePersistence <-  FALSE # TRUE !!:RW for now, persistence calculations not working - need to create _lag1 columns in input dataframe
output <- vector(mode = "list", length = length(yrsToTrain)) # We will save results to a list (a list of lists)

# If observations don't have a "year", "month" and "quarter" columns, add them. Will need a "date" column to get them
if(!"year" %in% colnames(dataframe)) {
  dataframe$year <- year(dataframe$date)  
}
if(!"month" %in% colnames(dataframe)) {
  dataframe$month <- month(dataframe$date)  
}
if(!"quarter" %in% colnames(dataframe)) {
  dataframe$quarter <- quarter(dataframe$date)  
}
for (j in 1:length(output)) {
  output[[j]] <- peelSDM(obs = dataframe, sdmType = sdmType, varNames = gbm.x, 
                         targetName = response_variable,
                         k = k, tc = tc, lr = lr, max.trees = max.trees, noYrs = yrsToTrain[j], 
                         yrsToForecast = yrsToForecast, includePersistence = includePersistence)
}

# First reshape the output list so is easier to work with
suppressWarnings(rm(sdmSkill1, sdmSkillAll))
for(b in 1:length(output)) { 
  sdmSkill1 <- output[[b]]$sdmSkill
  if(exists("sdmSkillAll")) {
    sdmSkillAll <- rbind(sdmSkillAll, sdmSkill1)
  } else {
    sdmSkillAll <- sdmSkill1
  }
}

# Aggregate so skills are averaged across different months/seasons. "year" is the test year
sdmSkillAllAgg <- aggregate(auc ~ time + space + noYrsTrain + terminalYr + year, sdmSkillAll, FUN = mean, na.rm = TRUE)
# Add a field converting year to forecast time horizon
sdmSkillAllAgg$forecastHorizon <- sdmSkillAllAgg$year - sdmSkillAllAgg$terminalYr

allAUCs$noYrsTrain <- 7:24
aucNLO <- allAUCs %>%
  pivot_longer(
    cols = c(
      "trainAUC",
      "1YrsOut",
      "2YrsOut",
      "3YrsOut",
      "4YrsOut",
      "5YrsOut"
    ),
    names_to = "forecastHorizon",
    values_to = "auc"
  )
aucNLO <- aucNLO %>% filter(forecastHorizon != "trainAUC") %>% 
            mutate(forecastHorizon = substr(forecastHorizon, start = 1, stop = 1),
                   forecastHorizon = as.numeric(forecastHorizon),
                   auth = "NLO")
aucBM <- sdmSkillAllAgg %>% filter(time == "year") %>% 
          dplyr::select(terminalYr, noYrsTrain, forecastHorizon, auc) %>%
          mutate(auth = "BM")
aucBoth <- rbind(aucNLO, aucBM)

# range of difference between Barb and Nerea's calculations for full domain-year AUCs
aucBoth %>% group_by(terminalYr, noYrsTrain, forecastHorizon) %>%
  summarize(deltaAUC = max(auc) - min(auc)) %>% pull(deltaAUC) %>% range(na.rm = TRUE)

# effect of training years (sample size)
aucBoth %>% ggplot() +
  geom_point(aes(x = noYrsTrain, y = auc, colour = auth)) + facet_wrap(~forecastHorizon)
aucBoth %>% ggplot() +
  geom_boxplot(aes(x = as.factor(noYrsTrain), y = auc, colour = auth))
# effect of  forecast horizon
aucBoth %>% ggplot(aes(x = as.character(forecastHorizon), y = auc, fill = auth)) +
  geom_boxplot()
# spatiotemporal resolution X training years
sdmSkillAllAgg %>% ggplot() +
  geom_boxplot(aes(x = as.factor(noYrsTrain), y = auc)) +
  facet_grid(rows = vars(space), cols = vars(time))
# spatiotemporal resolution X forecast horizon
sdmSkillAllAgg %>% ggplot() +
  geom_boxplot(aes(x = as.character(forecastHorizon), y = auc)) +
  facet_grid(rows = vars(space), cols = vars(time))

# find best scale for forecasts
sdmSkillAllAgg %>% group_by(space, time, forecastHorizon) %>%
  summarize(maxAUC = max(auc, na.rm = TRUE),
           medAUC = median(auc, na.rm = TRUE),
           meanAUC = mean(auc, na.rm = TRUE),
           minAUC = min(auc, na.rm = TRUE),
           yrN = n()) %>%
  # arrange(desc(medAUC)) %>% filter(medAUC >= 0.5) %>% print(n=35)
  arrange(desc(medAUC), desc(minAUC)) %>% filter(minAUC >= 0.45) # 
  