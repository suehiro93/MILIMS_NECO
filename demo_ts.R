# Read R files
# NOTE: If you don't use "CPLEX" optimizer with "Rcplex", please set "RCPLEX <- FALSE" in "setup.R"
source("setup.R")

# Read datasets:
italy_train_orign <- read_ts("./data4demo/ItalyPowerDemand/ItalyPowerDemand_TRAIN")
italy_test_orign <- read_ts("./data4demo/ItalyPowerDemand/ItalyPowerDemand_TEST")

# Set parameters:
param.list <- set_MILIMS_parameter() # use default parameters

# Set shapelet lengths
shapelet_lengths <- round(c(0.3, 0.4)*ncol(italy_train_orign$x))

# Run MILIMS for time-series classification and the shapelet-based classification model
model <- MILIMS4TS_script(italy_train_orign$x, italy_train_orign$y, shapelet_lengths, param.list)

# Calculate classification accuracy for test set
calc_accuracy4TS(italy_model,italy_test_orign$x,italy_test_orign$y)

# Get predicted labels
val <- evalFun_TS(model, italy_test_orign$x)
predicted_labels <- sign(val)