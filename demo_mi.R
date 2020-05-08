# Read R files
# If you don't use "CPLEX" optimizer with "Rcplex", please set "RCPLEX <- FALSE" in "setup.R"
source("setup.R")

# Read datasets:
musk1norm <- make_bags.list4MI_DATA("./data4demo/musk1norm_matlab.mat")
sample_num <- length(musk1norm$bag_labels)
train_ids <- sample(1:sample_num, round(sample_num*0.7))

# Set parameters:
param.list <- set_MILIMS_parameter() # use default parameters

# Run MILIMS for MI classification and the shapelet-based classification model
model <- MILIMS(musk1norm$bags.list[train_ids], musk1norm$bag_labels[train_ids], musk1norm$instance_nums[train_ids],param.list)

# Calculate classification accuracy for test set
calc_accuracy(model, musk1norm$bags.list[-train_ids], musk1norm$bag_labels[-train_ids])

# Get predicted labels
val <- evalFun(model, musk1norm$bags.list[-train_ids])
predicted_labels <- sign(val)