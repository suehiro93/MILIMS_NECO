# Shapelet-based Multiple-Insatnce Learning  
* This repository gives the implementation of the method published in [1]. This method generally works for Multiple-Instance Learning tasks and Shapelet-Learning tasks for time-series classification.
* This is the first shepelet-learning method with the theoretical generalization performance.
* Competitively performs with SOTA methods in practice.

1. Daiki Suehiro, Kohei Hatano, Eiji Takimoto, Shuji Yamamoto, Kenichi Bannai, Akiko Takeda, "Theory and Algorithms for Shapelet-based Multiple-Instance Learning", Neural Computation, to appear (arXiv version is ).

"MILIMS" comes from "Multiple-Instance Learning by Infinitely Many Shapelet-based classifiers".

---

# Requirement
* R (version 3.4 or later.)
* R Packages: "data.table", "kernlab", "qlcMatrix", "R.matlab", "MASS", "compiler", and "Rcplex" or "lpSolveAPI". For visualization, you need "zoo", "RColorBrewer", "gplots"  
* Software  
CPLEX (A Solver provided by IBM. If you use Rcplex).  

To know how to install CPLEX and Rcplex, the following pages may be helpful for you:
* http://www-stat.wharton.upenn.edu/~josezubi/INSTALL
* https://cran.r-project.org/web/packages/Rcplex/INSTALL
* https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/

---

# Usage

* Install required packages:  
`install.packages(c("data.table", "kernlab", "qlcMatrix", "R.matlab", "MASS"))`
* To solve LP problems, one of the following packages is required:  
`install.package("Rcplex") # recommended, but you need to install CPLEX software (free for academic user), and you need some setting for configuring Rcplex.`  
`install.package("lpSolveAPI") # You can intanstly use our method without installing external software and setting, but less efficient than Rcplex.`  

* Read R files for setting up:  
***NOTE: If you don't use "CPLEX" optimizer with "Rcplex", please set `"RCPLEX <- FALSE"` in "setup.R"***  
`source("setup.R")`

* Run demos:  
Shapelet-Learning tasks for time-series classification:  
`source("./demo_ts.R")`  
MIL tasks:  
`source("./demo_mi.R")`  
  
Below is the simple commands for time-series classification. For MIL tasks, see demo_mi.R

* Read datasets:  
`italy_train_origin <- read_ts("./data4demo/ItalyPowerDemand/ItalyPowerDemand_TRAIN")`  
`italy_test_origin <- read_ts("./data4demo/ItalyPowerDemand/ItalyPowerDemand_TEST")`  

* Set parameters:  
`param.list <- set_MILIMS_parameter() # use default parameters`  

* Set shapelet lengths:  
`shapelet_lengths <- round(c(0.3, 0.4)*ncol(italy_train_origin$x))`  

* Run MILIMS for time-series classification and the shapelet-based classification model:  
`model <- MILIMS4TS_script(italy_train_origin$x, italy_train_origin$y, shapelet_lengths, param.list)`  

* Calculate classification accuracy for test set:  
`calc_accuracy4TS(italy_model,italy_test_origin$x,italy_test_origin$y)`  

* Get predicted labels:  
`val <- evalFun_TS(model, italy_test_origin$x)`  
`predicted_labels <- sign(val)`  


# Inputs of main functions
MILIMS requires:  
* Training bags: **list**.  
Each bag (i.e., the element of the list) contains set of instances as a matrix (rownum: #sample, colnum: #dimension)
* The labels of the bags: **vector (+1 or -1)**.  
*Note: We are now implementing the code for multi-class classification. The current version is available for binary classification.*  
* The numbers of instances: **vector of integers**.   
* Parameters: **list** (details are later).  
* seed (default=1): **small integer**. We can change initial seed of kmeans function.  
* SHAPELET (default=FALSE): **boolean**. If TRUE, the algorithms roughly solve the weak learning problem (see Appendix B.1 in the paper). This mode may be useful for the rough hyper-parameter search.  

MILIMS4TS_script requires:  
* Time-series datasets: **matrix** (rownum: #samples, colnum: #time-series length).  
If you want to use time-series datasets that contain various lengths of time-series, 
please pad the shorter time series with NA.*  
* The labels of the time series: **vector (+1 or -1)**.  
*Note: We are now implementing the code for multi-class classification. The current version is available for binary classification.*  
* Parameters: **list** (details are later).  
* Lengths of shapelets (hyper-parameter): **vector of integers**.  
* seed (default=1): **small integer**. We can change initial seed of kmeans function.  
* SHAPELET (default=FALSE): **boolean**. If TRUE, the algorithms roughly solve the weak learning problem (see Appendix B.1 in the paper). This mode may be useful for the rough hyper-parameter search.  


# Outputs (learned model)
The classification model is a convex combination of shapelet-based classifiers. The fomula is g(B) in the end of Section 4 in the paper.  
* w: weights of shapelet-based classifiers.  
* b: bias term of the classification function.  
* alpha, Kx: \sum \alpha K(x,) (=u in the paper). large alpha implies the importance of x for classification.  
* kerneldot: kernel function.  
* sigma: parameter of kernel.  
* train_time: running time.  
* ells (only for MILIMS4TS): shapelet lengths.  
* ells_ids (only for MILIMS4TS): index of ells for each shapelet-based classifier.  

# About hyper-parameters
We can set the hyperparameters as a list using "set_MILIMS_parameter".  
The important parameters are:
* nu: Lower bound of the training error. That is, if nu=0.2, our the training error of the output model does not exceed 0.2. Note that very small nu induces the overfitting. If the margin (rho) is 0, the output model should be worthless because of the small nu.
* sigma: The parameter of the Gaussian kernel. Note that we use the Gaussian kernel implemented by kernlab package, i.e., K(x,z) = exp(-sigma||x-z||^2).
* Km: K of the k-means clustering algorithm (see Section 6.4 in the paper). If Km=NULL, the algorithm use all instances appearing in the all training bags. Our algorithm runs faster with smaller K, but very small K decrease the classification performance. If you have large number of instances in total, please set Km as some constant ingteger (e.g., 100). For time-series classificaiton, the total number of instances (i.e., the number of subsequences) is very large.  

If you cannot obtain the desired classification accuracy, we recommend to tune the above hyperparameters.  
  
If you run out of memory, please set the parameter "optimizer" as 2. The algorithm solves the weak learning problems by AdaGrad (we need to set a learning rate ETA_adagrad (defalut =10)).  

# Visualization
We give a sample code for the visualization of obtained shapelet-based classifiers.  

`source("./visualization_beta_ver.R")`

Get a model  
`model <- MILIMS4TS_script(italy_train_origin$x, italy_train_origin$y, shapelet_lengths, param.list)`  

Observe top-5 (shapelet-like) important short sequences (warm colored sequence contributes to positive, cold colored negative) in the classifier.  
`visualize_shapelets(italy_test_origin$x[1,], model, 5)`  

Observe maximizers in the bag (i.e., matched sequences corresponds to the shapelets)
for the top-10 shapelets in the classifier.  
`visualize_shapelets2(italy_test_origin$x[1,,drop=FALSE], model, 10, y_min=-3, y_max=3)`

# Acknowledgement  
We use some of the following datasets for the demos.
* https://www.cs.ucr.edu/~eamonn/time_series_data/
* http://www.cs.columbia.edu/~andrews/mil/datasets.html
