read_ts <- function(filename){
  library(data.table)
  tmp.mat <- as.matrix(fread(filename,header=FALSE))
  x <- tmp.mat[,-1]
  y <- tmp.mat[,1]
  uni_y <- sort(unique(y))
  if(length(uni_y)==2){
    bin <- c(-1,1)
    for(i in 1:length(uni_y)){
      y[y==uni_y[i]] <- bin[i]
    }
  }
  ret.list <- list(x=x, y=y)
  return(ret.list)
}

ts_data2MDmi_data <- function(ts.dat, ts.labels, ells){
  ell_num <- length(ells)
  bags.list.list <- vector("list", length=ell_num)
  instance_nums.list <- vector("list", length=ell_num)
  for(i in 1:ell_num){
    bags.dat <- ts_data2mi_data(ts.dat, ts.labels, ells[i])
    bags.list.list[[i]] <- bags.dat$bags.list
    instance_nums.list[[i]] <- bags.dat$instance_nums
  }
  return(list(bags.list.list=bags.list.list, bag_labels = ts.labels, instance_nums.list=instance_nums.list))
}

ts_data2mi_data <- function(ts_data, labels, ell, scaling=FALSE){
  ret.list <- gen_all_of_subsequences(ts_data, ell)
  m <- nrow(ts_data)
  instance_nums <- rep(0, m)
  for(i in 1:m){
    instance_nums[i] <- nrow(ret.list[[i]])
  }
  if(scaling == TRUE){
    for(i in 1:m){
      for(j in 1:instance_nums[i]){
        tmp <- znorm(ret.list[[i]][j,])
        ret.list[[i]][j,] <- tmp
      }
    }
  }
  return(list(bags.list = ret.list, bag_labels = labels, instance_nums = instance_nums, instance_labels = NULL))
}

make_bags.list4MI_DATA <- function(filename, SCALE=TRUE, NORM=FALSE){
  if(!(regexpr('\\.csv$', filename) < 0)){
    return(gen_mi_data4ticaret(filename, NORM))
  }
  dat <- readMat(filename)
  bag_ids <- as.numeric(dat$bag.ids)
  labels <- as.numeric(dat$labels)
  if(SCALE==FALSE){
    features <- as.data.frame(as.matrix(dat$features))
  }else{
    features <- scale(as.matrix(dat$features))
    features[is.nan(features)] <- 0
    features <- as.data.frame((features))
  }
  if(NORM==TRUE){
    features <- normalize_01(features)
  }
  dat.df <- cbind(bag_ids, labels, features)
  sp_dat.df <- split(dat.df, dat.df$bag_ids)
  bag_num <- length(sp_dat.df)
  bag_labels <- rep(-1, bag_num)
  ret.list <- vector("list", bag_num)
  instance_nums <- rep(0, bag_num)
  for(i in 1:bag_num){
    ret.list[[i]] <- as.matrix(sp_dat.df[[i]][,c(-1,-2)])
    instance_nums[i] <- nrow(sp_dat.df[[i]])
    if(is.element(1, sp_dat.df[[i]]$labels)){
      bag_labels[i] <- 1
    }
  }
  return(list(bags.list = ret.list, bag_labels = bag_labels, instance_nums = instance_nums, instance_labels = labels))
}

gen_mi_data4ticaret <- function(csv_file, NORM=FALSE){
  dat <- as.matrix(read.csv(csv_file, header = FALSE))
  bag_ids <- as.numeric(dat[,2])
  labels <- as.numeric(dat[,1])
  features <- as.data.frame(as.matrix(dat[,c(-1,-2)]))
  if(NORM==TRUE){
    features <- normalize_L2ball(features)
  }
  dat.df <- cbind(bag_ids, labels, features)
  sp_dat.df <- split(dat.df, dat.df$bag_ids)
  bag_num <- length(sp_dat.df)
  bag_labels <- rep(-1, bag_num)
  ret.list <- vector("list", bag_num)
  instance_nums <- rep(0, bag_num)
  for(i in 1:bag_num){
    browser()
    ret.list[[i]] <- as.matrix(sp_dat.df[[i]][,c(-1,-2)])
    instance_nums[i] <- nrow(sp_dat.df[[i]])
    if(is.element(1, sp_dat.df[[i]]$labels)){
      bag_labels[i] <- 1
    }
  }
  return(list(bags.list = ret.list, bag_labels = bag_labels, instance_nums = instance_nums, instance_labels = labels))
}

normalize_01 <- function(dat.mat){
  max_val <- max(dat.mat)
  min_val <- min(dat.mat)
  normalized_dat.mat <- dat.mat - min_val
  temp <- max_val - min_val
  normalized_dat.mat <- normalized_dat.mat/temp
  return(normalized_dat.mat)
}

normalize_L2ball <- function(dat.mat){
  max_norm <- 0
  for(i in 1:nrow(dat.mat)){
    dat.norm <- sqrt(sum((as.numeric(dat.mat[i,]))^2))
    if(dat.norm > max_norm){
      max_norm <- dat.norm
    }
  }
  return(dat.mat/max_norm)
}