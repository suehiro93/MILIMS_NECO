make_sliding_matrix <- function(n, Q){
  ret.mat <- matrix(0, Q*n, n)
  j <- 1
  for(i in 1:n){
    ret.mat[j:(j+Q-1), i] <- 1
    j <- j+Q
  }
  return(ret.mat)
}

make_sliding_matrix2 <- function(n, nrows){
  ret.mat <- matrix(0, sum(nrows), n)
  j <- 1
  for(i in 1:n){
    ret.mat[j:(j+nrows[i]-1), i] <- 1
    j <- j+nrows[i]
  }
  return(ret.mat)
}

make_sliding_matrix2 <- cmpfun(make_sliding_matrix2)

kernelBC <- function(kerneldot, exam.list, shape.mat){
  m <- length(exam.list)
  kernel.list <- vector("list", m)
  for(i in 1:m){
    kernel.list[[i]] <- kernelMatrix(kerneldot, exam.list[[i]], shape.mat)
  }
  return(kernel.list)
}

kernelBC <- cmpfun(kernelBC)

get_segment <- function(series, L, R){
  r <- L-R+1
  ret.mat <- matrix(NA,r,R)
  for(i in 1:r){
    ret.mat[i, ] <- series[i:(i+R-1)]
  }
  ret.mat <- na.omit(ret.mat)
  return(ret.mat)
}
get_segment <- cmpfun(get_segment)


gen_all_of_subsequences <- function(T.mat, ell){
  L <- ncol(T.mat)
  m <- nrow(T.mat)
  subseq.list <- vector("list", m)
  for(j in 1:m){
    series <- T.mat[j,]
    subseq.list[[j]] <- get_segment(series, L, ell)
  }
  return(subseq.list)
}
gen_all_of_subsequences <- cmpfun(gen_all_of_subsequences)

matlist2mat <- function(matlist){
  Q <- nrow(matlist[[1]])
  listnum <- length(matlist)
  rownum <- Q*listnum
  colnum <- ncol(matlist[[1]])
  ret.mat <- matrix(NA, nrow=rownum, ncol=colnum)
  for(i in 1:listnum){
    ret.mat[((i-1)*Q+1):(i*Q),] <- matlist[[i]]
  }
  return(ret.mat)
}

matlist2mat2 <- function(matlist, row_nums){
  total_rownum <- sum(row_nums)
  colnum <- ncol(matlist[[1]])
  ret.mat <- matrix(NA, nrow=total_rownum, ncol=colnum)
  row_id <- 1
  for(i in 1:length(matlist)){
    ret.mat[row_id:(row_id+row_nums[i]-1),] <- matlist[[i]]
    row_id <- row_id + row_nums[i]
  }
  return(ret.mat)
}

matlist2mat2 <- cmpfun(matlist2mat2)

set_MILIMS_parameter <- function(sigma=0.01, nu=0.1, Km=100, boostiter=500, kerneldot=rbfdot, optimizer=1){
  param.list <- list()
  param.list$kerneldot <- kerneldot
  param.list$sigma <- sigma
  param.list$boostiter <- boostiter
  param.list$Km <- Km
  param.list$nu <- nu
  param.list$optimizer <- optimizer
  return(param.list)
}