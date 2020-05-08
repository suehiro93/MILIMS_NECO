## Maximum number of restart for our algorithms to avoid the fail induced by the randomness of k-means
MAX_RESTART <- 30
## Maximum number of DCA iteration for weak learn
MAX_WLT <- 50
## Maximum number of the iteration for AdaGrad
MAX_GRAD_ITER <- 100 
## Initial learning rate of AdaGrad
ETA_adagrad <- 10
shape_num <- NULL
KM_SUBSEQ_LIM <- Inf


MILIMS4TS_script <- function(ts.dat, ts.labels, ells, param.list, seed=1, SHAPELET=FALSE){
  ell_num <- length(ells)
  bags.list.list <- vector("list", length=ell_num)
  instance_nums.list <- vector("list", length=ell_num)
  for(i in 1:ell_num){
    bags.dat <- ts_data2mi_data(ts.dat, ts.labels, ells[i])
    bags.list.list[[i]] <- bags.dat$bags.list
    instance_nums.list[[i]] <- bags.dat$instance_nums
  }
  ret <- MILIMS4TS(bags.list.list, ts.labels, instance_nums.list, param.list, ells, seed=seed, SHAPELET=SHAPELET)
  return(ret)
}

MILIMS4TS <- function(bags.list.list, bag_labels, instance_nums.list, 
                      param.list, ells, seed=1, SHAPELET=FALSE, restart_num=MAX_RESTART){
  time1 <- proc.time()
  ord_y <- order(bag_labels)
  y <- matrix(bag_labels[ord_y], nrow=1)
  optimizer <- param.list$optimizer
  if(is.null(optimizer)){
    optimizer <- 1
  }
  T <- param.list$boostiter
  threshold <- 1e-8*1
  nu <- param.list$nu
  Km <- param.list$Km
  m <- length(y)
  ypos_bool <- y==1
  yneg_bool <- y==-1
  ell_num <- length(ells)
  Kall.list.list <- vector("list", length=ell_num)
  Kpos.list.list <- vector("list", length=ell_num)
  Kneg.list.list <- vector("list", length=ell_num)
  pos_nums <- rep(NA, ell_num)
  neg_nums <- rep(NA, ell_num)
  max_vals.list.list <- vector("list", length=ell_num)
  neg_instance_nums.list <- vector("list", length=ell_num)
  shape_nums <- rep(NA, ell_num)
  shape.mat.list <- vector("list", length=ell_num)
  kerneldot <- param.list$kerneldot
  kerneldot_original <- kerneldot
  for(i in 1:ell_num){
    message("preprocessing:", i)
    bags.list <- bags.list.list[[i]]
    bags.list <- bags.list[ord_y]
    instance_nums <- instance_nums.list[[i]]
    instance_nums <- instance_nums[ord_y]
    pos_bags.list <- bags.list[ypos_bool]
    neg_bags.list <- bags.list[yneg_bool]
    ell <- ells[i]
    kerneldot <- kerneldot_original(param.list$sigma/ell)
    pos_shape.mat <- matlist2mat2(pos_bags.list, instance_nums[ypos_bool])
    neg_shape.mat <- matlist2mat2(neg_bags.list, instance_nums[yneg_bool])
    if(!is.null(Km)){
      set.seed(seed)
      if(Km!=0){
        if(nrow(pos_shape.mat)>Km){
          if(nrow(pos_shape.mat)>KM_SUBSEQ_LIM){
            sample_id <- sample(1:nrow(pos_shape.mat), KM_SUBSEQ_LIM)
            pos_shape.mat <- pos_shape.mat[sample_id,]
          }
          pos_shape.mat <- kmeans(pos_shape.mat, Km, iter.max = 100)$centers
        }
        if(nrow(neg_shape.mat)>Km){
          if(nrow(neg_shape.mat)>KM_SUBSEQ_LIM){
            sample_id <- sample(1:nrow(neg_shape.mat), KM_SUBSEQ_LIM)
            neg_shape.mat <- neg_shape.mat[sample_id,]
          }          
          neg_shape.mat <- kmeans(neg_shape.mat, Km, iter.max = 100)$centers
        }
      }
    }
    shape.mat <- rbind(pos_shape.mat, neg_shape.mat)
    shape_num <<- nrow(shape.mat)
    Kall.list <- kernelBC(kerneldot, bags.list, shape.mat)
    pos_num <- length(pos_bags.list)
    neg_num <- length(neg_bags.list)
    for(j in 1:m){
      Kall.list[[j]] <- cbind(Kall.list[[j]], -Kall.list[[j]])
    }
    Kpos.list <- Kall.list[ypos_bool]
    Kneg.list <- Kall.list[yneg_bool]
    max_vals.list <- get_closed_vals(Kpos.list, Kneg.list, shape_num)
    Kall.list.list[[i]] <- Kall.list
    Kpos.list.list[[i]] <- Kpos.list
    Kneg.list.list[[i]] <- Kneg.list
    pos_nums[i] <- pos_num
    neg_nums[i] <- neg_num
    max_vals.list.list[[i]] <- max_vals.list
    neg_instance_nums.list[[i]] <- instance_nums[yneg_bool]
    shape_nums[i] <- shape_num
    shape.mat.list[[i]] <- shape.mat
    gc(reset = TRUE)
    gc(reset = TRUE)
  }
  ### run LPBoost ###
  d <- rep(1/m, m)
  w <- rep(1, T) ## weight of weak classifiers
  alpha.list <- vector("list", T)
  dpos <- d[ypos_bool]
  dneg <- d[yneg_bool]
  gamma <- 0
  Kall_vals.mat <- matrix(NA, m, T)
  VAL_bef <- Inf  
  ell_ids <- rep(0, T)
  message("Preprocessing completed")
  for(iter in 1:T){
    message("#Iteration =", iter)
    ### run DCA
    sol <- MI_WeakLearn_MD(neg_instance_nums.list, Kpos.list.list, Kneg.list.list, 
                           dpos, dneg, max_vals.list.list, ell_num, shape_nums, threshold, SHAPELET, optimizer)
    alpha <- sol$alpha
    ell_ids[iter] <- sol$ell_id
    d_nonzero_id <- which(d!=0)
    Kall_vals <- rep(NA, m)
    alpha <- as.matrix(alpha, nrow=1)
    shape_num <- length(alpha)/2
    Kall.list <- Kall.list.list[[sol$ell_id]]
    for(i in 1:m){
      Kall_vals[i] <- max(tcrossprod(alpha,Kall.list[[i]]))
    }
    alpha <- matrix(alpha[1:shape_num] - alpha[(shape_num+1):(shape_num*2)], nrow=1)
    alpha.list[[iter]] <- alpha
    Kall_vals.mat[, iter] <- Kall_vals
    VAL <- crossprod(t(y),(d*Kall_vals))
    if(abs(VAL) < gamma){
      alpha.list <- alpha.list[1:(iter-1)]
      Kall_vals.mat <- Kall_vals.mat[ ,1:(iter-1)]
      ell_ids <- ell_ids[1:(iter-1)]
      break
    }
    if(VAL_bef == VAL){
      alpha.list <- alpha.list[1:(iter-1)]
      Kall_vals.mat <- Kall_vals.mat[ ,1:(iter-1)]
      ell_ids <- ell_ids[1:(iter-1)]
      break      
    }else{
      VAL_bef <- VAL
    }
    solution_LP <- run_LPBoost_dual(t(Kall_vals.mat[, 1:iter]), y, nu=nu*m)
    d <- solution_LP$d
    dpos <- d[ypos_bool]
    dneg <- d[yneg_bool]
    gamma <- solution_LP$gamma
  }
  solution_LP <- run_1normSVM(Kall_vals.mat, y, nu=nu*m)
  w <- solution_LP$w
  b <- solution_LP$b
  rho <- solution_LP$rho  
  train_time <- (proc.time()-time1)[3]
  message("rho=",rho)
  if(rho<0){
    restart_num <- restart_num - 1
    if(restart_num==0){
      message("The obtained solution is not meaningful.")
      return(NULL)
    }
    message("Restart with changing seed, maybe because the kmeans clustering failed.")
    return(MILIMS4TS(bags.list.list, bag_labels, instance_nums.list, 
                     param.list, ells, seed=seed*(restart_num+1), SHAPELET=SHAPELET, restart_num=restart_num))
  }
  return(list(alpha=alpha.list, Kx.list=shape.mat.list, w=w, b=b, rho=rho, kerneldot=kerneldot_original, sigma=param.list$sigma, ells = ells, ell_ids=ell_ids, train_time=train_time))
}


MILIMS <- function(bags.list, bag_labels, instance_nums, param.list, seed=1, SHAPELET=FALSE, restart_num=MAX_RESTART){
  original_bags.list <- bags.list
  original_instance_nums <- instance_nums
  time1 <- proc.time()
  ord_y <- order(bag_labels)
  bags.list <- bags.list[ord_y]
  instance_nums <- instance_nums[ord_y]
  y <- matrix(bag_labels[ord_y], nrow=1)
  T <- param.list$boostiter
  threshold <- 1e-8*1
  nu <- param.list$nu
  kerneldot <- param.list$kerneldot(param.list$sigma)
  Km <- param.list$Km
  m <- length(y)
  ypos_bool <- y==1
  yneg_bool <- y==-1
  pos_bags.list <- bags.list[ypos_bool]
  neg_bags.list <- bags.list[yneg_bool]
  ell <- ncol(bags.list[[1]])
  pos_shape.mat <- matlist2mat2(pos_bags.list, instance_nums[ypos_bool])
  neg_shape.mat <- matlist2mat2(neg_bags.list, instance_nums[yneg_bool])
  if(!is.null(Km)){
    set.seed(seed)
    if(Km!=0){
      if(nrow(pos_shape.mat)>Km){
        pos_shape.mat <- kmeans(pos_shape.mat, Km, nstart=30, iter.max = 100)$centers
      }
      if(nrow(neg_shape.mat)>Km){
        neg_shape.mat <- kmeans(neg_shape.mat, Km, nstart=30, iter.max = 100)$centers
      }
    }
  }
  shape.mat <- rbind(pos_shape.mat, neg_shape.mat)
  shape_num <<- nrow(shape.mat)
  Kall.list <- kernelBC(kerneldot, bags.list, shape.mat)
  Kpos.list <- Kall.list[ypos_bool]
  Kneg.list <- Kall.list[yneg_bool]
  pos_num <- length(Kpos.list)
  neg_num <- length(Kneg.list)
  for(i in 1:m){
    Kall.list[[i]] <- cbind(Kall.list[[i]], -Kall.list[[i]])
  }
  Kpos.list <- Kall.list[ypos_bool]
  Kneg.list <- Kall.list[yneg_bool]
  ### run LPBoost ###
  d <- rep(1/m, m)
  w <- rep(1, T) ## weight of weak classifiers
  alpha.list <- vector("list", T)
  dpos <- d[ypos_bool]
  dneg <- d[yneg_bool]
  gamma <- 0
  Kall_vals.mat <- matrix(NA, m, T)
  VAL_bef <- Inf  
  ell_ids <- rep(0, T)
  max_vals.list <- get_closed_vals(Kpos.list, Kneg.list, shape_num)
  neg_instance_nums <- instance_nums[yneg_bool]
  for(iter in 1:T){
    message("#Iteration =", iter)
    ### run DCA
    if(optimizer==1){
      sol <- MI_WeakLearn(neg_instance_nums, Kpos.list, Kneg.list, dpos, dneg, max_vals.list, shape_num, threshold, SHAPELET)
    }else{
      sol <- MI_WeakLearn_AdaGrad(neg_instance_nums, Kpos.list, Kneg.list, dpos, dneg, max_vals.list, shape_num, threshold, SHAPELET)      
    }
    alpha <- sol$alpha
    d_nonzero_id <- which(d!=0)
    Kall_vals <- rep(NA, m)
    alpha <- as.matrix(alpha, nrow=1)
    for(i in 1:m){
      Kall_vals[i] <- max(tcrossprod(alpha,Kall.list[[i]]))
    }
    alpha <- matrix(alpha[1:shape_num] - alpha[(shape_num+1):(shape_num*2)], nrow=1)
    alpha.list[[iter]] <- alpha
    Kall_vals.mat[, iter] <- Kall_vals
    VAL <- crossprod(t(y),(d*Kall_vals))
    if(abs(VAL) < gamma){
      message("LPBoost finished!")
      alpha.list <- alpha.list[1:(iter-1)]
      Kall_vals.mat <- Kall_vals.mat[ ,1:(iter-1)]
      break
    }
    if(VAL_bef == VAL){
      alpha.list <- alpha.list[1:(iter-1)]
      Kall_vals.mat <- Kall_vals.mat[ ,1:(iter-1)]
      break      
    }else{
      VAL_bef <- VAL
    }
    solution_LP <- run_LPBoost_dual(t(Kall_vals.mat[, 1:iter]), y, nu=nu*m)
    d <- solution_LP$d
    dpos <- d[ypos_bool]
    dneg <- d[yneg_bool]
    gamma <- solution_LP$gamma
  }
  solution_LP <- run_1normSVM(Kall_vals.mat, y, nu=nu*m)
  w <- solution_LP$w
  b <- solution_LP$b
  rho <- solution_LP$rho  
  train_time <- (proc.time()-time1)[3]
  message("rho=",rho)
  if(rho<0){
    restart_num <- restart_num - 1
    if(restart_num==0){
      message("The obtained solution is not meaningful.")
      return(NULL)
    }
    message("Restart with changing seed, maybe because the kmeans clustering failed.")
    return(MILIMS(original_bags.list, bag_labels, original_instance_nums, 
                  param.list, seed=seed*(restart_num+1), SHAPELET=SHAPELET, restart_num=restart_num))
  }
  return(list(alpha=alpha.list, Kx=shape.mat, w=w, b=b, rho=rho, ell=ell, kernel=kerneldot,train_time=train_time))
}


get_closed_vals <- function(Kpos.list, Kneg.list, shape_num){
  pos_num <- length(Kpos.list)
  neg_num <- length(Kneg.list)
  temp_Kpos_vals.mat <- matrix(NA, pos_num, shape_num*2)
  temp_Kneg_vals.mat <- matrix(NA, neg_num, shape_num*2)  
  temp_alpha_old <- rep(1, shape_num*2)
  temp_id <- 1
  for(i in 1:pos_num){
    temp_Kpos_vals.mat[i,] <- as.vector(rowMax(temp_alpha_old*t(Kpos.list[[i]])))
    temp_id <- temp_id + 1
  }  
  temp_id <- 1
  for(i in 1:neg_num){
    temp_Kneg_vals.mat[i,] <- as.vector(rowMax(temp_alpha_old*t(Kneg.list[[i]])))
    temp_id <- temp_id + 1
  }
  return(list(pos = temp_Kpos_vals.mat, neg = temp_Kneg_vals.mat))
}

get_closed_vals <- cmpfun(get_closed_vals)

MI_WeakLearn_MD <- function(neg_instance_nums.list, Kpos.list.list, Kneg.list.list, 
                            dpos, dneg, max_vals.list.list, ell_num, shape_nums, threshold, SHAPELET, optimizer){
  val <- Inf
  final_sol <- list()
  for(i in 1:ell_num){
    if(optimizer==1){
      sol <- MI_WeakLearn(neg_instance_nums.list[[i]],Kpos.list.list[[i]], Kneg.list.list[[i]], 
                          dpos, dneg, max_vals.list.list[[i]], shape_nums[i], threshold, SHAPELET)
    }else if(optimizer==2){
      sol <- MI_WeakLearn_AdaGrad(Kpos.list.list[[i]], Kneg.list.list[[i]], 
                                  dpos, dneg, max_vals.list.list[[i]], shape_nums[i], SHAPELET)
    }
    if(sol$optval < val){
      final_sol$alpha <- sol$alpha
      final_sol$ell_id <- i
      val <- sol$optval
    }
  }
  return(final_sol)
}

MI_WeakLearn <- function(neg_instance_nums,Kpos.list, Kneg.list, dpos, dneg, max_vals.list, shape_num, threshold, SHAPELET){
  pos_num <- length(dpos)
  neg_num <- length(dneg)
  ## Good Initialization
  dneg_nonzero_id <- which(dneg!=0)
  dpos_nonzero_id <- which(dpos!=0)  
  dneg_nonzero_logic <- (dneg!=0)
  dpos_nonzero_logic <- (dpos!=0)
  dpos_nonzero_num <- length(dpos_nonzero_id)
  dneg_nonzero_num <- length(dneg_nonzero_id)
  dpos_nonzero <- matrix(dpos[dpos_nonzero_id], nrow=1)
  dneg_nonzero <- matrix(dneg[dneg_nonzero_id], nrow=1)
  max_val <- -Inf
  vals <- ((crossprod(t(dpos_nonzero), max_vals.list$pos[dpos_nonzero_id,,drop=FALSE]))-(crossprod(t(dneg_nonzero),max_vals.list$neg[dneg_nonzero_id,,drop=FALSE])))
  max_id <- which.max(vals)
  max_val <- vals[max_id]
  alpha_old <- matrix(rep(0, shape_num*2), nrow=1)
  alpha_old[max_id] <- 1
  if(SHAPELET==TRUE){
    return(list(alpha=alpha_old, optval = -max_val))
  }
  Kpos_vals.mat <- matrix(NA, dpos_nonzero_num, shape_num*2)
  ## set up optimization problem
  mat2 <- c(rep(1, shape_num*2), rep(0, dneg_nonzero_num))
  if(ncol(dpos_nonzero)==0){
    return(list(alpha=alpha_old, optval = -max_val))
  }
  if(ncol(dneg_nonzero)==0){
    return(list(alpha=alpha_old, optval = -max_val))
  }
  mat1 <- make_sliding_matrix2(length(dneg_nonzero_id), neg_instance_nums[dneg_nonzero_id])
  Kneg.mat <- matlist2mat2(Kneg.list[dneg_nonzero_id], neg_instance_nums[dneg_nonzero_id])
  Amat <- (rbind(cbind(Kneg.mat, -mat1), mat2))
  dir <- c(rep("L", nrow(Amat)-1),"L")
  xopt_old <- 10000
  lprec <- make.lp(0, shape_num*2+dneg_nonzero_num)
  const_num <- nrow(Amat)
  for(i in 1:(const_num-1)){
    add.constraint(lprec, Amat[i,], "<=", 0)
  }  
  add.constraint(lprec, Amat[const_num,], "<=", 1)
  set.bounds(lprec, lower = c(rep(0, (shape_num*2)), rep(-Inf, length(dneg_nonzero_id))), 
             columns = c(1:(shape_num*2+dneg_nonzero_num)))
  for(k in 1:MAX_WLT){
    temp_id <- 1
    for(i in dpos_nonzero_id){
      Kpos_vals.mat[temp_id, ] <- Kpos.list[[i]][which.max((tcrossprod(alpha_old,Kpos.list[[i]]))),]
      temp_id <- temp_id + 1
    }
    set.objfn(lprec, c(-crossprod(t(dpos_nonzero),Kpos_vals.mat), dneg_nonzero))
    solve(lprec)
    solution <- get.variables(lprec) 
#    cvec <- c(-crossprod(t(dpos_nonzero),Kpos_vals.mat), dneg_nonzero)
#    solution <- Rcplex(cvec=cvec, Amat=Amat, bvec=bvec,
#                       lb=lb, control=list(trace=0, method=2))
    objval <- get.objective(lprec)
    if(xopt_old - objval <= threshold){
      alpha <- matrix(solution[1:(shape_num*2)], nrow=1)
      break
    }else{
      xopt_old <- objval
      alpha_old <- matrix(solution[1:(shape_num*2)], nrow=1)
    }
    if(k==MAX_WLT){
      alpha <- alpha_old
    }
  }
  return(list(alpha=alpha, optval = objval))
}

MI_WeakLearn <- cmpfun(MI_WeakLearn)

# projection to l1-ball provided by Boyd et al.
projsplx4AdaGrad = function(v, a, b=1){
  v[is.nan(v)] <- 0
  v[v<0] <- 0
  if(sum(abs(v))<1){
    return(v)
  }
  ord_id <- order(v/a, decreasing = TRUE)
  sv1 <- cumsum((a*v)[ord_id])
  sv2 <- (v/a)[ord_id]
  sv3 <- cumsum((a^2)[ord_id])
  rho <- which((sv1-(sv2*sv3)) < b)
  rho <- rho[length(rho)]
  theta <- (sv1[rho] - b) / sv3[rho] 
  w = v - theta*a
  w[w<0] <- 0
  ret <- w*a
  ret[is.nan(ret)] <- 0
  return(ret)
}

projsplx4AdaGrad <- cmpfun(projsplx4AdaGrad)


MI_WeakLearn_AdaGrad <- function(Kpos.list, Kneg.list, dpos, dneg, max_vals.list, shape_num, SHAPELET){
  eta <- ETA_adagrad
  pos_num <- length(dpos)
  neg_num <- length(dneg)
  dneg_nonzero_id <- which(dneg!=0)
  dpos_nonzero_id <- which(dpos!=0)  
  dneg_nonzero_logic <- (dneg!=0)
  dpos_nonzero_logic <- (dpos!=0)
  dpos_nonzero_num <- length(dpos_nonzero_id)
  dneg_nonzero_num <- length(dneg_nonzero_id)
  dpos_nonzero <- matrix(dpos[dpos_nonzero_id], nrow=1)
  dneg_nonzero <- matrix(dneg[dneg_nonzero_id], nrow=1)
  max_val <- -Inf
  vals <- ((crossprod(t(dpos_nonzero), max_vals.list$pos[dpos_nonzero_id,,drop=FALSE]))-(crossprod(t(dneg_nonzero),max_vals.list$neg[dneg_nonzero_id,,drop=FALSE])))
  max_id <- which.max(vals)
  max_val <- vals[max_id]
  alpha <- matrix(rep(0, shape_num*2), nrow=1)
  alpha[max_id] <- 1
  if(SHAPELET==TRUE){
    return(list(alpha=alpha, optval = -max_val))
  }
  fval_old <- Inf
  fvals <- rep(Inf, MAX_WLT+1)
  alpha.list <- vector("list", MAX_WLT+1)
  alpha.list[[1]] <- alpha
  s_seq <- rep(0, MAX_WLT)
  g <- g_k <- rep(0, shape_num*2)
  for(k in 1:MAX_GRAD_ITER){
    kpos_vec <- 0
    fval <- 0
    for(i in dpos_nonzero_id){
      vals <- (tcrossprod(alpha,Kpos.list[[i]]))
      max_id <- which.max(vals)
      fval <- fval + (vals[max_id] *dpos[i])
      kpos_vec <- kpos_vec + dpos[i] * Kpos.list[[i]][max_id,]
    }
    kneg_vec <- 0
    for(j in dneg_nonzero_id){
      vals <- (tcrossprod(alpha,Kneg.list[[j]]))
      max_id <- which.max(vals)
      fval <- fval - (vals[max_id] *dneg[j])
      kneg_vec <- kneg_vec - dneg[j] * Kneg.list[[j]][max_id,]
    }
    s_k <- -(kpos_vec + kneg_vec)
    g_k <- (s_k*s_k)
    g <- g + g_k
    H_k <- sqrt(g)
    v <- sqrt(H_k) - ((eta*s_k)/sqrt(H_k))
    fvals[k] <- -fval
    alpha.list[[k+1]] <- alpha <- projsplx4AdaGrad(v, 1/sqrt(H_k))
  }
  fval <- 0
  for(i in dpos_nonzero_id){
    vals <- (tcrossprod(alpha,Kpos.list[[i]]))
    max_id <- which.max(vals)
    fval <- fval + (vals[max_id] *dpos[i])
    kpos_vec <- kpos_vec + dpos[i] * Kpos.list[[i]][max_id,]
  }
  kneg_vec <- 0
  for(j in dneg_nonzero_id){
    vals <- (tcrossprod(alpha,Kneg.list[[j]]))
    max_id <- which.max(vals)
    fval <- fval - (vals[max_id] *dneg[j])
    kneg_vec <- kneg_vec - dneg[j] * Kneg.list[[j]][max_id,]
  }
  fvals[k+1] <- -fval
  fval_min_id <- which.min(fvals)
  alpha <- matrix(alpha.list[[fval_min_id]], nrow=1)
  return(list(alpha=alpha, optval=fvals[fval_min_id]))
}

MI_WeakLearn_AdaGrad <- cmpfun(MI_WeakLearn_AdaGrad)




### solution = (d1,....,dm), gamma
run_LPBoost_dual <- function(x, y, nu){
  x <- as.matrix(x)
  x <- rbind(x, -x)
  const_num <- nrow(x)
  m <- ncol(x)
  lprec <- make.lp(0, m+1)
  set.objfn(lprec, c(rep(0, m),  1))
  mat <- t(as.numeric(y)*t(x))
  for(i in 1:const_num){
    add.constraint(lprec, c(mat[i,],-1), "<=", 0)
  }
  add.constraint(lprec, c(rep(1,m),0), "=", 1)
  set.bounds(lprec, lower = c(rep(0, m), -Inf), columns = c(1:(m+1)))
  set.bounds(lprec, upper = c(rep((1/nu), m), Inf), columns = c(1:(m+1)))
  solve(lprec)
  solution <- get.variables(lprec)
  d <- solution[1:m]
  gamma <- solution[m+1]
  return(list(d=d, gamma=gamma))
}

run_LPBoost_dual <- cmpfun(run_LPBoost_dual)

run_1normSVM <- function(x, y, nu){
  x <- as.matrix(x)
  x <- cbind(x, -x)
  d <- ncol(x)
  m <- nrow(x)
  lprec <- make.lp(0, d+m+2)
  set.objfn(lprec, c(rep(0, d), 0, 1, rep(-1/nu, m)))
  lp.control(lprec,sense="max")  
  mat <- cbind(as.numeric(y)*x, matrix(y, m, 1), matrix(-1, m, 1), diag(1, m))
  for(i in 1:nrow(mat)){
    add.constraint(lprec, mat[i,], ">=", 0)    
  }
  add.constraint(lprec, c(rep(1,d), rep(0, m+2)), "=", 1)
  set.bounds(lprec, lower = c(rep(0, d), -Inf, -Inf, rep(0,m)), columns = c(1:(d+m+2)))
  set.bounds(lprec, upper = c(rep(1, d), Inf, Inf, rep(Inf,m)), columns = c(1:(d+m+2))) 
  solve(lprec)
  sol <- get.variables(lprec)
  w_plus <- sol[1:(d/2)]
  w_minus <- sol[((d/2)+1):d]
  w <- w_plus - w_minus
  b <- sol[d+1]
  rho <- sol[d+2]
  xi <- sol[(d+3):(d+2+m)]
  return(list(w=w, b=b, rho=rho))
}
