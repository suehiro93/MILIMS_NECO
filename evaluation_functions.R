evalFun <- function(hypos.list, test_bags.list){
  w <- hypos.list$w
  b <- hypos.list$b
  alpha.list <- hypos.list$alpha
  Kx <- hypos.list$Kx
  kerneldot <- hypos.list$kernel
  hypo_num <- length(alpha.list)
  m <- length(test_bags.list)
  Kall_vals <- rep(0, m)
  if(is.numeric(w)==FALSE){
    return(NULL)
  }
  Kall.list <- kernelBC(kerneldot, test_bags.list, Kx)
  for(j in 1:hypo_num){
    w_temp <- w[j]
    if(w_temp==0){
      next
    }
    alpha_temp <- alpha.list[[j]]
    for(i in 1:m){
      Kall_vals[i] <- Kall_vals[i] + w_temp*max(tcrossprod(alpha_temp, (Kall.list[[i]])))
    }
  }
  Kall_vals <- Kall_vals+b
  return(Kall_vals)
}

evalFun <- cmpfun(evalFun)

evalFun_TS <- function(hypos.list, X_test){
  w <- hypos.list$w
  b <- hypos.list$b
  alpha.list <- hypos.list$alpha
  Kx.list <- hypos.list$Kx.list
  kerneldot_original <- hypos.list$kerneldot
  ell_ids <- hypos.list$ell_ids
  ells <- hypos.list$ells
  hypo_num <- length(alpha.list)
  m <- nrow(X_test)
  noneed_id <- setdiff(1:(length(ells)), ell_ids)
  sigma <- hypos.list$sigma
  if(length(noneed_id)==0){
    noneed_id <- 0
  }
  ell_num <- length(ells)
  Kall.list.list <- vector("list", ell_num)
  for(r in 1:ell_num){
    kerneldot <- kerneldot_original(sigma/ells[r])
    if(is.element(r,noneed_id)){
      next
    }
    bags.list <- ts_data2mi_data(X_test, NULL, ells[r])$bags.list
    Kall.list.list[[r]] <- kernelBC(kerneldot, bags.list, Kx.list[[r]])
  }
  Kall_vals <- rep(0, m)
  if(is.numeric(w)==FALSE){
    return(NULL)
  }
  for(j in 1:hypo_num){
    w_temp <- w[j]
    if(w_temp==0){
      next
    }
    alpha_temp <- alpha.list[[j]]
    Kall.list <- Kall.list.list[[ell_ids[j]]]
    for(i in 1:m){
      Kall_vals[i] <- Kall_vals[i] + w_temp*max(tcrossprod(alpha_temp,Kall.list[[i]]))
    }
  }
  Kall_vals <- Kall_vals+b
  return(Kall_vals)
}

evalFun_TS <- cmpfun(evalFun_TS)

calc_accuracy4TS <- function(hypos.list, x_test, y_test){
  m <- length(y_test)
  val <- evalFun_TS(hypos.list, x_test)
  if(is.null(val)){
    return(0)
  }
  acc <- length(which((val)*(y_test)>0)==TRUE)/length(y_test)
  return(acc)
}



calc_accuracy <- function(hypos.list, test_bags.list, test_bag_labels){
  m <- length(test_bag_labels)
  val <- evalFun(hypos.list, test_bags.list)
  if(is.null(val)){
    return(0)
  }
  acc <- length(which((val)*(test_bag_labels)>0)==TRUE)/m
  return(accuracy=acc)
}

