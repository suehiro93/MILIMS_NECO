library(zoo)
library(RColorBrewer)
library(gplots)

ALPHA_ADD <- 0.15

br250 <- gplots::colorpanel(n = 256, low = "skyblue", mid = "white", high = "red")
my_cm.colors <- function(num, alpha){
  adjustcolor(br250[1:128][1:num], alpha=alpha)
}

my_heat.colors <- function(num, alpha){
  adjustcolor(br250[256:129][1:num], alpha=alpha)
}

which.max.matrix <- function(mat){
  max_val <- -Inf
  argmax_row <- 0
  argmax_col <- 0  
  for(i in 1:ncol(mat)){
    for(j in 1:nrow(mat))
      if(mat[j,i] > max_val){
        max_val <- mat[j,i]
        argmax_col <- i
        argmax_row <- j
      }
  }
  return(list(row=argmax_row, col=argmax_col, val=max_val))
}

check_sparsity <- function(model){
  hypo_num <- length(which(model$w !=0))
  total_alpha_num <- hypo_num*ncol(model$alpha)
  nonzero_alpha_num <- length(which(model$alpha[which(model$w !=0),] !=0))
  message("hypo_num = ", hypo_num, " total_alpha_num = ", total_alpha_num, " nonzero_alpha_num = ", nonzero_alpha_num)
}

ret.topk.shape <- function(model.list, topk){
  Kx <- model.list$Kx
  alpha <- model.list$alpha
  w <- model.list$w
  mat <- w*alpha
  ord <- order(abs(mat),decreasing = TRUE)
  weight <- mat[ord[1:topk]]
  s <- vector("list", topk)
  for(i in 1:topk){
    maxrow <- ord[i]%%nrow(mat)
    maxcol <- ceiling(ord[i]/nrow(mat))
    message(mat[maxrow, maxcol])
    s[[i]] <- Kx[maxcol,]
  }
  return(list(w=weight, s=s))
}

ret.topk.shape2 <- function(model, topk){
  if(is.null(model$ell_ids)){
    model$ells <- model$ell#rep(model$ell, length(model$w))
    model$ell_ids <- rep(1, length(model$w))
    model$Kx.list <- list(model$Kx)
    model$sigma <- model$kernel@kpar$sigma
  }
  Kx <- model$Kx.list
  alpha <- matlist2mat(model$alpha)
  w <- model$w
  ell_ids <- model$ell_ids
  mat <- w*alpha
  ord <- order(abs(mat),decreasing = TRUE)
  weight <- mat[ord[1:topk]]
  s <- vector("list", topk)
  for(i in 1:topk){
    maxrow <- ceiling(ord[i]%%nrow(mat))
    if(maxrow==0){
      maxrow <- nrow(mat)
    }
    maxcol <- ceiling(ord[i]/nrow(mat))
    message(mat[maxrow, maxcol])
    s[[i]] <- Kx[[ell_ids[maxrow]]][maxcol,]
  }
  return(list(w=weight, s=s))
}

visualize_shapelets <- function(series, model.list, topk, filename=NULL, y_min=NULL, y_max=NULL, leg=TRUE){
  model <- ret.topk.shape2(model.list, topk)
  NAomit_id <- !is.na(series)
  series <- series[NAomit_id]
  s <- model$s
  w <- round(model$w, 3)
  ord <- order(model$w,decreasing = TRUE)  
  s <- s[ord]
  pos_num <- length(which(w>=0))
  neg_num <- length(w)-pos_num
  my_colour <- c(heat.colors(pos_num, alpha = 0.7), cm.colors(neg_num, alpha = 0.7))
  target_s <- s[1:topk]
  target_w <- (w[ord])[1:topk]
  if(is.null(y_min)){
    y_min <- min(c(unlist(target_s), series))
  }
  if(is.null(y_max)){
    y_max <- max(c(unlist(target_s), series))
  }
  x_min <- 1
  x_max <- length(series)
  Rs <- lengths(s)
  if(!is.null(filename)){
    if(leg==TRUE){
      pdf(filename,width = 12, height = 9, pointsize = 24)
    }else{
      pdf(filename,width = 10, height = 10, pointsize = 24)
    }
  }
  if(leg==TRUE){
    par(mar=c(3,3,3,3), ps=18)
  }else{
    par(mar=c(3,3,3,3), ps=18)
  }
  plot(series, type="l",lty=1, xlim=c(x_min, x_max), ylim=c(y_min, y_max), lwd=5, xlab="",ylab="")#xlab="time",ylab="value")
  par(new=TRUE)
  for(i in 1:length(target_s)){
    subsequences.mat <- rollapply(series, Rs[i], c)
    d <- sqrt(colSums((t(subsequences.mat) - target_s[[i]])^2))
    match_idx <- which.min(d)
    x <- match_idx:(match_idx+ length(target_s[[i]]) -1)
    plot(x, target_s[[i]],type="l", col=my_colour[i],xlim=c(x_min, x_max), ylim=c(y_min, y_max), lwd=7, xlab="",ylab="")
    par(new=TRUE, ps=18)
  }
  par(xpd=T)
  if(leg==TRUE){
    legend(par()$usr[1] + 0.2, par()$usr[4], legend=as.character(signif(target_w, digits = 3)),lty=1,lwd=7, col = my_colour, y.intersp=0.88)
  }
  if(!is.null(filename)){
    graphics.off()
  }
}


visualize_shapelets2 <- function(series, model, topk, filename=NULL, y_min=NULL, y_max=NULL, leg=TRUE){
  if(is.null(model$ell_ids)){
    model$ells <- model$ell#rep(model$ell, length(model$w))
    model$ell_ids <- rep(1, length(model$w))
    model$Kx.list <- list(model$Kx)
    model$sigma <- model$kernel@kpar$sigma
  }
  vals.list <- EvalFun4visualize(model, series, topk)
  ord_vals <- order(vals.list$vals, decreasing = TRUE)
  vals <- vals.list$vals[ord_vals]
  max_ids <- vals.list$max_ids[ord_vals]
  ells <- vals.list$ells[ord_vals]
  pos_num <- length(which(vals>=0))
  neg_num <- length(vals)-pos_num
  colours <- make_colours.mat(vals)
  x_min <- 1
  x_max <- length(series)
  if(!is.null(filename)){
    if(leg==TRUE){
      pdf(filename,width = 10, height = 7, pointsize = 24)
    }else{
      pdf(filename,width = 10, height = 10, pointsize = 24)
    }
  }
  if(leg==TRUE){
    par(mar=c(3,3,3,3), ps=18)#, bg="grey88")
  }else{
    par(mar=c(3,3,3,3), ps=18)#, bg="grey88")
  }
  plot(as.numeric(series), type="l",lty=1, xlim=c(x_min, x_max), ylim=c(y_min, y_max), lwd=5, xlab="time",ylab="value")
  par(new=TRUE)
  fin_id <- 1
  for(i in 1:length(vals)){
    if(is.na(vals[i])){
      fin_id <- i-1
      break
    }
    fin_id <- i
    x <- max_ids[i]:(max_ids[i]+ells[i]-1)
    y <- series[,max_ids[i]:(max_ids[i]+ells[i]-1)]
    plot(x, y,type="l", col=colours[i],xlim=c(x_min, x_max), ylim=c(y_min, y_max), lwd=10,xlab="time",ylab="value")
    par(new=TRUE, ps=18)
  }
  par(xpd=T)
  if(leg==TRUE){
    legend(par()$usr[1], par()$usr[4], legend=as.character(signif(vals[1:fin_id],digits=3)),lty=1,lwd=7, col = colours[1:fin_id])
  }
  if(!is.null(filename)){
    graphics.off()
  }
}

EvalFun4visualize <- function(hypos.list, X_test, topk){
  w <- hypos.list$w
  ord_w <- order(w, decreasing = TRUE)#[1:topk]
  w <- w#[ord_w]
  b <- hypos.list$b
  alpha.list <- hypos.list$alpha#[ord_w]
  Kx.list <- hypos.list$Kx.list#[ord_w]
  ell_ids <- hypos.list$ell_ids#[ord_w]
  ell_set <- hypos.list$ells
  hypo_num <- length(alpha.list)
  m <- nrow(X_test)
  noneed_id <- setdiff(1:(length(ell_set)), ell_ids)
  sigma <- hypos.list$sigma
  if(length(noneed_id)==0){
    noneed_id <- 0
  }
  ell_num <- length(ell_set)
  Kall.list.list <- vector("list", ell_num)
  for(r in 1:ell_num){
    kerneldot <- rbfdot(sigma/ell_set[r])
    if(is.element(r,noneed_id)){
      next
    }
    subseq_all.list <- gen_all_of_subsequences(X_test, ell_set[r])
    Kall.list.list[[r]] <- kernelBC(kerneldot, subseq_all.list, Kx.list[[r]])
  }
  Kall_vals <- rep(0, m)
  if(is.numeric(w)==FALSE){
    return(NULL)
  }
  #val.list < vector("list", length(w))
  vals <- rep(NA, hypo_num)
  max_ids <- rep(NA, hypo_num)
  ells <- rep(NA, hypo_num)
  nonzero_num <- 0
  for(j in 1:hypo_num){
    w_temp <- w[j]
    if(w_temp==0){
      next
    }
    alpha_temp <- alpha.list[[j]]
    Kall.list <- Kall.list.list[[ell_ids[j]]]
    #  browser()
    vals[j] <- w_temp*max(tcrossprod(alpha_temp,(Kall.list[[1]])))
    max_ids[j] <- which.max(tcrossprod(alpha_temp,(Kall.list[[1]])))
    ells[j] <- ell_set[ell_ids[j]]
    nonzero_num <- nonzero_num + 1
  }
  #    browser()
  vals <- vals + b/nonzero_num
  #vals <- vals
  ord_vals <- order(vals, decreasing = TRUE, na.last = TRUE)[1:topk]
  vals <- vals[ord_vals] #aaa#######################
  max_ids <- max_ids[ord_vals]
  ells <- ells[ord_vals]
  return(list(vals=vals[!is.na(vals)], max_ids=max_ids[!is.na(vals)], ells=ells[!is.na(vals)]))
}


EvalFun_TSMD2_4visualize <- function(hypos.list, X_test, topk){
  w <- hypos.list$w
  b <- hypos.list$b
  alpha.list <- hypos.list$alpha
  Kx.list <- hypos.list$Kx.list
  shape_nums <- hypos.list$shape_nums
  ells <- hypos.list$ells
  WS <- hypos.list$WS
  hypo_num <- length(alpha.list)
  m <- nrow(X_test)
  sigma <- hypos.list$sigma
  ell_num <- length(ells)
  bags.list.list <- vector("list", ell_num)
  instance_nums_kind <- rep(NA, ell_num)
  ell_all <- NULL
  for(r in 1:ell_num){
    mi.list <- ts_data2mi_data(X_test, NULL, ells[r])
    bags.list.list[[r]] <- mi.list$bags.list
    instance_nums_kind[r] <- mi.list$instance_nums[1]
    ell_all <- c(ell_all, rep(ells[r], nrow(Kx.list[[r]])))
  }
  Kall.list <- kernelBC2(sigma, bags.list.list, Kx.list, m, instance_nums_kind, shape_nums, WS)
  Kall_vals <- rep(0, m)
  if(is.numeric(w)==FALSE){
    return(NULL)
  }
  vals <- rep(NA, hypo_num)
  max_ids <- rep(NA, hypo_num)
  ells2 <- rep(NA, hypo_num)
  nonzero_num <- 0
  for(j in 1:hypo_num){
    w_temp <- w[j]
    if(w_temp==0){
      next
    }
    alpha_temp <- alpha.list[[j]]
    vals[j] <- w_temp*max(cppdot(alpha_temp,t(Kall.list[[1]])))
    max_ids[j] <- which.max(cppdot(alpha_temp,t(Kall.list[[1]])))
    ells2[j] <- ell_all[max_ids[j]]
    nonzero_num <- nonzero_num + 1
  }
  vals <- vals + b/nonzero_num
  ord_vals <- order(vals, decreasing = TRUE, na.last = TRUE)[1:topk]
  vals <- vals[ord_vals] 
  max_ids <- max_ids[ord_vals]
  ells2 <- ells2[ord_vals]
  return(list(vals=vals[!is.na(vals)], max_ids=max_ids[!is.na(vals)], ells=ells2[!is.na(vals)]))
}


make_colours.mat <- function(vals, alpha_add=ALPHA_ADD){
  pos_num <- length(which(vals>=0))
  neg_num <- length(vals)-pos_num
  absvals <- abs(vals[!is.na(vals)])
  norm_vals <- vals/(sum(absvals))
  colours.mat <- matrix(NA, length(vals), length(vals))
  for(i in 1:nrow(colours.mat)){
    if(is.na(norm_vals[i])){
      break
    }
    if(pos_num == 0){
      colours.mat[i,] <- c(my_cm.colors(neg_num, alpha = min(c(abs(norm_vals[i])+alpha_add, 1))))
    }else if(neg_num == 0){
      colours.mat[i,] <- c(my_heat.colors(pos_num, alpha = (min(c(abs(norm_vals[i])+alpha_add, 1)))))
    }else{
      colours.mat[i,] <- c(my_heat.colors(pos_num, alpha = (min(c(abs(norm_vals[i])+alpha_add, 1)))), my_cm.colors(neg_num, alpha = min(c(abs(norm_vals[i])+alpha_add, 1))))
    }
  }
  colours <- diag(colours.mat)
  colours <- colours[!is.na(colours)]
  return(colours)
}