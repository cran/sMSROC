auc_ci_boot <- function(marker,
                        outcome,
                        status,
                        observed.time,
                        left,
                        right,
                        time,
                        data_type,
                        meth,
                        grid,
                        probs,
                        ci.cl,
                        ci.nboots,
                        parallel,
                        ncpus,
                        all){
  bootstrap_auc <- function (i,
                             marker,
                             outcome,
                             status,
                             observed.time,
                             left,
                             right,
                             time,
                             data_type,
                             meth,
                             grid,
                             probs,
                             all){
  boots.idx <- sample(1:length(marker), length(marker), replace = T)
    if (data_type == "binout"){
      boots.auc <- sMS_binout(marker[boots.idx],
                              status[boots.idx],
                              meth,
                              grid,
                              probs,
                              all)
    } else if (data_type == "timerc"){
      boots.auc <- sMS_timerc(marker[boots.idx],
                              status[boots.idx],
                              observed.time[boots.idx],
                              outcome[boots.idx],
                              time,
                              meth,
                              grid,
                              probs,
                              all)
    } else {
      boots.auc <- sMS_timeic(marker[boots.idx],
                              left[boots.idx],
                              right[boots.idx],
                              outcome[boots.idx],
                              time,
                              meth,
                              grid,
                              probs,
                              all)
    }
    ret <- boots.auc$auc
  }
  auc.all <-NULL
  i <- NULL
  if (as.logical(parallel)){
    clust <- makeCluster(ncpus)
    auc.all <- foreach(i = 1:ci.nboots,
                       .combine=rbind,
                       .multicombine = TRUE,
                       .packages = "parallel") %dopar%{
    auc.all[i] <- bootstrap_auc(i,
                                marker = marker,
                                outcome = outcome,
                                status = status,
                                observed.time = observed.time,
                                left = left,
                                right = right,
                                time = time,
                                data_type = data_type,
                                meth = meth,
                                grid = grid,
                                probs = probs,
                                all = all)
    }
    stopCluster(clust)
    } else {
      auc.all <- sapply(1:ci.nboots,
                        bootstrap_auc,
                        marker = marker,
                        outcome = outcome,
                        status = status,
                        observed.time = observed.time,
                        left = left,
                        right = right,
                        time = time,
                        data_type = data_type,
                        meth = meth,
                        grid = grid,
                        probs = probs,
                        all = all)
      }
      ic.l     <- round(quantile(auc.all, ci.cl),3)
      ic.u     <- round(quantile(auc.all, 1 - ci.cl),3)
      ret      <- list()
      ret$ic.l <- ic.l
      ret$ic.u <- ic.u
      ret
}
auc_ci_empr <- function(SE, SP, auc, probs, controls, cases, ci.cl){
  fuNint <- function(a, b){
    nint <- sapply(2:length(a), function(i){(a[i] - a[i-1]) * b[i-1]})
    nint.all <- sum(nint)
  }
  mean.probs <- mean(probs)
  var.empr   <- (fuNint(SP, SE^2) - (fuNint(SP, SE))^2 +
                ((1 - mean.probs)/mean.probs) * (fuNint(1-SE, SP^2) - fuNint(1-SE, SP)^2))^ 0.5 /
                ((1 - mean.probs) * (cases + controls)) ^ 0.5
  ic.l <- round((auc + var.empr * qnorm(ci.cl)), 3)
  ic.u <- round((auc - min((var.empr * qnorm(ci.cl)),1)), 3)
  ret      <- list()
  ret$ic.l <- ic.l
  ret$ic.u <- ic.u
  ret
}
auc_ci_nvar <- function(marker,
                        outcome,
                        status,
                        observed.time,
                        left,
                        right,
                        time,
                        meth,
                        data_type,
                        grid,
                        probs,
                        sd.probs,
                        ci.cl,
                        nboots,
                        SE,
                        SP,
                        auc,
                        parallel,
                        ncpus,
                        all){
  if (missing(sd.probs)){
    sd.probs   <- variance_probs(marker,
                                 outcome,
                                 status,
                                 observed.time,
                                 left,
                                 right,
                                 time,
                                 meth,
                                 data_type,
                                 grid,
                                 probs,
                                 nboots,
                                 parallel,
                                 ncpus,
                                 all)
    sd.matr    <- splinefun(sort(marker), sd.probs$sd.probs)(sort(marker))
  } else {
    sd.matr    <- sd.probs
  }
  mean.probs <- mean(probs)
  E1 <- (1- mean.probs)^(-1) * (SE - auc) * (1 - probs)
  E2 <-  mean.probs^(-1) * (SP - auc) * probs
  E3 <- (mean.probs^(-1) * (SP - auc) - (1- mean.probs)^(-1) * (SE - auc)) * sd.matr
  var.nvar <- (mean((E1-E2) ^ 2 ) + mean(E3)^2 )^0.5 / sqrt(length(marker))
  ic.l <- round(auc + qnorm(ci.cl) * var.nvar, 3)
  ic.u <- round(auc - min(qnorm(ci.cl) * var.nvar, 1), 3)
  ret      <- list()
  ret$ic.l <- ic.l
  ret$ic.u <- ic.u
  ret
}
variance_probs <- function(marker,
                           outcome,
                           status,
                           observed.time,
                           left,
                           right,
                           time,
                           meth,
                           data_type,
                           grid,
                           probs,
                           ci.nboots,
                           parallel,
                           ncpus,
                           all){
  if (missing(all)){
    all <- "T"
  }
  sd.all <- NULL
  i <- NULL
  bootstrap_var <- function (i,
                             marker,
                             outcome,
                             status,
                             observed.time,
                             left,
                             right,
                             time,
                             meth,
                             data_type,
                             grid,
                             probs,
                             all){
  boots.var <- NULL
  boots.idx <- sample(1:length(marker), length(marker), replace = T)
    if (data_type == "binout"){
      boots.var <- sMS_binout(marker[boots.idx],
                              outcome[boots.idx],
                              meth,
                              grid,
                              probs,
                              all)
    } else if (data_type == "timerc"){
      boots.var <- sMS_timerc(marker[boots.idx],
                              status[boots.idx],
                              observed.time[boots.idx],
                              outcome[boots.idx],
                              time,
                              meth,
                              grid,
                              probs,
                              all)
    } else {
      boots.var <- sMS_timeic(marker[boots.idx],
                              left[boots.idx],
                              right[boots.idx],
                              outcome[boots.idx],
                              time,
                              meth,
                              grid,
                              probs,
                              all)
    }
    Iv   <- which(!is.na(boots.var$probs))
    if (length(Iv) > 0) {
      fu_P <- approxfun(boots.var$marker[Iv], boots.var$probs[Iv])(sort(marker))
    } else {
      fu_P <- rep(NA, length(marker))
    }
  }
  if (parallel == "T"){
    clust  <- makeCluster(ncpus)
    sd.all <- foreach(i = 1:ci.nboots,
                      .combine=rbind,
                      .multicombine = TRUE,
                      .packages = "parallel") %dopar%{
      sd.all[i] <- bootstrap_var(i,
                                 marker = marker,
                                 outcome = outcome,
                                 status = status,
                                 observed.time = observed.time,
                                 left = left,
                                 right = right,
                                 time =time,
                                 data_type = data_type,
                                 meth = meth,
                                 grid = grid,
                                 probs = probs,
                                 all = all)
      }
      stopCluster(clust)
      } else {
        sd.all  <- sapply(1:ci.nboots,
                          bootstrap_var,
                          marker = marker,
                          outcome = outcome,
                          status = status,
                          observed.time = observed.time,
                          left = left,
                          right = right,
                          time = time,
                          data_type = data_type,
                          meth = meth,
                          grid = grid,
                          probs = probs,
                          all = all)
        sd.all = t(sd.all)
      }
      if (!is.null(sd.all)){
        ret <- list()
        ret$sd.probs <- apply(sd.all, 2, sd, na.rm = TRUE)
        ret
      }
}
