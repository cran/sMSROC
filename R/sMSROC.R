sMSROC <- function (marker = NULL,
                    status = NULL,
                    observed.time = NULL,
                    left = NULL,
                    right = NULL,
                    time = 1,
                    meth = c("L", "S", "E"),
                    grid,
                    probs,
                    sd.probs,
                    conf.int = c("F", "T"),
                    ci.cl,
                    ci.meth = c("E", "V", "B"),
                    ci.nboots = 500,
                    parallel = c("F", "T"),
                    ncpus = 1,
                    all = c("T", "F") ){
  data_meth     <- match.arg(meth)
  data_meth     <- check_meth(data_meth, probs)$meth
  data_conf.int <- match.arg(conf.int)
  data_ci.meth  <- match.arg(ci.meth)
  data_parallel <- match.arg(parallel)
  data_type     <- check_type_outcome(status, observed.time, left, right)
  data_grid     <- check_grid(grid)
  data_all      <- match.arg(all)
  data_ci       <- check_conf_int(data_conf.int,
                                  ci.cl,
                                  data_ci.meth,
                                  ci.nboots,
                                  data_parallel,
                                  ncpus)
  if (data_type$type.outcome == "binout"){
    data_time <- NULL
    data_sample <- check_marker_binout(marker, status, probs, sd.probs)
    sMSROC <- sMS_binout(data_sample$marker,
                         data_sample$outcome,
                         data_meth,
                         data_grid$grid,
                         data_sample$probs,
                         data_all)
    } else{
      data_time <- check_tim(time)$time
        if (data_type$type.outcome == "timerc"){
          data_sample <- check_marker_timerc(marker,
                                             status,
                                             observed.time,
                                             data_time,
                                             probs,
                                             sd.probs)
          sMSROC <- sMS_timerc(data_sample$marker,
                               data_sample$status,
                               data_sample$observed.time,
                               data_sample$outcome,
                               data_time,
                               data_meth,
                               data_grid$grid,
                               data_sample$probs,
                               data_all)
          } else if (data_type$type.outcome == "timeic"){
            data_sample <- check_marker_timeic(marker,
                                               left,
                                               right,
                                               data_time,
                                               probs,
                                               sd.probs)
            sMSROC <- sMS_timeic(data_sample$marker,
                                 data_sample$left,
                                 data_sample$right,
                                 data_sample$outcome,
                                 data_time,
                                 data_meth,
                                 data_grid$grid,
                                 data_sample$probs,
                                 data_all)
            }
      }
      if(!is.null(data_grid$message)){
        data_sample$message <- rbind(data_sample$message, data_grid$message)
      }
      if(!is.null(data_ci$message)){
        data_sample$message <- rbind(data_sample$message, data_ci$message)
      }
      if (data_conf.int == "T"){
        ci.cl <- (1 - data_ci$ci.cl) / 2
        if (data_ci$ci.meth == "E"){
          auc.ci <- auc_ci_empr(sMSROC$SE,
                                sMSROC$SP,
                                sMSROC$auc,
                                sMSROC$probs,
                                data_sample$controls,
                                data_sample$cases,
                                ci.cl)
        } else if (data_ci$ci.meth == "B"){
          auc.ci <- auc_ci_boot(marker = data_sample$marker,
                                outcome = data_sample$outcome,
                                status = data_sample$status,
                                observed.time = data_sample$observed.time,
                                left = data_sample$left,
                                right = data_sample$right,
                                time = data_time,
                                data_type = data_type$type.outcome,
                                meth = data_meth,
                                grid = data_grid$grid,
                                probs = data_sample$probs,
                                ci.cl = ci.cl,
                                ci.nboots = data_ci$nboots,
                                parallel = data_parallel,
                                ncpus = data_ci$ncpus,
                                all = data_all)
        } else {
          auc.ci <- auc_ci_nvar(marker = data_sample$marker,
                                outcome = data_sample$outcome,
                                status = data_sample$status,
                                observed.time = data_sample$observed.time,
                                left = data_sample$left,
                                right = data_sample$right,
                                time = data_time,
                                meth = data_meth,
                                data_type = data_type$type.outcome,
                                grid = data_grid$grid,
                                ci.cl = ci.cl,
                                nboots = data_ci$nboots,
                                SE = sMSROC$SE,
                                SP = sMSROC$SP,
                                auc = sMSROC$auc,
                                probs = sMSROC$probs,
                                parallel = data_parallel,
                                ncpus = data_ci$ncpus,
                                all = data_all)
        }
      } else {
        auc.ci <- NULL
      }
      if (!is.null(sMSROC$SE)){
        data <- list()
        data$type <- data_type
        data$grid <- data_grid$grid
        data$marker    <- data_sample$marker
        data$outcome   <- data_sample$outcome
        data$ncpus     <- data_ci$ncpus
        data$ci.nboots <- data_ci$ci.nboots
        data$parallel  <- data_parallel
        data$meth      <- data_meth
        if (data_type == "timerc"){
          data$status <- data_sample$status
          data$observed.time <- data_sample$observed.time
        } else {
          data$status <- NULL
          data$observed.time <- NULL
        }
        if (data_type == "timeic"){
          data$left  <- data_sample$left
          data$right <- data_sample$right
        } else {
          data$left  <- NULL
          data$right <- NULL
        }
        ret <- list()
        ret$thres    <- sMSROC$marker
        ret$SE       <- sMSROC$SE
        ret$SP       <- sMSROC$SP
        ret$probs    <- sMSROC$probs
        ret$u        <- sMSROC$u
        ret$ROC      <- sMSROC$ROC
        ret$auc      <- sMSROC$auc
        ret$auc.ci.l <- auc.ci$ic.l
        ret$auc.ci.u <- auc.ci$ic.u
        ret$ci.cl    <- data_ci$ci.cl
        ret$ci.meth  <- ci.meth
        ret$time     <- data_time
        ret$data     <- data
        ret$message  <- data_sample$message
        class(ret) <- "sMSROC"
        ret
      } else {
        stop(message("Non results to be shown"))
     }
}


