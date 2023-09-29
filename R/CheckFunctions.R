# Type of outcome
check_type_outcome <- function(status = NULL,
                               observed.time = NULL,
                               left = NULL,
                               right = NULL){
  type.outcome <- NULL
  if (missing(status) | is.null(status)){
    if((missing(left)) | (missing(right))){
      type.outcome <- "unknow"
    } else{
      type.outcome <- "timeic"
    }
  } else if (missing(observed.time) | is.null(observed.time)){
      type.outcome <- "binout"
  } else {
      type.outcome <- "timerc"
  }
  ret <- list()
  ret$type.outcome <- type.outcome
  ret
}
check_marker_binout <- function (marker, status, probs, sd.probs){
  mt <- NULL
  if (missing(status) | is.null(status) | sum(!is.na(status))== 0){
    m <- "Response data should be indicated."
    stop(m)
  } else {
    level <- (levels(as.factor(status)))
    l <- length(level)
    if (l < 2){
      m <- "Response data must contain at least two different values."
      stop(m)
    } else if (l > 2){
      m <- "There are more than two different values in response data.
            Only the two lowest will be considered."
      mt <- rbind(mt, c("outcome", m))
    }
  }
  if (missing(marker) | is.null(marker) | sum(!is.na(marker))== 0){
    m <- "Marker data should be indicated."
    stop(m)
  } else if (!is.numeric(marker)){
    m <- "Marker data are not numeric."
    stop(m)
  }
  if(length(status) != length(marker)){
    m <- "Response and marker vectors should have the same length."
    stop(m)
  }
  if (!missing(probs)){
    if (!(is.numeric(probs))){
      m <- "Probabilities should be a numeric vector."
      stop(m)
    }
    if (length(probs) != length(marker)){
      m <- "The predictive model and marker vectors should have the same length."
      stop(m)
    }
    if ((sum(probs > 1) + sum (probs < 0)) > 0){
      m <- "The predictive model vector contains non valid values."
    }
    if (!missing(sd.probs)){
      if (length(probs) != length(sd.probs)){
        m <- "The predictive model and the deviation vectors should have the same length."
        stop(m)
      }
    }
  }
  Im <- which(!is.na(marker))
  marker.nm  <- marker[Im]
  status.nm  <- status[Im]
  if (length(marker) > length(marker.nm)){
    m  <- "Observations with missing Marker values have been removed."
    mt <- rbind (mt, c("marker",m))
  }
  neg <- split(marker.nm, f = status.nm)[[level[1]]]
  pos <- split(marker.nm, f = status.nm)[[level[l]]]
  Ims <- which(is.na(status.nm))
  n0  <- length(neg)
  n1  <- length(pos)
  nm  <- length(Ims)
  l0  <- ifelse (level[1] == 0, level[1],0)
  l1  <- ifelse (level[l] == 1, level[l],1)
  marker.bm  <- c(neg, pos)
  status.bm <- c(rep(l0, n0), rep(l1, n1))
  if(!missing(probs)){
    probs.nm <- probs[Im]
    negp  <- split(probs.nm, f = status.nm)[[level[1]]]
    posp  <- split(probs.nm, f = status.nm)[[level[l]]]
    probs <- c(negp, posp, probs.nm[Ims])
    probs <- as.numeric(probs)
      if (!missing(sd.probs)){
        sd.probs.nm <- sd.probs[Im]
        sd.negp  <- split(sd.probs.nm, f = status.nm)[[level[1]]]
        sd.posp  <- split(sd.probs.nm, f = status.nm)[[level[l]]]
        sd.probs <- c(sd.negp, sd.posp, sd.probs.nm[Ims])
        sd.probs <- as.numeric(sd.probs)
      }
    } else {
      probs <- NULL
      sd.probs <- NULL
    }
    marker <- c(marker.bm, marker.nm[Ims])
    status <- c(status.bm, rep(-1, length = nm))
    if (missing(probs)){
      probs <- NULL
    }
    if (missing(sd.probs)){
      sd.probs <- NULL
    }
  ret <- list()
  ret$marker   <- as.numeric(marker)
  ret$outcome  <- as.numeric(status)
  ret$probs    <- probs
  ret$sd.probs <- sd.probs
  ret$controls <- n0
  ret$cases    <- n1
  ret$misout   <- nm
  ret$message  <- mt
  ret
}
check_marker_timerc <- function (marker, status, observed.time, time, probs, sd.probs){
  mt <- NULL
  if (missing(status) | is.null(status) | sum(!is.na(status))== 0){
    m <- "Status data should be indicated."
    stop(m)
  } else{
    level <- (levels(as.factor(status)))
    l <- length(level)
    if (l < 2){
      m <- "Status data must contain at least two different values."
      stop(m)
    } else if (l > 2){
      m <- "There are more than two different values in response data.
            Only the two lowest will be considered."
      mt <- rbind(mt, c("Status", m))
    }
  }
  if (missing(marker) | is.null(marker) | sum(!is.na(marker))== 0){
    m <- "Marker data should be indicated."
    stop(m)
  } else if(!is.numeric(marker)){
    m <- "Marker data are not numeric."
    stop(m)
  }
  if (length(status) != length(marker)){
    m <- "Status and marker vectors should have the same length."
    stop(m)
  }
  if (length(marker) != length(observed.time)){
    m <- "Observed time vector length should be the same as marker and status vectors."
    stop(m)
  }
  Im <- which(!is.na(marker))
  marker.bm <- marker[Im]
  status.bm <- status[Im]
  observed.time.bm <- observed.time[Im]
  if (length(observed.time.bm) < length(observed.time)){
    m <- "Observations with missing Marker values have been removed."
    mt <- rbind (mt, c("Marker",m))
  }
  if (missing(time) | is.null(time) | sum(!is.na(time)) == 0){
    m <- "The point of time should be indicated."
    stop(m)
  } else if (!is.numeric(time)){
    m <- "Marker data are not numeric."
    stop(m)
  }
  if (!missing(probs)){
    probs.bm <- probs[Im]
    if (!missing(sd.probs)){
      sd.probs.bm <- sd.probs[Im]
    }
  }
  # Missing status values
  Ims       <- which(!is.na(status.bm))
  status.bm <- status.bm[Ims]
  marker.bm <- marker.bm[Ims]
  observed.time.bm <- observed.time.bm[Ims]
  if(!missing(probs)){
    probs.bm <- probs.bm[Ims]
    if (!missing(sd.probs)){
      sd.probs.bm <- sd.probs.bm[Ims]
    }
  }
  # Cases, controls and censored observations
  neg   <-  split(marker.bm, f = status.bm)[[level[1]]]
  pos   <-  split(marker.bm, f = status.bm)[[level[l]]]
  negt  <-  split(observed.time.bm, status.bm)[[level[1]]]
  post  <-  split(observed.time.bm, status.bm)[[level[l]]]
  l0    <-  ifelse (level[1] == 0, level[1],0)
  l1    <-  ifelse (level[l] == 1, level[l],1)
  marker.bm <- c(neg, pos)
  status.bm <- c(rep(l0, length(neg)), rep(l1, length(pos)))
  observed.time.bm <- c(negt, post)
  if(!missing(probs)){
    negp <- split(probs.bm, f = status.bm)[[level[1]]]
    posp <- split(probs.bm, f = status.bm)[[level[l]]]
    probs.bm <- c(negp, posp)
    M    <- data.frame(marker.bm, status.bm, observed.time.bm, probs.bm)
    if (!missing(sd.probs)){
      negsdp <- split(sd.probs.bm, f = status.bm)[[level[1]]]
      possdp <- split(sd.probs.bm, f = status.bm)[[level[l]]]
      sd.probs.bm <- c(negsdp, possdp)
      M    <- data.frame(marker.bm, status.bm, observed.time.bm, probs.bm, sd.probs.bm)
    }
  } else {
    M <- data.frame(marker.bm, status.bm, observed.time.bm)
  }
  Ineg <- which(M[,3] > time)
  Ipos <- which(M[,3] <= time & M[,2] == l1)
  Icen <- which(M[,3] <= time & M[,2] == l0)
  marker <- as.numeric(c(M[Ineg,1], M[Ipos,1], M[Icen, 1]))
  status <- as.numeric(c(M[Ineg,2], M[Ipos,2], M[Icen, 2]))
  observed.time <- as.numeric(c(M[Ineg, 3], M[Ipos,3], M[Icen, 3]))
  if(!missing(probs)){
    probs <- as.numeric(c(M[Ineg, 4], M[Ipos,4], M[Icen, 4]))
      if (!missing(sd.probs)){
        sd.probs <- as.numeric(c(M[Ineg, 5], M[Ipos,5], M[Icen, 5]))
      }
    }
  outcm  <- as.numeric(c(rep(0, length(Ineg)), rep(1, length(Ipos)), rep(-1, length(Icen))))
  if (missing(probs)){
    probs <- NULL
  }
  if (missing(sd.probs)){
    sd.probs <- NULL
  }
  ret <- list()
  ret$marker <- marker
  ret$status <- status
  ret$observed.time <- observed.time
  ret$probs    <- probs
  ret$sd.probs <- sd.probs
  ret$outcome  <- outcm
  ret$controls <- length(Ineg)
  ret$cases    <- length(Ipos)
  ret$misout   <- length(Icen)
  ret$message  <- mt
  ret
}
check_marker_timeic <- function (marker, left, right, time, probs, sd.probs){
  mt <- NULL
  if (missing(left) | is.null(left) | sum(!is.na(left)) == 0){
    m <- "Left/Right interval edges should be indicated."
    stop(m)
  } else {
    if (missing(right) | is.null(right)){
      m <- "Right/Left interval edges should be indicated."
      stop(m)
    } else {
      if (sum(left > (ifelse(is.na(right), (max(left) + 1 ), right))) > 0){
        m <- "Non-valid data. Interval left side greater than right side"
        stop(m)
      } else {
        if(length(left) != length(right)){
          m <- "Left and right edges vectors should have the same length."
          stop(m)
        }
      }
    }
  }
  if (missing(marker) | is.null(marker) | sum(!is.na(marker))== 0){
    m <- "Marker data should be indicated."
    stop(m)
  } else if(!is.numeric(marker)){
    m <- "Marker data are not numeric."
    stop(m)
  }
  if(length(left) != length(marker)){
    m <- "Intervals and marker vectors should have the same length."
    stop(m)
  }
  Im     <- which(!is.na(marker))
  marker <- marker[Im]
  left   <- left[Im]
  right  <- right[Im]
  if (!missing(probs)){
    probs <- probs[Im]
    M     <- cbind(left, right, marker, probs)
    if (!missing(sd.probs)){
      sd.probs<- sd.probs[Im]
      M     <- cbind(M, sd.probs)
    }
  } else {
    M <- cbind(left, right, marker)
  }
  if (!(length(Im) == length (marker))){
    m    <- "Observations with missing Marker values have been removed."
    mt   <- rbind (mt, c("Marker",m))
  }
  Ineg   <- which(M[,1]  >= time)
  Ipos   <- which(M[,2]  <= time)
  Iund   <- which((M[,1] < time) & ( time < M[,2] | is.na(M[,2])))
  marker <- as.numeric(c(M[Ineg,3], M[Ipos,3], M[Iund,3]))
  left   <- as.numeric(c(M[Ineg,1], M[Ipos,1], M[Iund,1]))
  right  <- as.numeric(c(M[Ineg,2], M[Ipos,2], M[Iund,2]))
  outcm  <- c(rep(0, length(Ineg)), rep(1, length(Ipos)), rep(-1, length(Iund)))
  if (!missing(sd.probs)){
    probs <- as.numeric(c(M[Ineg,4], M[Ipos,4], M[Iund,4]))
    if (!missing(sd.probs)){
      probs <- as.numeric(c(M[Ineg,5], M[Ipos,5], M[Iund,5]))
    }
  } else {
    probs <- NULL
    sd.probs <- NULL
  }
  ret <- list()
  ret$marker   <- marker
  ret$left     <- left
  ret$right    <- right
  ret$probs    <- probs
  ret$sd.probs <- sd.probs
  ret$outcome  <- outcm
  ret$controls <- length(Ineg)
  ret$cases    <- length(Ipos)
  ret$misout   <- length(Iund)
  ret$message  <- mt
  ret
}
check_grid <- function(grid){
  mt <- NULL
  if (missing(grid)){
    grid <- 1000
    m <- "No grid indicated. 1000 grid is assumed."
    mt <- rbind(mt, c("grid", m))
  } else if (!((is.numeric(grid)) | grid > 0)){
    m <- paste(grid, "-", "Invalid grid selection. ")
    stop(m)
  }
  ret <- list()
  ret$grid    <- as.numeric(grid)
  ret$message <- mt
  ret
}
check_conf_int <- function(conf.int, ci.cl, ci.meth, ci.nboots, parallel, ncpus){
  mt <- NULL
  if (conf.int == "T"){
    data_ci.cl <- check_ci_cl(ci.cl)
    ci.cl <- data_ci.cl$ci.cl
    mt    <- rbind(mt, data_ci.cl$mt)
    data_ci.nboots <- check_nboots(ci.nboots)
    ci.nboots  <- data_ci.nboots$nboots
    mt <- rbind(mt, data_ci.nboots$mt)
    if (parallel == "T"){
      check_ncpus(ncpus)
    }
  } else {
    m <- "No Confidence Intervals will be computed."
    mt <- rbind(mt, c("conf.int",m))
    ci.cl     <- NULL
    ci.meth   <- NULL
    ci.nboots <- NULL
    ncpus     <- NULL
  }
  ret <- list()
  ret$ci.cl   <- ci.cl
  ret$ci.meth <- ci.meth
  ret$nboots  <- ci.nboots
  ret$ncpus   <- ncpus
  ret$message <- mt
  ret
}
check_tim <- function (time){
  if (missing(time) | is.null(time) | sum(!is.na(time))== 0){
    m <- "Time data should be indicated."
    stop(m)
  } else if(!is.numeric(time)){
    m <- "Time data are not numeric."
    stop(m)
  } else if(time < 0){
    m <- "Invalid time."
    stop(m)
  }
  ret <- list()
  ret$time <- time
  ret
}
check_meth <- function (meth, probs){
  if (missing(probs)){
    meth <- meth
  } else {
    meth <- "M"
  }
  ret <- list()
  ret$meth <- meth
  ret
}
check_ncpus <- function(ncpus){
  if(!(is.integer(ncpus) | ncpus > 0)){
    m <- "Invalid ncpus selection."
    stop(m)
  } else if ((ncpus > 2)){
    m <- "Maximum number of CPUs should be 2."
    stop(m)
  }
}
check_nboots <- function(nboots){
  mt <- NULL
  if (missing(nboots)){
    nboots <- 500
    m <- "500 bootstrap samples will be computed."
    mt <- rbind(mt, c("ci.nboots", m))
  } else if(!(is.numeric(nboots) | nboots > 0)){
    m <- "Invalid ci.nboots selection."
    stop(m)
  }
  ret <- list()
  ret$nboots  <- nboots
  ret$message <- mt
  ret
}
check_ci_cl <- function(ci.cl){
  mt <- NULL
  if(missing(ci.cl)){
    ci.cl <- 0.95
    m <- "Confidence Intervals will be computed at 0.95 confidence level."
    mt <- rbind(mt, c("ci.cl",m))
  } else if(!(is.numeric(ci.cl))){
    m <- paste(ci.cl, "-", "Should be a numerical value between 0 and 1.")
    stop(m)
  } else if(!( (0 <= ci.cl) & (ci.cl <= 1))){
    m <- "Invalid confidence level."
    stop(m)
  }
  ret         <- list()
  ret$ci.cl   <- ci.cl
  ret$message <- mt
  ret
}
