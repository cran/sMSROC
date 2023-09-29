compute_ROC <- function(marker, probs, grid){
  SE <- sapply(1:length(marker),
  function(i){sum(probs[which(marker > marker[i])])/sum(probs)})
  SP <- sapply(1:length(marker),
  function(i){sum(1 - probs[which(marker <= marker[i])])/sum(1 - probs)})
  u   <- seq(0, 1, 1/grid)
  dt  <- cbind.data.frame(marker = marker, SE = SE, SP = SP)
  dtu <- unique(dt)
  ROC <- approxfun(c(0, 1-dtu$SP, 1), c(0, dtu$SE, 1), ties = max)(u)
  auc <- sum(ROC) * 1/grid
  ret      <- list()
  ret$SE   <- SE
  ret$SP   <- SP
  ret$u    <- u
  ret$ROC  <- ROC
  ret$auc  <- auc
  ret$marker <- marker
  ret$probs  <- probs
  ret
}
sMS_binout <- function(marker, status, meth, grid, probs, all){
  if (meth == "E"){
    p.model <- pred_model_emp(marker, status)
  } else {
    if (meth == "M"){
      MP  <- data.frame(marker = marker, probs = probs)
      MP  <- MP[order(MP$marker),]
      p.model <- list()
      p.model$marker <- MP$marker
      p.model$probs  <- MP$probs
    } else {
      p.model <- pred_model_binout(marker, status, meth)
    }
  }
  if (all == "F"){
    Index <- which(p.model$status >= 0)
    p.model$probs[Index] == p.model$status
  }
  bin.ROC  <- compute_ROC(p.model$marker, p.model$probs, grid)
    ret      <- list()
    ret$SE   <- bin.ROC$SE
    ret$SP   <- bin.ROC$SP
    ret$u    <- bin.ROC$u
    ret$ROC  <- bin.ROC$ROC
    ret$auc  <- bin.ROC$auc
    ret$marker <- p.model$marker
    ret$probs  <- p.model$probs
    ret
}
sMS_timerc <- function(marker,
                       status,
                       observed.time,
                       outcome,
                       time,
                       meth,
                       grid,
                       probs,
                       all){
  if (meth == "E"){
    Ims     <- which(outcome > -1)
    p.model <- pred_model_emp(marker[Ims], outcome[Ims])
  } else {
    if (meth == "M"){
      MP  <- data.frame(marker = marker, probs = probs)
      MP  <- MP[order(MP$marker),]
      p.model <- list()
      p.model$marker <- MP$marker
      p.model$probs  <- MP$probs
    } else {
      p.model <- pred_model_timerc(marker, status, observed.time, outcome, time, meth)
    }
  }
  if (all == "F"){
    Index <- which(p.model$outcome >= 0)
    p.model$probs[Index] <- p.model$outcome[Index]
  }
  timrc.ROC  <- compute_ROC(p.model$marker, p.model$probs, grid)
  ret      <- list()
  ret$SE   <- timrc.ROC$SE
  ret$SP   <- timrc.ROC$SP
  ret$u    <- timrc.ROC$u
  ret$ROC  <- timrc.ROC$ROC
  ret$auc  <- timrc.ROC$auc
  ret$marker  <- p.model$marker
  ret$outcome <- p.model$outcome
  ret$probs   <- p.model$probs
  ret
}
sMS_timeic <- function(marker, left, right, outcome, time, meth, grid, probs, all){
  if (meth == "E"){
    p.model <- pred_model_emp(marker, outcome)
  } else {
    if (meth == "M"){
      MP  <- data.frame(marker = marker, probs = probs)
      MP  <- MP[order(MP$marker),]
      p.model <- list()
      p.model$marker <- MP$marker
      p.model$probs  <- MP$probs
    } else {
      p.model <- pred_model_timeic(marker, left, right, outcome, time, meth)
    }
  }
  if (all == "F"){
    Index <- which(p.model$outcome >= 0)
    p.model$probs[Index] == p.model$outcome
  }
  timeic.ROC <- compute_ROC(p.model$marker, p.model$probs, grid)
  ret      <- list()
  ret$SE   <- timeic.ROC$SE
  ret$SP   <- timeic.ROC$SP
  ret$u    <- timeic.ROC$u
  ret$ROC  <- timeic.ROC$ROC
  ret$auc  <- timeic.ROC$auc
  ret$marker <- p.model$marker
  ret$oucome <- p.model$outcome
  ret$probs  <- p.model$probs
  ret
}
