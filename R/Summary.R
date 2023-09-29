probs_pred <- function (sMS,
                        var =  c("F", "T"),
                        nboots = 500,
                        parallel = c("F", "T"),
                        ncpus = 1){
  data_var      <- match.arg(var)
  data_parallel <- match.arg(parallel)
  sd.matr <- NULL
  x <- NULL
  y <- NULL
  if (data_parallel == "T"){
    if(!(is.integer(ncpus) | ncpus > 0)){
      m <- "Invalid ncpus selection."
      stop(m)
    } else if (ncpus > 2){
      m <- "Maximum number of CPUs should be 2."
      stop(m)
    }
  }
  if (!(inherits(sMS, "sMSROC"))){
    m <- paste("Object", sMS, "should be 'sMSROC' object")
    stop(m)
  }
  df <- data.frame(x = sMS$thres, y = sMS$probs)
  p  <- ggplot2::ggplot(data = df, aes(x, y)) +
          geom_point(size = 1.5, alpha = 0.6, colour = "#003366") +
          xlab("Marker") +
          ylab("Probabilities") +
          ggtitle("Predictive Model") +
          theme_classic() +
          theme(panel.border = element_blank(),
                plot.title   = element_text(face="bold", size = 15, hjust = 0.5),
                axis.text.x  = element_text( size=rel(1.5)),
                axis.text.y  = element_text( size=rel(1.5)),
                axis.title.x = element_text(face="bold",size = 12),
                axis.title.y = element_text(face="bold",size = 12))
  if (data_var == "T"){
    var <- variance_probs(marker = sMS$data$marker,
                          outcome = sMS$data$outcome,
                          status = sMS$data$status,
                          observed.time = sMS$data$observed.time,
                          left = sMS$data$left,
                          right = sMS$data$right,
                          meth = sMS$data$meth,
                          data_type = sMS$data$type,
                          time = sMS$time,
                          parallel = data_parallel,
                          ci.nboots = nboots,
                          grid = sMS$data$grid,
                          ncpus = ncpus)
    sd.matr    <- splinefun(sort(sMS$data$marker), var$sd.probs)(sort(sMS$data$marker))
    ymin <- (sMS$probs - 1.96 * sd.matr)
    ymax <-  pmin(sMS$probs + 1.96 * sd.matr, rep(1,length(sMS$thres)))
    df1  <- data.frame(x = sMS$thres, y = sMS$probs, ymin = ymin, ymax = ymax)
    p    <- ggplot2::ggplot(data = df1, aes(x, y)) +
              geom_ribbon(data = df1,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "#f0f8ff"),
              alpha = 0.8, size =1.2, color = "#f0f8ff", show.legend = FALSE ) +
              scale_fill_manual(values = "#f0f8ff") +
              geom_line(size = 1.5, colour = "#003366")+
              xlab("Marker") +
              ylab("Probability") +
              ggtitle("Predictive Model") +
              theme_classic() +
              theme(panel.border = element_blank(),
                    plot.title   = element_text(face="bold", size = 15, hjust = 0.5),
                    axis.text.x  = element_text( size=rel(1.5)),
                    axis.text.y  = element_text( size=rel(1.5)),
                    axis.title.x = element_text(face="bold",size = 12),
                    axis.title.y = element_text(face="bold",size = 12))
  }
  ret <- list()
  ret$plot     <- p
  ret$thres    <- sMS$thres
  ret$probs    <- sMS$pros
  ret$sd.probs <- sd.matr
  ret
}
evol_auc <- function(marker,
                     status,
                     observed.time,
                     left,
                     right,
                     time = 1,
                     meth = c("L", "S", "E"),
                     grid = 500){
  e.time  <- time
  e.auc <- NULL
  if (missing(status)){status <- NULL}
  if (missing(observed.time)){observed.time <- NULL}
  if (missing(left)){left <- NULL}
  if (missing(right)){right <- NULL}
  l <- length(e.time)
  if (l == 1){
    m <- "A single point of time was indicated"
    warning(m)
  }
  e.time <- sapply(1:l, function(i){check_tim(e.time[i])$time})
  e.auc  <- sapply(1:l, function(i){sMSROC(marker = marker,
                                           status = status,
                                           observed.time = observed.time,
                                           left = left,
                                           right = right,
                                           time = e.time[i],
                                           meth = meth,
                                           grid = grid)$auc})
  df     <- data.frame(time = e.time, e.auc = e.auc)

  plot   <- ggplot2::ggplot(data = df, aes(time, e.auc)) +
              geom_line(linewidth = 2.5, alpha = 0.6, colour = "#a2b285") +
              xlab("Follow-up ") +
              ylab("AUC") +
              ggtitle("Evolution of the AUC") +
              theme_classic() +
              theme(panel.border = element_blank(),
              plot.title  = element_text(face="bold", size = 15, hjust = 0.5),
              axis.text.x = element_text( size=rel(1.5)),
              axis.text.y = element_text( size=rel(1.5)),
              axis.title.x = element_text(face="bold",size = 12),
              axis.title.y = element_text(face="bold",size = 12))
  ret <- list()
  ret$evol.auc  <- plot
  ret$time <- time
  ret$auc <- e.auc
  ret
}
sMSROC_plot <- function(sMS, m.value){
  dp <- function(x){
    if ((x %% 1) != 0){
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  if (!(inherits(sMS, "sMSROC"))){
    m <- paste("Object", sMS, "should be 'sMSROC' object")
    stop(m)
  }
  SP  <- NULL; SP  <- sMS$SP
  TPR <- NULL; TPR <- c(0, sMS$SE, 1)
  FPR <- NULL; FPR <- c(0, 1-SP, 1)
  thr <- NULL; thr <- c(sMS$thres[length(sMS$thres)], sMS$thres, sMS$thres[1])
  x   <- NULL; y   <- NULL;  m   <- NULL
  df  <- unique(data.frame(thr = thr, FPR = FPR, TPR = TPR))
  df  <- df[order(df$FPR),]
  df1 <- data.frame(x = c(0,1), y = c(0,1))
  p.basic  <- ggplot2::ggplot(df, aes(FPR, TPR, label = thr)) +
                geom_roc(stat = "identity") +
                style_roc()
  roc.plot <- ggplot2::ggplot(df, aes(FPR, TPR)) +
                geom_line(linewidth = 2.5, alpha = 0.6, colour = "#003366") +
                xlab("1 - Specificity") +
                ylab("Sensitivity") +
                ggtitle("ROC Curve") +
                theme_classic() +
                theme(panel.border = element_blank(),
                plot.title  = element_text(face="bold", size = 15, hjust = 0.5),
                axis.text.x = element_text( size=rel(1.5)),
                axis.text.y = element_text( size=rel(1.5)),
                axis.title.x = element_text(face="bold",size = 12),
                axis.title.y = element_text(face="bold",size = 12)) +
                scale_x_continuous(limits = c(0,1),
                                   breaks = seq(0,1,0.2)) +
                scale_y_continuous(limits = c(0,1),
                                   breaks = seq(0,1,0.2)) +
                geom_line(data = df1, aes(x, y), size = 1.5, colour = "gray", linetype = "twodash")
  ret <- list()
  ret$basicplot <- p.basic
  ret$rocplot   <- roc.plot
  ct <- NULL
  if (! missing(m.value)){
    if(!(is.numeric(m.value))){
      m <- "Non-valid marker value"
    } else {
      d <- dp(m.value)
      df$thr <- round(df$thr, digits = d)
      m.pos  <- which(df$thr == m.value)
      if (length(m.pos) > 0){
        ct <- cbind(df$FPR[m.pos], df$TPR[m.pos], df$thr[m.pos])
      } else {
        if (m.value >= max(df$thr)){
          ct <- cbind(df$FPR[length(df$thr)],
                      df$TPR[length(df$thr)],
                      df$thr[length(df$thr)])
        } else {
          if (m.value <= min(df$thr)){
            ct <- cbind(df$FPR[1], df$TPR[1], df$thr[1])
          } else {
            m.pos <- which(df$thr > m.value)
            ct <- rbind(ct, cbind(df$FPR[m.pos[1]], df$TPR[m.pos[1]], df$thr[m.pos[1]]))
            ct <- rbind(ct, cbind(df$FPR[m.pos[1] - 1], df$TPR[m.pos[1] - 1], df$thr[m.pos[1] - 1]))
          }
        }
      }
    }
  }
  if (!is.null(ct)){
    if (dim(ct)[1] == 1){
      df2 <- data.frame(x = ct[,1], y = ct[,2], m = ct[,3])
      roc.plotp <- roc.plot +
                   geom_point(data = df2, aes(x, y, shape = "18"), size = 5, show.legend = FALSE) +
                   geom_label(data = df2,
                              aes(x, y, label = (round(m,2))),
                              size = 5,
                              vjust = 1,
                              hjust = -0.5,
                              label.size = 0)
    } else {
      df2 <- data.frame(x = ct[1,1], y = ct[1,2], m = ct[1,3])
      df3 <- data.frame(x = ct[2,1], y = ct[2,2], m = ct[2,3])
      roc.plotp <- roc.plot +
                   geom_point(data = df2,
                              aes(x, y, shape = "18"),
                              size = 5,
                              show.legend = FALSE) +
                   geom_label(data = df2,
                              aes(x, y, label = (round(m,2))),
                              size = 5,
                              vjust = 1,
                              hjust = -0.5,
                              label.size = 0) +
                   geom_point(data = df3,
                              aes(x, y,  shape = "18"),
                              size = 5,
                              show.legend = FALSE) +
                   geom_label(data = df3,
                              aes(x, y, label = (round(m,2))),
                              size = 5,
                              vjust = 1,
                              hjust= -0.5,
                              label.size = 0)
      }
      ret$rocplot <-roc.plotp
  }
  ret
}
conf_int_print <- function(sMS){
  if (!(inherits(sMS, "sMSROC"))){
    m <- paste("Object", sMS, "should be 'sMSROC' object")
    stop(m)
  }
  print(paste0("AUC: ",
               round(sMS$auc,2),
               "; ",
               sMS$ci.cl,
               "% C.I.", "[",
               round(sMS$auc.ci.l,2),
               ", ",
               round(sMS$auc.ci.u,2),
               "]"))
}
print.sMSROC <- function(x, ...){
  if (!(inherits(x, "sMSROC"))){
    m <- 'Object should have class sMSROC'
    stop(m)
  }
  pos <- length(which(x$data$outcome == 1))
  neg <- length(which(x$data$outcome == 0))
  und <- length(which(x$data$outcome == -1))
  if (x$data$meth == 'M'){
    pos <- length(which(x$probs == 1))
    neg <- length(which(x$probs == 0))
    und <- length(which(x$probs > 0 & x$probs < 1))
  }
  cat(paste0("The AUC is ", round(x$auc,3), '.'), "\n")
  if (x$data$type$type.outcome == 'binout'){
    if (x$data$meth == 'L'){
      cat(paste0("Predictive model computed through a logistic regression model,"))
    } else if (x$data$meth == 'S'){
      cat(paste0("Predictive model computed through a smooth logistic regression model,"))
    } else if (x$data$meth == 'M'){
      cat(paste0("Predictive model externally computed,"))
    } else if (x$data$meth == 'E'){
      cat(paste0("Predictive model computed removing participants whose status in unknonw,"))
    }
  } else if (x$data$type$type.outcome == 'timerc') {
      if (x$data$meth == 'L'){
        cat(paste0("Predictive model computed through a Cox PH regression model,"))
      } else if (x$data$meth == 'S'){
        cat(paste0("Predictive model computed through a smooth Cox PH regression model,"))
      } else if (x$data$meth == 'M'){
        cat(paste0("Predictive model externally computed,"))
      } else if (x$data$meth == 'E'){
        cat(paste0("Predictive model computed removing participants whose status in unknonw,"))
      }
  } else if (x$data$type$type.outcome == 'timeic') {
      if (x$data$meth == 'L'){
        cat(paste0("Predictive model computed through a D.Finkelstein PH regression model, accounting to the length of the observed intervals,"))
      } else if (x$data$meth == 'S'){
        cat(paste0("Predictive model computed through a D.Finkelstein PH regression model, not taking into account the length of the observed intervals,"))
      } else if (x$data$meth == 'M'){
        cat(paste0("Predictive model externally computed," ))
      } else if (x$data$meth == 'E'){
        cat(paste0("Predictive model computed removing participants whose status in unknonw,"))
      }
  } else {cat("Predictive model not computed.")}
  cat(" based on", pos, "positive,", neg, "negative and", und, "undefined (mixed) subjects.", "\n")
}




