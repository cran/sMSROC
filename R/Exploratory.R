explore_table <- function(marker = NULL,
                          status = NULL,
                          observed.time = NULL,
                          left = NULL,
                          right = NULL,
                          time = 1,
                          d = 2,
                          ...){
  M  <- data.frame()
  rn <- NULL
  data_type <- check_type_outcome(status, observed.time, left, right)
    if (data_type$type.outcome == "binout"){
      data <- check_marker_binout(marker, status)
    } else{
      data_time <- check_tim(time)$time
        if (data_type$type.outcome == "timerc"){
          data <- check_marker_timerc(marker, status, observed.time, data_time)
        } else if (data_type$type.outcome == "timeic"){
          data <- check_marker_timeic(marker, left, right, data_time)
        } else if (data_type$type.outcome == "unknow"){
          m <- "Not able to determine the distributions with the current input data."
          stop(m)
        }
      }
      if(!is.numeric(d)){
        stop ("Round decimals should be an integer number.")
      } else if (!(floor(d)==d)){
        stop ("Round decimals should be an integer number.")
      }
      pos <- data$marker[which(data$outcome == 1)]
      neg <- data$marker[which(data$outcome == 0)]
      mso <- data$marker[which(data$outcome == -1)]
      if (length(pos) > 0){
        M <- rbind(M, c(data$cases,
                        round(min(pos), d), round(max(pos), d),
                        round(mean(pos), d), round(sd(pos), d),
                        round(sd(pos) ^ 2, d),
                        round(quantile(pos, probs = c(0.25, 0.5, 0.75)), d)))
        rn <- rbind(rn, "Positive")
      }
      if (length(neg) > 0){
        M <- rbind(M, c(data$controls,
                        round(min(neg), d), round(max(neg), d),
                        round(mean(neg), d), round(sd(neg), d),
                        round(sd(neg) ^ 2, d),
                        round(quantile(neg, probs = c(0.25, 0.5, 0.75)), d)))
        rn <- rbind(rn, "Negative")
      }
      if (length(mso) > 0){
        M <- rbind(M, c(data$misout,
                        round(min(mso), d),
                        round(max(mso), d),
                        round(mean(mso), d),
                        round(sd(mso), d),
                        round(sd(mso) ^ 2, d),
                        round(quantile(mso, probs = c(0.25, 0.5, 0.75)), d)))
        rn <- rbind(rn, "Miss/Cens/Und")
      }
      if ((length(pos) | length(neg) | length(mso)) > 0 ){
        M <- rbind(M, c(length(data$marker),
                        round(min(data$marker), d),
                        round(max(data$marker), d),
                        round(mean(data$marker), d),
                        round(sd(data$marker), d),
                        round(sd(data$marker)^2, d),
                        round(quantile(data$marker, probs = c(0.25, 0.5, 0.75)),d)))
        rn <- rbind(rn, "Total")
      }
      M  <- cbind(rn, M)
      hd <- c("Sample", "Size", "Minimun", "Maximun", "Mean",
              "Sd", "Variance", "Q1", "Median", "Q3")
      colnames(M) <- hd
      ret <- list()
      ret$summary <- M
      ret$table   <- flextable(M, ...)
      ret
}
explore_plot <- function(marker = NULL,
                         status = NULL,
                         observed.time = NULL,
                         left = NULL,
                         right = NULL,
                         time = 1){
  outcome <- NULL
  data_type <- check_type_outcome(status, observed.time, left, right)
  if (data_type$type.outcome == "binout"){
    data <- check_marker_binout(marker, status)
  } else {
    data_time <- check_tim(time)
    if (data_type$type.outcome == "timerc"){
      data <- check_marker_timerc(marker, status, observed.time, data_time$time)
    } else if (data_type$type.outcome == "timeic"){
      data <- check_marker_timeic(marker, left, right, data_time$time)
    } else if (data_type$type.outcome == "unknow"){
      m <- "Not able to determine distributions with the current input data."
      stop(m)
    }
  }
  neg <- data$marker[which(data$outcome== 0)]
  pos <- data$marker[which(data$outcome == 1)]
  dt1 <- data.frame(marker = neg, outcome  = rep("0", length(neg)))
  dt2 <- data.frame(marker = pos, outcome  = rep("1", length(pos)))
  dt  <- rbind(dt1, dt2)
  p   <- ggplot2::ggplot(dt, aes(x = marker, color = outcome, fill = outcome)) +
           geom_density(alpha = 0.7, linewidth = 1.2) +
           scale_color_manual(values = c("#a2b285", "#5499C7" ),
                              labels = c("Negative", "Positive"),
                              name = "         ") +
           scale_fill_manual(values=  c("#a2b285", "#5499C7"),
                              labels = c("Negative", "Positive"),
                              name = "         ") +
           xlab("Marker") +
           ylab("Densities") +
           ggtitle("         ") +
           theme_classic() +
           theme(panel.border = element_blank(),
                 plot.title   = element_text(face="bold", size = 22, hjust = 0.5),
                 axis.text.x  = element_text( size=rel(2)),
                 axis.text.y  = element_text( size=rel(2)),
                 axis.title.x = element_text(face="bold", size = 22),
                 axis.title.y = element_text(face="bold", size = 22),
                 legend.text  = element_text(face="bold", size = 20),
                 legend.position = "top")
  ret <- list()
  ret$plot <- p
  ret$neg  <- neg
  ret$pos  <- pos
  ret
}

