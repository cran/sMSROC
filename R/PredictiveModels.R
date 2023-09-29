pred_model_binout <- function(marker, status, meth){
   Ics <- which(status >= 0)
   if (meth == "L"){
      mod <- glm(status[Ics] ~ marker[Ics], family=binomial)
   } else if (meth == "S"){
      mod <- glm(status[Ics] ~ rcs(marker[Ics]), family=binomial)
   }
   dt  <- cbind.data.frame(marker = marker[Ics], model = as.numeric(-predict(mod)))
   dtu <- unique(dt)
   dtu <- dtu[order(dtu$marker),]
   vmi <- min(marker)
   vma <- max(marker)
   if (length(which(dtu$marker == vmi)) > 0){
      pmi <- dtu$model[which(dtu$marker == vmi)]
   } else {
      I   <- which(dtu$marker > vmi)
      pmi <- dtu$model[I[1]]
   }
   if (length(which(dtu$marker == vma)) > 0){
      pma <- dtu$model[which(dtu$marker == vma)]
   } else {
      I   <- which(dtu$marker < vma)
      J   <- which(dtu$marker == max(dtu$marker[I]))
      pma <- dtu$model[J[1]]
   }
   fuT <- approxfun(c(vmi, dtu$marker, vma), c(pmi,dtu$model, pma))(marker)
   MP  <- cbind(marker, as.numeric(1 / (1 + exp(fuT))))
   MP  <- MP[order(marker),]
   ret <- list()
   ret$marker <- MP[,1]
   ret$probs  <- MP[,2]
   ret
}
pred_model_emp <- function(marker, status){
   Ics <- which(status >= 0)
   PR  <- status[Ics]
   MP  <- cbind(marker[Ics], status[Ics], PR)
   MP  <- na.omit(MP[order(marker[Ics]),])
   ret <- list()
   ret$marker  <- MP[,1]
   ret$outcome <- MP[,2]
   ret$probs   <- MP[,3]
   ret
}
pred_model_timerc <- function(marker, status, observed.time, outcome, time, meth){
   fuT <- function(i, basehtime, basehhazard, predtime, t){
      fuT <- approxfun(basehtime, 1 - exp(-basehhazard * exp(predtime[i])))(t)
   }
   if (meth == "L"){
      mod <- coxph(Surv(observed.time, status) ~ marker)
   } else {
      mod <- coxph(Surv(observed.time, status) ~ pspline(marker))
   }
   BP  <- predict(mod, type = "lp")
   BH  <- survfit(mod, newdata = data.frame(marker = 0))
   P   <- sapply(1:length(status), fuT, BH$time, BH$cumhaz, BP, time)
   MP  <- cbind(marker, P, outcome)
   MP  <- MP[order(marker),]
   ret <- list()
   ret$marker  <- MP[,1]
   ret$probs   <- MP[,2]
   ret$outcome <- MP[,3]
   ret
}
pred_model_timeic <- function(marker, left, right, outcome, time, meth){
   fuP <- function(point, Tbull.ints, Scurve){
      In <- which((Tbull.ints[,1] <= point) & (Tbull.ints[,2] >= point))
      if (length(In)==0){
         pi <- ifelse(length(which(Tbull.ints[,1] <= point)) > 0 ,
         Scurve[[1]][max(which((Tbull.ints[,1] <= point)))], 1)
      } else {
         if (Tbull.ints[In,1] == Tbull.ints[In,2]){
            pi <- Scurve[[1]][In]
         } else{
            pi <- approxfun(c(Tbull.ints[In,1], Tbull.ints[In,2]),
                            c(Scurve[[1]][In-1], Scurve[[1]][In]))(point)
         }
      }
   }
   fuA <- function(i, point, Tbull.ints, Scurve, left, right){
      vp <- ifelse(is.na(right[i]),0, fuP(right[i], Tbull.ints, Scurve[i]))
      up <- fuP(left[i], Tbull.ints, Scurve[i])
      tp <- fuP(point, Tbull.ints, Scurve[i])
      if (left[i] > point){
         P <-  1 - ((tp - vp)/(tp - vp))
      } else {
         if(is.na(right[i])){
            P <- 1 - ((tp - vp)/(up - vp))
         } else if (right[i] <= point){
            P <- 1 - ((tp - tp)/(up - tp))
         } else {
            P <- 1 - ((tp - vp)/(up - vp))
         }
      }
   }
   dat <- data.frame(marker = marker, left = left, right = right)
   mod <- ic_sp(Surv(left, right, type="interval2") ~ marker, data = dat, B = c(0,1), model = "ph")
   new <- data.frame(marker = c(as.numeric(marker)))
   row.names(new) <- sapply(1:length(marker), function(i){paste0("grp",i)})
   SurvCurves     <- getSCurves(mod,new)
   if (meth == "L"){
      PR <- sapply(1:length(marker), fuA, time, SurvCurves$Tbull_ints, SurvCurves$S_curves, left, right)
   } else{
      PR <- sapply(1:length(marker), function(i){1 - fuP(time, SurvCurves$Tbull_ints, SurvCurves$S_curves[i])})
   }
   MP  <- cbind(marker, outcome, PR)
   MP  <- MP[order(marker),]
   ret <- list()
   ret$marker  <- MP[,1]
   ret$outcome <- MP[,2]
   ret$probs   <- MP[,3]
   ret
}
