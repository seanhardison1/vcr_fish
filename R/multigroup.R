# This script implements a fix for the piecewiseSEM::multigroup function written
# by Jon Lefcheck. The fix was introduced by Elisa Lonso Aller here: 
# https://github.com/jslefche/piecewiseSEM/issues/209,  and I build on it to allow for multiple
# random intercept terms. I also included a step for simulating model residuals from SEM
# submodels that are fitted to subsets of data in the multigroup analysis process. 

multigroup <- 
  function (modelList, group, standardize = "scale", standardize.type = "latent.linear", 
          test.type = "III", model_sim = T) 
{
  name <- deparse(match.call()$modelList)
  data <- modelList$data
  modelList <- piecewiseSEM:::removeData(modelList, formulas = 1)
  intModelList <- lapply(modelList, function(i) {
    rhs2 <- ifelse(any(grepl("merMod", class(i)) == T), paste(paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group),
                                                               collapse = " + "),
                                                         paste("+ (", 
                                                               as.character(findbars(formula(i))),")", 
                                                               collapse = " ")),
                   paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group), collapse = " + "))
    message(rhs2)
    i <- update(i, formula(paste(". ~ ", rhs2)))
    # print(summary(i))
    if (model_sim) {s <- DHARMa::simulateResiduals(i,
                                                   n = 1000);DHARMa:::plot.DHARMa(s)}
    return(i)
  })
  newModelList <- lapply(unique(data[, group]), function(i) {
    if (model_sim){
      message(paste("Fitting model with data subsetted by",i))
      j <- update(as.psem(modelList), 
             data = data[data[, group] == i, ])
      for (l in 1:(length(j)-1)){
        s <- DHARMa::simulateResiduals(j[[l]],n = 1000)
        message(formula(j[[l]]))
        oldpar <- par(mfrow = c(1,2), oma = c(0,1,2,1))
        on.exit(par(oldpar))
        DHARMa:::plotQQunif(s)
        DHARMa:::plotResiduals(s)
      }
    }
    return(j)
    })
  
  names(newModelList) <- unique(data[, group])
  
  # extract coefficients fitted to data from different groups
  coefsList <- lapply(newModelList, coefs, standardize, standardize.type, 
                      test.type)
  names(coefsList) <- unique(data[, group])
  
  # extract coefficients from global model fit
  coefTable <- coefs(modelList, standardize, standardize.type, 
                     test.type)
  
  # use ANOVA to test for significance of interactions
  anovaTable <- anova(as.psem(intModelList))[[1]]
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), 
  ]
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Response", 
                                                   "Predictor")]
  global$Predictor <- sub(":", "\\1", sub(group, "\\1", global$Predictor))
  if (nrow(global) == nrow(anovaInts)) 
    newCoefsList <- list(global = coefTable)
  else {
    newCoefsList <- lapply(names(coefsList), function(i) {
      ct <- as.matrix(coefsList[[i]])
      idx <- which(apply(ct[, 1:2], 1, paste, collapse = "___") %in% 
                     apply(global[, 1:2], 1, paste, collapse = "___"))
      ct[idx, ] <- as.matrix(coefTable[idx, ])
      ct <- cbind(ct, ifelse(1:nrow(ct) %in% idx, "c", 
                             ""))
      for (j in 1:nrow(ct)) {
        if (ct[j, ncol(ct)] == "c") {
          model <- modelList[[which(sapply(piecewiseSEM:::listFormula(modelList), 
                                           function(x) piecewiseSEM:::all.vars.merMod(x)[1] == ct[j, 
                                                                                   "Response"]))]]
          data. <- data[data[, group] == i, ]
          sd.x <- piecewiseSEM:::GetSDx(model, modelList, data., standardize)
          sd.x <- sd.x[which(names(sd.x) == ct[j, "Predictor"])]
          sd.y <- piecewiseSEM:::GetSDy(model, data., standardize, standardize.type)
          new.coef <- as.numeric(ct[j, "Estimate"]) * 
            (sd.x/sd.y)
          ct[j, "Std.Estimate"] <- ifelse(length(new.coef) > 
                                            0, round(as.numeric(new.coef), 4), "-")
        }
      }
      ct <- as.data.frame(ct)
      ct[is.na(ct)] <- "-"
      names(ct)[(ncol(ct) - 1):ncol(ct)] <- ""
      return(ct)
    })
    names(newCoefsList) <- names(coefsList)
  }
  if (nrow(global) == nrow(anovaInts)) 
    gof <- fisherC(modelList)
  else {
    b <- basisSet(modelList)
    cf <- coefTable[coefTable$Response %in% global$Response & 
                      coefTable$Predictor %in% global$Predictor, ]
    b <- lapply(b, function(i) {
      for (j in 3:length(i)) {
        value <- cf[cf$Response == i[2] & cf$Predictor == 
                      i[j], "Estimate"]
        if (length(value) != 0) 
          i[j] <- paste0("offset(", value, "*", i[j], 
                         ")")
      }
      return(i)
    })
    if (length(b) == 0) 
      b <- NULL
    gof <- fisherC(modelList, basis.set = b)
  }
  ret <- list(name = name, group = group, global = global, 
              anovaInts = anovaInts, group.coefs = newCoefsList, Cstat = gof)
  class(ret) <- "multigroup.psem"
  return(ret)
}
