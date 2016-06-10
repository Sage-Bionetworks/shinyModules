library("nplr")
#require(memoise)

get_drugResponse_stats <- function(conc,viability,...){
  #res <- memoise(function(conc, viability,...) nplr(conc, viability,...))
  res <- nplr(conc, viability,...)
  results <- getAUC(res)
  #results['goodNess_of_fit'] <- getGoodness(res)
  #results['stdErr'] <- getStdErr(res)
  ICx_est = getEstimates(res, targets= c(.10,.20,.30,.40,.50,.60,.70,.80,.90))
  results['IC50'] = ICx_est[5,'x']
  results['maxEfficacy'] = max(getYcurve(res)) #get the maximum efficacy of the drug
  #results['bottom_asymptote'] = res@pars['bottom']
  #results['top_asymptote'] = res@pars['top']
  #results['hillSlope'] =  res@pars['scal']
  
  fittedVals <- data.frame(fittedX = getXcurve(res),
                           fittedY = getYcurve(res))
  results <- cbind(results,fittedVals)
  results
}


tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$conc, df$normViability, useLog=T)  
  },error=function(e){
    print(dim(df))
    print(df$conc)
    print(df$normViability)
    print(unique(df$sample))
    print(unique(df$drug))
    stop('stopped')
  })
}

QC_plot_name <- function(data){
  result <- "maxResp"
  if(!all(is.na(data$AUC))){
    result <- c(result, "AUC")
  }
  if(!all(is.na(data$AC50))){
    result <- c(result, "AC50")
  }
  return(result)
}
