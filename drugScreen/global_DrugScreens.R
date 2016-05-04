library("nplr")

get_drugResponse_stats <- function(conc,viability,...){
  res <- nplr(conc, viability,...)
  results <- getAUC(res)
  results['goodNess_of_fit'] <- getGoodness(res)
  results['stdErr'] <- getStdErr(res)
  ICx_est = getEstimates(res, targets= c(.10,.20,.30,.40,.50,.60,.70,.80,.90))
  results['IC10'] = ICx_est[1,'x']
  results['IC20'] = ICx_est[2,'x']
  results['IC30'] = ICx_est[3,'x']
  results['IC40'] = ICx_est[4,'x']
  results['IC50'] = ICx_est[5,'x']
  results['IC60'] = ICx_est[6,'x']
  results['IC70'] = ICx_est[7,'x']
  results['IC80'] = ICx_est[8,'x']
  results['IC90'] = ICx_est[9,'x']
  results['maxEfficacy'] = max(getYcurve(res)) #get the maximum efficacy of the drug
  results['bottom_asymptote'] = res@pars['bottom']
  results['top_asymptote'] = res@pars['top']
  results['hillSlope'] =  res@pars['scal']
  
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
    print(unique(df$cellLine))
    print(unique(df$drug))
    print(unique(df$experiment))
    stop('stopped')
  })
}
