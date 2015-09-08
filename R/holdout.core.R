##
## Core functions for a holdout calculations
##
##
## Author: Tomasz Konopka
##




#######################################################################
## Generic helper functions


## functions to compute correlations
## pearsonR - computes pearson correlation the proper way (no tricks)
pearsonR = function(x, y) {
  xres = x-mean(x)
  yres = y-mean(y)
  TA = sum(xres*yres)
  TB = sqrt(sum(xres*xres)*sum(yres*yres))
  return(TA/TB)
}
## pearsonRfast - computes pearson correlation, but takes shortcuts
## assumes xres and yres both have mean zero.
## Also assumes the variances are one.
pearsonRfast = function(xres, yres) {
  TA = sum(xres*yres)
  return(TA/length(xres))
}

spearmanRho = function (x, y) {
  return(pearsonR(rank(x), rank(y)))
}



#######################################################################
## Functions for building classifiers

##' Compute correlations between a data matrix and a signal vector
##'
##' Evaluates correlations between features in a data matrix and a
##' signal vector. When only one data and signal object is provided,
##' the output is a vector of straightforward correlations. When a
##' background dataset and background signal vector are also provided,
##' the function treats the primary data as a reusable holdout: it
##' will output correlations from either the holdout data or the
##' background set.
##' 
##' @param dat - a data matrix with S samples in columns and F features
##' rows
##' @param signal - a numeric vector. The function will compute
##' correlations between rows in the data matrix and this signal vector.
##' @param dat.bg - background dataset
##' @param signal.bg - background signal
##' @param min.cor - minimal expected correlation (function
##' will output zero if the actuall correlation is below
##' threshold)
##' @param tolerance.factor - one of the penalties used in the
##' reusable holdout proposal
##' @param threshold.factor - one of the penalties used in the
##' reusable holdout proposal. 
##' 
##' @export
getFeatureCorrelations = function(dat, signal,
  min.cor=1/sqrt(length(signal)),
  dat.bg=NULL, signal.bg=NULL,
  tolerance.factor=1, threshold.factor=4) {
  
  if (length(signal) != ncol(dat)) {
    stop(paste0("data matrix and signal vector don't match\n",
                "(length(signal) != ncol(dat))\n"));
  }
  if (nrow(dat)<1) {
    return(c());
  }
  
  ## compute the correlations on main data and background data
  AAfindCors = function(dd, pp) {
    if (!is.null(dd) & !is.null(pp)) {
      return(apply(dd, 1, function(x) { pearsonR(x, pp) }))
    }    
    return(NULL)
  }
  dat.cors = AAfindCors(dat, signal)
  names(dat.cors) = rownames(dat)

  ## use background set if available
  if (!is.null(dat.bg)) {
    siglen = length(signal);
    sqsl = 1/sqrt(siglen);
    
    dat.bg.cors = AAfindCors(dat.bg, signal.bg)

    diff.cors = abs(dat.bg.cors-dat.cors)
    more.cors = diff.cors > (threshold.factor*sqsl) +
      rnorm(siglen, 0, tolerance.factor*sqsl);    
    nummore = sum(more.cors)

    ## features which are less correlated get bg values
    if (siglen-nummore>0) {
      dat.cors[!more.cors] = dat.bg.cors[!more.cors];
    }
    ## features which are more correlated get dataset values
    ## with some noise
    if (nummore>0) {
      dat.cors[more.cors] = dat.cors[more.cors] +
        rnorm(nummore, 0, tolerance.factor*sqsl)
    }        
  }

  ## only return the correlations above some threshold
  dat.cors[abs(dat.cors)<min.cor] = 0
  
  return(dat.cors)  
}


##' Attribute weights to each feature according to correlations values
##'
##' Function outputs a vector with values in {-1, 0, 1} determining
##' binary weights for each element in the cors vector.
##' Technically, the weights are not binary, bu
##' 
##' @param cors - vector or matrix with correlations
##' If it is a vector, it is intepreted as containign correlations for
##' each feature. If a matrix, it is intepreted as having features as rows,
##' and various estimates in columns. The function will use the mean
##' consensus correlation among the columns. 
##' @param k - number of features to call
##
##' @export
getBinaryFeatureWeights = function(cors, k) {
  
  ## if cors is a matrix, get a consensus correlation value using the mean
  if (class(cors)=="matrix") {
    numcols = ncol(cors)
    isbad = apply(sign(cors), 1, function(x) {abs(sum(x))<numcols} )
    cors = apply(cors, 1, mean)
    cors[isbad] = 0;
  }
  
  ## get the top k features
  topk = order(abs(cors), decreasing=T)[1:k]
  ans = rep(0, length(cors))
  ans[topk] = sign(cors[topk])
  names(ans) = names(cors);
  
  return(ans)  
}



##' Simple prediction of a +1,-1,0 class based on a data matrix
##' and feature weights
##' 
##' @param dat - a data matrix with samples in columns and features
##' in rows
##' @param weights - a numeric vector with weights for each feature.
##' 
##' @export
predictFromWeights = function(dat, weights) {  
  if (nrow(dat)!=length(weights)) {
    stop("number of rows must match weights vector")
  }
  ans = dat*weights;
  return(sign(apply(ans, 2, sum)))
}



##' Compute the accuracy comparing a set of truth values and predictions
##'
##' Compare a vector of truth classes and prediction classes.
##' The function is actually generic. It works for +1/-1/0 classes,
##' by any other labels work as well. 
##'
##' @param truth - a vector holding true labels
##' @param prediction - a vector holding labels output by some
##' prediction method
##' 
##' @export
evaluatePredictionAccuracy = function(truth, prediction) {
  if (length(truth)!=length(prediction)) {
    stop("inputs do not match")
  }
  if (length(truth)==0) {
    return(NA);
  }  
  return(sum(truth==prediction)/length(truth))  
}




###############################################################
## Functions for creating synthetic datasets


## two helper functions used in createRandomSet
createOneRandom = function(num.features, num.samples) {
  ans = matrix(rnorm(num.features*num.samples),
    ncol=num.samples, nrow=num.features)
  rownames(ans) = paste0("F", 1:num.features)
  return(ans)
}
createOneRandomWithSignal = function(num.features, num.samples,
  signal=sign(rnorm(num.samples)),
  nbiased=20, bias=6/sqrt(num.features)) {  
  ans = createOneRandom(num.features, num.samples)
  bb = bias*signal;
  for (i in 1:nbiased) {    
    ans[i,] = ans[i,]+bb;    
  }  
  return(ans)
}



##' Create a set of historical data and test data.
##'
##' This function will output a list with four objects.
##' Two objects will be data matrices representing historical
##' and test data. Two objects will be signals or sample classes
##' that describe the samples in both datasets.
##' 
##' @param num.features - number of features in each dataset
##' @param num.samples - number of samples in each dataset
##' @param nbiased - number of features with a true correlation
##' between feature value and signal value
##'
##' @param bias - a factor that determines the strength of
##' correlations for the nbiased features
##'
##' @export
createRandomSet = function(num.features, num.samples,
  nbiased=20, bias=6/sqrt(num.samples)) {
  
  train.signal = sign(rnorm(num.samples))
  test.signal = sign(rnorm(num.samples))
  ans = list(
    train = createOneRandomWithSignal(num.features, num.samples,
      nbiased, bias,
      signal=train.signal),
    test = createOneRandomWithSignal(num.features, num.samples,
      nbiased, bias,
      signal=test.signal),
    train.signal=train.signal,
    test.signal=test.signal)
  return(ans)  
}




