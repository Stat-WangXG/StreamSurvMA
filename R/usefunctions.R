Stx <- function(t,new11,w,time,status){   ####the survival function of weighted covariates
  Phiw <- as.matrix(new11)%*%w
  H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=t))
  Stx <- exp( -H0t*exp(Phiw) )
  return(Stx)
}

scorefunc <- function( batchData, beta, dims){      ###score function

  batchData <- batchData[order(-batchData$time), ]

  scores <- apply(batchData[, 1:(dims)], 1,
                  function(element) exp(element %*% beta) )

  nominator <- apply(batchData[,1:dims], 2,
                     function(element) cumsum(scores*element) )

  denominator <- cumsum(scores)

  partial_sum <- (batchData[,1:dims] - nominator/denominator)*batchData[, "status"]

  U_batch <- colSums(partial_sum)
  return(U_batch)

}


transcox_loglik <- function(mu,new11,time,status){  #### the log likelihood function of the parameter "mu".

  mu_1 <- append(x = mu, 0, after = length(mu))
  # mu_1 <- c(mu, 0)
  # mu_2 <- exp(mu_1)
  # mu_3 <- sum(mu_2)
  # MU<-mu_2/mu_3
  MU<-exp(mu_1)/sum(exp(mu_1))
  p1=dim(new11)[1]
  p2=dim(new11)[2]
  g<-c(2:(p2+1))
  Phiw <- as.matrix(new11)%*%MU
  Lw <- sum( status*(Phiw - sapply(time,function(timei){log(sum((time>=timei)*exp(Phiw)))}) ) )
  weight_trans <- Lw - (1/2)*log(p1)*as.vector(t(MU)%*%g)
  return(weight_trans)
}


Demu<-function(mu){

  mu_1 <- append(x = mu, 0, after = length(mu))

  num<-length(mu_1)
  for(ip in 1:num){
    assign(paste("M",ip,sep=""),mu_1[ip])
  }

  Ms <- noquote(paste0("M", 1:num))


  DMU <- matrix(0,nrow=num,ncol=num-1)
  HeMU <- matrix(0,nrow=num-1,ncol=num*(num-1))

  for(i in 1:num){
    Mi <- noquote(paste0("M", i))
    eMi<-noquote(paste0("exp(", Mi,")"))

    f=(noquote(paste(noquote(paste0(eMi,"/(", paste0("exp(", Ms, ")", collapse = "+"),")")),sep = "")))
    df <- eval(stats::deriv3(parse(text = f), paste("M", 1:(num-1), sep = "")))

    redf =attr(df,'gradient')
    Hedf =attr(df,'hessian')[1,,]

    DMU[i,]<-redf
    HeMU[,((num-1)*(i-1)+1):((num-1)*i)]<-Hedf

  }

  Dre<-list(DMU,HeMU)

  return(Dre)
}

bczph.cuee <-
  function(formula,
           cumA = NA,
           cumU = NA,
           prevBeta = NA,
           prevInvVar = NA,
           prevCumScore = NA,
           prevCumProd = NA,
           newdata,
           transform = "km")
  {
    ## local fit to get betahat
    newFit <-
      do.call(survival::coxph, list(formula = formula, data = newdata))
    ########################################################################
    ## scaled schoenfeld residual of the local fit, not very useful
    ## just need event times, betahat and its inverse variance
    ## schoenfeld residual
    varnames <- names(newFit$coefficients)
    nvar <- length(varnames)
    ssresid <- stats::resid(newFit, "schoenfeld")
    sresid <- ssresid %*% newFit$var
    ndead <- length(ssresid) / nvar
    times <- as.numeric(rownames(ssresid))
    if (is.character(transform)) {
      tname <- transform
      ttimes <- switch(
        transform,
        'identity' = times,
        'rank'    = rank(times),
        'log'     = log(times),
        'km' = {
          temp <- survival::survfitKM(factor(rep(1, nrow(
            newFit$y
          ))), newFit$y,
          se.fit = FALSE)
          # A nuisance to do left cont KM
          t1 <- temp$surv[temp$n.event > 0]
          t2 <- temp$n.event[temp$n.event > 0]
          km <- rep(c(1, t1), c(t2, 0))
          if (is.null(attr(sresid, 'strata')))
            1 - km
          else
            (1 - km[sort.list(sort.list(times))])
        },
        stop("Unrecognized transform")
      )
    } else {
      tname <- deparse(substitute(transform))
      if (length(tname) > 1)
        tname <- 'user'
      ttimes <- transform(times)
    }
    xx <- ttimes - mean(ttimes)
    ########################################################################
    betaHat <- newFit$coefficients
    invVarHat <- solve(newFit$var)
    ## if this is the first block, initialize all components to zero
    if (sum(is.na(prevBeta))) {
      cumA <- matrix(0, nvar, nvar)
      cumU <- matrix(0, nvar, 1)
      prevBeta <- matrix(0, nvar, 1)
      prevInvVar <- matrix(0, nvar, nvar)
      prevCumScore <- matrix(0, nvar, 1)
      prevCumProd <- matrix(0, nvar, 1)
    }
    ## get betaCheck, equation 21 in technometrics paper
    betaCheck <- solve(prevInvVar + invVarHat,
                       prevCumProd + invVarHat %*% betaHat)

    ## get the score vector evaluated at betaCheck
    ctrl <- survival::coxph.control(iter.max = 0)
    checkFit <- do.call(survival::coxph, list(formula = formula, data = newdata,
                                              init = betaCheck, control = ctrl))
    score <- prevCumScore + colSums(survival::coxph.detail(checkFit)$score)
    invVarCheck <- solve(checkFit$var)

    ## betaTilde
    betaTilde <- solve(prevInvVar + invVarCheck,
                       prevCumProd + invVarCheck %*% betaCheck + score)

    ## tildeFit to get the A and U parts of the cumulative test statistic
    tildeFit <- do.call(survival::coxph, list(formula = formula, data = newdata,
                                              init = betaTilde, control = ctrl))
    V <- tildeFit$var
    cA <- sum(xx ^ 2) * V / ndead
    nssresid <- resid(tildeFit, "schoenfeld")
    ## schoenfeld residual
    nsresid <- nssresid %*% V
    cU <- c(xx %*% nsresid)
    cumA <- cumA + cA
    cumU <- cumU + cU
    cumZ <- c(t(cumU) %*% solve(cumA, cumU))
    output <- list(
      cumA = cumA,
      cumU = cumU,
      beta = betaTilde,
      invVar = prevInvVar + invVarCheck,
      cumScore = score,
      cumProd = prevCumProd + invVarCheck %*% betaCheck,
      Stat = cumZ,
      p = 1 - pchisq(cumZ, ncol(sresid))
    )
    class(output) <- "bczph.cuee"
    return(output)
  }




bczph.cee <-
  function(formula,
           cumA = 0,
           cumU = 0,
           prevBeta = 0,
           prevInvVar = 0,
           newdata,
           transform = "km")
  {
    newfit <- do.call(survival::coxph, list(formula = formula, data = newdata))
    ## scaled schoenfeld residual
    ssresid <- stats::resid(newfit, "schoenfeld")
    ## schoenfeld residual
    sresid <- ssresid %*% newfit$var
    varnames <- names(newfit$coefficients)
    nvar <- length(varnames)
    ndead <- length(ssresid) / nvar
    times <- as.numeric(rownames(ssresid))
    if (is.character(transform)) {
      tname <- transform
      ttimes <- switch(
        transform,
        'identity' = times,
        'rank'    = rank(times),
        'log'     = log(times),
        'km' = {
          temp <- survival::survfitKM(factor(rep(1, nrow(
            newfit$y
          ))), newfit$y,
          se.fit = FALSE)
          # A nuisance to do left cont KM
          t1 <- temp$surv[temp$n.event > 0]
          t2 <- temp$n.event[temp$n.event > 0]
          km <- rep(c(1, t1), c(t2, 0))
          if (is.null(attr(sresid, 'strata')))
            1 - km
          else
            (1 - km[sort.list(sort.list(times))])
        },
        stop("Unrecognized transform")
      )
    } else {
      tname <- deparse(substitute(transform))
      if (length(tname) > 1)
        tname <- 'user'
      ttimes <- transform(times)
    }
    xx <- ttimes - mean(ttimes)
    cumInvVar <- prevInvVar + solve(newfit$var)
    cumVar <- solve(cumInvVar)
    cumBeta <- cumVar %*% (prevInvVar %*% prevBeta +
                             solve(newfit$var) %*% newfit$coefficients)
    cumFit <-
      do.call(
        survival::coxph,
        list(
          formula = formula,
          data = newdata,
          init = cumBeta,
          control = survival::coxph.control(iter.max = 0)
        )
      )
    V <- cumFit$var
    cA <- sum(xx ^ 2) * V / ndead
    nssresid <- stats::resid(cumFit, "schoenfeld")
    ## schoenfeld residual
    nsresid <- nssresid %*% V
    cU <- c(xx %*% nsresid)
    cumA <- cumA + cA
    cumU <- cumU + cU
    cumZ <- c(t(cumU) %*% solve(cumA, cumU))
    output <- list(
      cumA = cumA,
      cumU = cumU,
      beta = cumBeta,
      invVar = cumInvVar,
      Stat = cumZ,
      p = 1 - pchisq(cumZ, ncol(sresid))
    )
    class(output) <- "bczph.cee"
    return(output)
  }


#' @title  Read data
#' @description The function that reads the data used in the example
#' @param file_name  The file name of data
#'
#' @return The data of a example
#' @export

packagedata <- function(file_name) {
  file_path <-system.file("extdata", file_name, package = "StreamSurvDataMA")
  read.csv(file_path)
}


