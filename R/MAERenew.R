
#' @title The extension of renewable model averaging method
#' @description  The function provides the results of the extended renewable model averaging method in stream datasets.
#' @param par_ori   The intitial parameter value
#' @param data_train  The training dataset
#' @param data_test   The test dataset
#' @param p  The dimension of covarites used to construct the submodels
#' @param e  The dimension of important covarites contained in every submodel
#' @param n_block The number of individuals at every stream data batch
#'
#' @return The extended renewable parameter estimators of submodels, the extended renewable model averaging weights and the predictive hazard ratio

MA_ERenew <- function(par_ori,data_train,data_test,p,e,n_block){

  n_b<-n_block
  max_iter <- 100
  tol=1e-8

  g_1<-c(2:(p+1))

  updatebeta <- list()
  updateJCB <- list()
  updateweight <- matrix(0, nrow = (max(data_train$block)), ncol = p)
  updatenew11<- list()
  pupdatenew11<- list()

  for (i_r in 1 : (max(data_train$block))){

    batchData <- data_train[data_train$block == i_r, ]

    MAJCB <- list()
    MAbeta <- list()
    new11_coxph <- matrix(0, nrow = n_block, ncol = p)
    prenew11_coxph<-matrix(0, nrow = dim(data_test)[1], ncol = p)

    for (m in 1:p){

      dims_renew <- m+e
      formula_m=paste("Surv(time,status)~",paste(colnames(batchData)[1:dims_renew],sep="",collapse = "+"),sep = "")
      formula_m=as.formula(formula_m)

      if(i_r==1){

        fit_m <- do.call(survival::coxph, list(formula = formula_m, data = batchData))
        beta <- as.matrix(fit_m$coefficients)
        JCB <- solve(fit_m$var)   #/n_b   #####negative Hessian matrix

      }else{

        beta_bef  <- updatebeta[[i_r-1]][[m]]
        JCB_bef <- updateJCB[[i_r-1]][[m]]

        Fit <-
          do.call(
            survival::coxph,
            list(
              formula = formula_m,
              data = batchData,
              init = beta_bef,
              control = coxph.control(iter.max = 100)
            )
          )

        Heinv <- solve(solve(Fit$var) + JCB_bef)
        score <- as.matrix(do.call(scorefunc,list(batchData=batchData, beta=beta_bef, dims=dims_renew)))   #/n_b ##### score function

        beta<-beta_bef+Heinv%*%score

        cumFit <-
          do.call(
            survival::coxph,
            list(
              formula = formula_m,
              data = batchData,
              init = beta,
              control = coxph.control(iter.max = 100)
            )
          )

        JCB <- JCB_bef + solve(cumFit$var)   #/n_b
        beta <- beta

      }

      MAJCB[[m]] <- JCB
      MAbeta[[m]] <- beta
      new11_coxph[,m] <- as.matrix(batchData[,1:dims_renew])%*%as.vector(beta)
      prenew11_coxph[,m] <- as.matrix(data_test[,1:dims_renew])%*%beta

    }

    colnames(new11_coxph) <- colnames(new11_coxph, do.NULL = FALSE, prefix = "R.")
    data_renew <- data.frame(new11_coxph,time=batchData$time,status=batchData$status)

    wformula_m_renew=as.formula(paste("Surv(time,status)~",paste(colnames(data_renew)[1:p],sep="",collapse = "+"),sep = ""))

    if(i_r==1){

      res_tay<-stats::optim(par=par_ori,fn=transcox_loglik,
                            control=list(maxit=1e4,fnscale=-1),method = "BFGS",
                            new11=new11_coxph,time=batchData$time,status=batchData$status)

      Wbeta <- res_tay$par

      w_renew <- append(x = exp(Wbeta)/(sum(exp(Wbeta))+1), 1/(sum(exp(Wbeta))+1), after = p-1)

      WFit <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew,
            data = data_renew,
            init = w_renew,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWJCB <- solve(WFit$var)    #/n_b
      oriscore <- as.matrix(do.call(scorefunc,list( batchData=data_renew, beta=w_renew, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron<-kronecker(t(oriscore),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron<-kronecker(t(g_1),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WJCB <- t(Demu(Wbeta)[[1]])%*%oriWJCB%*%Demu(Wbeta)[[1]]-Scokron%*%t(Demu(Wbeta)[[2]])+(1/2)*log(n_b)*(Pekron%*%t(Demu(Wbeta)[[2]]))

    }else{

      Wbeta_bef  <- Wbeta   ### mu in the last batch
      WJCB_bef <- WJCB

      orWbeta_bef  <- w_renew  ### w in the last batch

      WFit <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew,
            data = data_renew,
            init = orWbeta_bef,
            control = coxph.control(iter.max = 100)
          )
        )

      orWHeinv <- solve(WFit$var)    #/n_b

      orWscore <- as.matrix(do.call(scorefunc,list( batchData=data_renew, beta=orWbeta_bef, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron<-kronecker(t(orWscore),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron<-kronecker(t(g_1),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WHeinv <- solve(WJCB_bef+t(Demu(Wbeta_bef)[[1]])%*%orWHeinv%*%Demu(Wbeta_bef)[[1]]-Scokron%*%t(Demu(Wbeta_bef)[[2]])+(1/2)*log(n_b)*(Pekron%*%t(Demu(Wbeta_bef)[[2]])))
      Wscore <- t(Demu(Wbeta_bef)[[1]])%*%orWscore

      Wbeta<-Wbeta_bef+WHeinv%*%(Wscore-(1/2)*log(n_b)*t(Demu(Wbeta_bef)[[1]])%*%g_1)
      w_renew <- append(x = exp(Wbeta)/(sum(exp(Wbeta))+1), 1/(sum(exp(Wbeta))+1), after = p-1)

      WcumFit <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew,
            data = data_renew,
            init = w_renew,
            control = coxph.control(iter.max = 100)
          )
        )

      orWHeinv_now <- solve(WcumFit$var)    #/n_b
      orWscore_now <- as.matrix(do.call(scorefunc,list( batchData=data_renew, beta=w_renew, dims=p)))    #/n_b ##### score function

      Scokron_now<-kronecker(t(orWscore_now),I)  ### the part of 2 order deveriate of phi(mu)
      Pekron_now<-kronecker(t(g_1),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WJCB <-  WJCB_bef + t(Demu(Wbeta)[[1]])%*%orWHeinv_now%*%Demu(Wbeta)[[1]]-Scokron_now%*%t(Demu(Wbeta)[[2]])+(1/2)*log(n_b)*(Pekron_now%*%t(Demu(Wbeta)[[2]]))
    }

    updateweight[i_r,] <- t(w_renew)

    updatebeta[[i_r]] <- MAbeta
    updateJCB[[i_r]] <- MAJCB

    updatenew11[[i_r]] <- new11_coxph
    pupdatenew11[[i_r]] <- prenew11_coxph

  }

  ### output
  out <- list(Weight_Erenew=updateweight,
              Para_Erenew=updatebeta,
              Trainvarible_Erenew=updatenew11,
              Predictvarible_Erenew=pupdatenew11

  )
  return(out)
}

