
#' @title The renewable model averaging method
#' @description     The function provides the results of the renewable model averaging method in stream datasets.
#' @param par_ori    The intitial parameter value
#' @param data_train The training dataset
#' @param data_test  The test dataset
#' @param p  The dimension of covarites used to construct the submodels
#' @param e  The dimension of important covarites contained in every submodel
#' @param n_block The number of individuals at every stream data batch
#'
#' @return The renewable parameter estimators of submodels, the renewable model averaging weights and the predictive hazard ratio

MA_Renew <- function(par_ori,data_train,data_test,p,e,n_block){
  n_b<-n_block
  max_iter <- 100
  tol=1e-8

  g<-c(2:(p+1))

  updatebeta_renew <- list()
  updateJCB_renew <- list()
  updateweight_renew <- matrix(0, nrow = (max(data_train$block)), ncol = p)
  updatenew11_renew<- list()
  pupdatenew11_renew<- list()


  for (i_s in 1 : (max(data_train$block))){

    batchData_renew <- data_train[data_train$block == i_s, ]

    MAJCB_renew <- list()
    MAbeta_renew <- list()
    new11_renew <- matrix(0, nrow = n_block, ncol = p)
    prenew11_renew<-matrix(0, nrow = dim(data_test)[1], ncol = p)

    for (m_s in 1:p){

      dims_renew2 <- m_s + e
      formula_m_renew=paste("Surv(time,status)~",paste(colnames(batchData_renew)[1:dims_renew2],sep="",collapse = "+"),sep = "")
      formula_m_renew=as.formula(formula_m_renew)

      if(i_s==1){

        fit_m_renew <- do.call(survival::coxph, list(formula = formula_m_renew, data = batchData_renew))
        beta_renew <- as.matrix(fit_m_renew$coefficients)
        JCB_renew <- solve(fit_m_renew$var)   #/n_b   #####negative Hessian matrix

      }else{

        beta_bef_renew  <- updatebeta_renew[[i_s-1]][[m_s]]
        JCB_bef_renew <- updateJCB_renew[[i_s-1]][[m_s]]

        Fit_renew <-
          do.call(
            survival::coxph,
            list(
              formula = formula_m_renew,
              data = batchData_renew,
              init = beta_bef_renew,
              control = coxph.control(iter.max = 100)
            )
          )

        Heinv_renew <- solve(solve(Fit_renew$var) + JCB_bef_renew)     #/n_b

        rr <- 1

        beta_old_renew <- beta_bef_renew

        repeat{

          score_renew <- as.matrix(do.call(scorefunc,list( batchData=batchData_renew, beta=beta_old_renew, dims=dims_renew2)))  #/n_b ##### score function
          va_renew <- Heinv_renew%*%(JCB_bef_renew%*%(beta_bef_renew-beta_old_renew)+score_renew)
          beta_renew <- beta_old_renew + as.vector(va_renew)
          if( max(abs(beta_renew-beta_old_renew))>tol & rr<max_iter ){
            beta_old_renew <- beta_renew
            rr <- rr + 1
          }else{
            break
          }
        }

        cumFit_renew <-
          do.call(
            survival::coxph,
            list(
              formula = formula_m_renew,
              data = batchData_renew,
              init = beta_renew,
              control = coxph.control(iter.max = 100)
            )
          )

        JCB_renew <- JCB_bef_renew + solve(cumFit_renew$var)     #/n_b
        beta_renew <- beta_renew

      }

      MAJCB_renew[[m_s]] <- JCB_renew
      MAbeta_renew[[m_s]] <- beta_renew
      new11_renew[,m_s] <- as.matrix(batchData_renew[,1:dims_renew2])%*%as.vector(beta_renew)

      prenew11_renew[,m_s] <- as.matrix(data_test[,1:dims_renew2])%*%beta_renew

    }

    colnames(new11_renew) <- colnames(new11_renew, do.NULL = FALSE, prefix = "R.")
    data_renew2 <- data.frame(new11_renew,time=batchData_renew$time,status=batchData_renew$status)

    wformula_m_renew2=as.formula(paste("Surv(time,status)~",paste(colnames(data_renew2)[1:p],sep="",collapse = "+"),sep = ""))

    if(i_s==1){

      res<-stats::optim(par=par_ori,fn=transcox_loglik,
                        control=list(maxit=1e4,fnscale=-1), method = "BFGS",
                        new11=new11_renew,time=batchData_renew$time,status=batchData_renew$status)

      Wbeta_renew <- res$par
      w_renew2 <- append(x = exp(Wbeta_renew)/(sum(exp(Wbeta_renew))+1), 1/(sum(exp(Wbeta_renew))+1), after = p-1)

      WFit_renew <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew2,
            data = data_renew2,
            init = w_renew2,
            control = coxph.control(iter.max = 100)
          )
        )
      oriWJCB2 <- solve(WFit_renew$var)
      oriscore2 <- as.matrix(do.call(scorefunc,list( batchData=data_renew2, beta=w_renew2, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron_renew<-kronecker(t(oriscore2),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_renew<-kronecker(t(g),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WJCB_renew <- t(Demu(Wbeta_renew)[[1]])%*%oriWJCB2%*%Demu(Wbeta_renew)[[1]]-Scokron_renew%*%t(Demu(Wbeta_renew)[[2]])+(1/2)*log(n_b)*(Pekron_renew%*%t(Demu(Wbeta_renew)[[2]]))

    }else{

      orWbeta_bef_renew  <- w_renew2  ### w
      WJCB_bef_renew <- WJCB_renew   ### information matrix of mu
      Wbeta_bef_renew  <- Wbeta_renew ### mu

      WFit_renew <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew2,
            data = data_renew2,
            init = orWbeta_bef_renew,
            control = coxph.control(iter.max = 100)
          )
        )

      orWHeinv_renew <- solve(WFit_renew$var)    #/n_b
      orWscore_renew <- as.matrix(do.call(scorefunc,list( batchData=data_renew2, beta=orWbeta_bef_renew, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron_renew<-kronecker(t(orWscore_renew),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_renew<-kronecker(t(g),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WHeinv_renew <- solve(WJCB_bef_renew+t(Demu(Wbeta_bef_renew)[[1]])%*%orWHeinv_renew%*%Demu(Wbeta_bef_renew)[[1]]-Scokron_renew%*%t(Demu(Wbeta_bef_renew)[[2]])+(1/2)*log(n_b)*(Pekron_renew%*%t(Demu(Wbeta_bef_renew)[[2]])))

      rr_2 <- 1

      orWbeta_old_renew <- orWbeta_bef_renew
      Wbeta_old_renew <- Wbeta_bef_renew

      repeat{

        Wscore_renew <-  t(Demu(Wbeta_old_renew)[[1]])%*%as.matrix(do.call(scorefunc,list( batchData=data_renew2, beta=orWbeta_old_renew, dims=p)))   #/n_b ##### score function
        Wva_renew <- WHeinv_renew%*%(WJCB_bef_renew%*%(Wbeta_bef_renew-Wbeta_old_renew)+Wscore_renew-(1/2)*log(n_b)*t(Demu(Wbeta_old_renew)[[1]])%*%g)
        Wbeta_renew <- Wbeta_old_renew + as.vector(Wva_renew)
        if( max(abs(Wbeta_renew-Wbeta_old_renew))>tol & rr_2<max_iter ){
          Wbeta_old_renew <- Wbeta_renew
          orWbeta_old_renew <- append(x = exp(Wbeta_renew)/(sum(exp(Wbeta_renew))+1), 1/(sum(exp(Wbeta_renew))+1), after = p-1)
          rr_2 <- rr_2 + 1
        }else{
          break
        }
      }

      Wbeta_renew<-Wbeta_renew
      w_renew2<-append(x = exp(Wbeta_renew)/(sum(exp(Wbeta_renew))+1), 1/(sum(exp(Wbeta_renew))+1), after = p-1)

      WcumFit_renew <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_renew2,
            data = data_renew2,
            init = w_renew2,
            control = coxph.control(iter.max = 100)
          )
        )

      orWHeinv_renewnow <- solve(WcumFit_renew$var)    #/n_b

      orWscore_renewnow <- as.matrix(do.call(scorefunc,list( batchData=data_renew2, beta=w_renew2, dims=p)))    #/n_b ##### score function

      Scokron_renewnow<-kronecker(t(orWscore_renewnow),I)  ### the part of 2 order deveriate of phi(mu)
      Pekron_renewnow<-kronecker(t(g),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WJCB_renew <-  WJCB_bef_renew + t(Demu(Wbeta_renew)[[1]])%*%orWHeinv_renewnow%*%Demu(Wbeta_renew)[[1]]-Scokron_renewnow%*%t(Demu(Wbeta_renew)[[2]])+(1/2)*log(n_b)*(Pekron_renewnow%*%t(Demu(Wbeta_renew)[[2]]))

    }

    updateweight_renew[i_s,] <- t(w_renew2)

    updatebeta_renew[[i_s]] <- MAbeta_renew
    updateJCB_renew[[i_s]] <- MAJCB_renew

    updatenew11_renew[[i_s]] <- new11_renew
    pupdatenew11_renew[[i_s]] <- prenew11_renew

  }

  ### output
  out <- list(Weight_renew=updateweight_renew,
              Para_renew=updatebeta_renew,
              Trainvarible=updatenew11_renew,
              Predictvarible=pupdatenew11_renew

  )
  return(out)

}
