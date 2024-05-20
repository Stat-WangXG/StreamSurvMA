
#' @title Online updating CUEE model averaging method
#' @description  The function provides the results of online updating CUEE model averaging method in stream datasets.
#' @param par_ori  The intitial parameter value
#' @param data_train The training dataset
#' @param data_test  The test dataset
#' @param p  The dimension of covarites used to construct the submodels
#' @param e  The dimension of important covarites contained in every submodel
#' @param n_block The number of individuals at every stream data batch
#'
#' @return The CUEE parameter estimators of submodels, the CUEE model averaging weights and the predictive hazard ratio

MA_Cuee <- function(par_ori,data_train,data_test,p,e,n_block){

  n_b<-n_block
  g_cuee<-c(2:(p+1))

  updatebeta_cuee <- list()
  updatecuma_cuee <- list()
  updatecumu_cuee <- list()
  updateprevinvvar_cuee <- list()
  updatecumscore_cuee <- list()
  updatecumprod_cuee <- list()

  updateweight_cuee <- matrix(0, nrow = (max(data_train$block)), ncol = p)
  updatenew11_cuee <- list()
  pupdatenew11_cuee <- list()

  for (i_1 in 1:(max(data_train$block))){

    blockDat_cuee <- data_train[data_train$block == i_1, ]

    MAbeta_cuee <- list()
    MAcuma_cuee <- list()
    MAcumu_cuee <- list()
    MAprevinvvar_cuee <- list()
    MAcumscore_cuee <- list()
    MAcumprod_cuee <- list()

    new11_cuee <- matrix(0, nrow = n_block, ncol = p)
    prenew11_cuee<-matrix(0, nrow = dim(data_test)[1], ncol = p)

    for (m_1 in 1:p){

      dims_cuee=m_1+e

      formula_m_cuee=paste("Surv(time,status)~",paste(colnames(blockDat_cuee)[1:dims_cuee],sep="",collapse = "+"),sep = "")
      formula_m_cuee=stats::as.formula(formula_m_cuee)

      if(i_1==1){

        cuma_cuee <- previnvvar_cuee <- matrix(0, dims_cuee, dims_cuee)
        prevbeta_cuee <- cumu_cuee <- prevcumscore_cuee <- prevcumprod_cuee <- matrix(0, dims_cuee, 1)

        bczphFit_cuee <- bczph.cuee(formula_m_cuee,
                                    cumA = cuma_cuee, cumU = cumu_cuee, prevBeta = prevbeta_cuee,
                                    prevInvVar = previnvvar_cuee, prevCumScore = prevcumscore_cuee,
                                    prevCumProd = prevcumprod_cuee, newdata = blockDat_cuee,
                                    transform = "km")

        cumA_cuee <- bczphFit_cuee$cumA
        cumU_cuee <- bczphFit_cuee$cumU
        beta_cuee <- bczphFit_cuee$beta
        invvar_cuee <- bczphFit_cuee$invVar
        cumscore_cuee <- bczphFit_cuee$cumScore
        cumprod_cuee <- bczphFit_cuee$cumProd

      }else{

        cuma_cuee <- updatecuma_cuee[[i_1-1]][[m_1]]
        cumu_cuee <- updatecumu_cuee[[i_1-1]][[m_1]]
        prevbeta_cuee <- updatebeta_cuee[[i_1-1]][[m_1]]
        previnvvar_cuee <- updateprevinvvar_cuee[[i_1-1]][[m_1]]
        prevcumscore_cuee <- updatecumscore_cuee[[i_1-1]][[m_1]]
        prevcumprod_cuee <- updatecumprod_cuee[[i_1-1]][[m_1]]

        bczphFit_cuee <- bczph.cuee(formula_m_cuee,
                                    cumA = cuma_cuee, cumU = cumu_cuee, prevBeta = prevbeta_cuee,
                                    prevInvVar = previnvvar_cuee, prevCumScore = prevcumscore_cuee,
                                    prevCumProd = prevcumprod_cuee, newdata = blockDat_cuee,
                                    transform = "km")

        cumA_cuee <- bczphFit_cuee$cumA
        cumU_cuee <- bczphFit_cuee$cumU
        beta_cuee <- bczphFit_cuee$beta
        invvar_cuee <- bczphFit_cuee$invVar
        cumscore_cuee <- bczphFit_cuee$cumScore
        cumprod_cuee <- bczphFit_cuee$cumProd

      }

      MAbeta_cuee[[m_1]] <- beta_cuee
      MAcuma_cuee[[m_1]] <- cumA_cuee
      MAcumu_cuee[[m_1]] <- cumU_cuee
      MAprevinvvar_cuee[[m_1]] <- invvar_cuee
      MAcumscore_cuee[[m_1]] <- cumscore_cuee
      MAcumprod_cuee[[m_1]] <- cumprod_cuee

      new11_cuee[,m_1]<-as.matrix(blockDat_cuee[,1:dims_cuee])%*%as.vector(beta_cuee)

      prc_m_cuee=as.matrix(cbind(data_test[,1:dims_cuee]))
      prenew11_cuee[,m_1]<-prc_m_cuee%*%beta_cuee

    }

    colnames(new11_cuee) <- LETTERS[1:p]
    data_cuee <- as.data.frame(cbind(new11_cuee=new11_cuee,time=blockDat_cuee$time,status=blockDat_cuee$status))

    wformula_m_cuee=paste("Surv(time,status)~",paste(colnames(data_cuee)[1:p],sep="",collapse = "+"),sep = "")
    wformula_m_cuee=as.formula(wformula_m_cuee)

    if(i_1==1){

      WprevInvVar_cuee <- matrix(0, (p-1), (p-1))
      WprevCumScore_cuee <- prevCumProd <- matrix(0, (p-1), 1)

      res_cuee <- stats::optim(par=par_ori,fn=transcox_loglik,
                               control=list(maxit=1e4,fnscale=-1), method = "BFGS",
                               new11=new11_cuee,time=blockDat_cuee$time,status=blockDat_cuee$status)

      w_cuee <- res_cuee$par

      WbetaHat1 <- append(x = exp(w_cuee)/(sum(exp(w_cuee))+1), 1/(sum(exp(w_cuee))+1), after = p-1)

      WFit_cuee <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cuee,
            data = data_cuee,
            init = WbetaHat1,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWinvVarHat_cuee <- solve(WFit_cuee$var)    #/n_b
      oriscore_cuee <- as.matrix(do.call(scorefunc,list( batchData=data_cuee, beta=WbetaHat1, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron_cuee<-kronecker(t(oriscore_cuee),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_cuee<-kronecker(t(g_cuee),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WinvVarHat_cuee <- t(Demu(w_cuee)[[1]])%*%oriWinvVarHat_cuee%*%Demu(w_cuee)[[1]]-Scokron_cuee%*%t(Demu(w_cuee)[[2]])+(1/2)*log(n_b)*(Pekron_cuee%*%t(Demu(w_cuee)[[2]]))

      WbetaCheck <- solve(WprevInvVar_cuee + WinvVarHat_cuee,
                          prevCumProd + WinvVarHat_cuee %*% w_cuee)

      oriWbetaCheck <- append(x = exp(WbetaCheck)/(sum(exp(WbetaCheck))+1), 1/(sum(exp(WbetaCheck))+1), after = p-1)

      WFit_cuee2 <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cuee,
            data = data_cuee,
            init = oriWbetaCheck,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWinvVarHat_cuee2 <- solve(WFit_cuee2$var)    #/n_b
      oriscore_cuee2 <- as.matrix(do.call(scorefunc,list( batchData=data_cuee, beta=oriWbetaCheck, dims=p)))    #/n_b ##### score function

      Scokron_cuee2<-kronecker(t(oriscore_cuee2),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_cuee2<-kronecker(t(g_cuee),I)

      WinvVarCheck <- t(Demu(WbetaCheck)[[1]])%*%oriWinvVarHat_cuee2%*%Demu(WbetaCheck)[[1]]-Scokron_cuee2%*%t(Demu(WbetaCheck)[[2]])+(1/2)*log(n_b)*(Pekron_cuee2%*%t(Demu(WbetaCheck)[[2]]))

      Wscore_cuee <- WprevCumScore_cuee + t(Demu(WbetaCheck)[[1]])%*%oriscore_cuee2

      ## betaTilde
      WbetaTilde <- solve(WprevInvVar_cuee + WinvVarCheck,
                          prevCumProd + WinvVarCheck %*% WbetaCheck + Wscore_cuee-(1/2)*i_1*log(n_b)*t(Demu(WbetaCheck)[[1]])%*%g_cuee)

      wbeta_cuee <- append(x = exp(WbetaTilde)/(sum(exp(WbetaTilde))+1), 1/(sum(exp(WbetaTilde))+1), after = p-1)

      WinvVar = WprevInvVar_cuee + WinvVarCheck
      cumScore = Wscore_cuee
      cumProd = prevCumProd + WinvVarCheck %*% WbetaCheck

    }else{

      WprevInvVar_cuee <- WinvVar
      WprevCumScore_cuee <- cumScore
      prevCumProd <- cumProd

      res_cuee <- stats::optim(par=par_ori,fn=transcox_loglik,
                               control=list(maxit=1e4,fnscale=-1), method = "BFGS",
                               new11=new11_cuee,time=blockDat_cuee$time,status=blockDat_cuee$status)

      w_cuee <- res_cuee$par

      WbetaHat1 <- append(x = exp(w_cuee)/(sum(exp(w_cuee))+1), 1/(sum(exp(w_cuee))+1), after = p-1)

      WFit_cuee <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cuee,
            data = data_cuee,
            init = WbetaHat1,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWinvVarHat_cuee <- solve(WFit_cuee$var)    #/n_b
      oriscore_cuee <- as.matrix(do.call(scorefunc,list( batchData=data_cuee, beta=WbetaHat1, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron_cuee<-kronecker(t(oriscore_cuee),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_cuee<-kronecker(t(g_cuee),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WinvVarHat_cuee <- t(Demu(w_cuee)[[1]])%*%oriWinvVarHat_cuee%*%Demu(w_cuee)[[1]]-Scokron_cuee%*%t(Demu(w_cuee)[[2]])+(1/2)*log(n_b)*(Pekron_cuee%*%t(Demu(w_cuee)[[2]]))

      WbetaCheck <- solve(WprevInvVar_cuee + WinvVarHat_cuee,
                          prevCumProd + WinvVarHat_cuee %*% w_cuee)

      oriWbetaCheck <- append(x = exp(WbetaCheck)/(sum(exp(WbetaCheck))+1), 1/(sum(exp(WbetaCheck))+1), after = p-1)

      WFit_cuee2 <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cuee,
            data = data_cuee,
            init = oriWbetaCheck,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWinvVarHat_cuee2 <- solve(WFit_cuee2$var)    #/n_b
      oriscore_cuee2 <- as.matrix(do.call(scorefunc,list( batchData=data_cuee, beta=oriWbetaCheck, dims=p)))    #/n_b ##### score function

      Scokron_cuee2<-kronecker(t(oriscore_cuee2),I)  ### the part of 2 order deveriate of phi(mu)
      Pekron_cuee2<-kronecker(t(g_cuee),I)

      WinvVarCheck <- t(Demu(WbetaCheck)[[1]])%*%oriWinvVarHat_cuee2%*%Demu(WbetaCheck)[[1]]-Scokron_cuee2%*%t(Demu(WbetaCheck)[[2]])+(1/2)*log(n_b)*(Pekron_cuee2%*%t(Demu(WbetaCheck)[[2]]))
      Wscore_cuee <- WprevCumScore_cuee + t(Demu(WbetaCheck)[[1]])%*%oriscore_cuee2

      WbetaTilde <- solve(WprevInvVar_cuee + WinvVarCheck,
                          prevCumProd + WinvVarCheck %*% WbetaCheck + Wscore_cuee-(1/2)*i_1*log(n_b)*t(Demu(WbetaCheck)[[1]])%*%g_cuee)

      wbeta_cuee <- append(x = exp(WbetaTilde)/(sum(exp(WbetaTilde))+1), 1/(sum(exp(WbetaTilde))+1), after = p-1)

      WinvVar = WprevInvVar_cuee + WinvVarCheck
      cumScore = Wscore_cuee
      cumProd = prevCumProd + WinvVarCheck %*% WbetaCheck

    }

    updateweight_cuee[i_1,] <- t(wbeta_cuee)

    updatebeta_cuee[[i_1]] <- MAbeta_cuee
    updatecuma_cuee[[i_1]] <- MAcuma_cuee
    updatecumu_cuee[[i_1]] <- MAcumu_cuee
    updateprevinvvar_cuee[[i_1]] <- MAprevinvvar_cuee
    updatecumscore_cuee[[i_1]] <- MAcumscore_cuee
    updatecumprod_cuee[[i_1]] <- MAcumprod_cuee

    updatenew11_cuee[[i_1]] <- new11_cuee
    pupdatenew11_cuee[[i_1]] <- prenew11_cuee

  }

  ### output
  out <- list(Weight_Cuee=updateweight_cuee,
              Para_Cuee=updatebeta_cuee,
              Trainvarible_Cuee=updatenew11_cuee,
              Predictvarible_Cuee=pupdatenew11_cuee

  )
  return(out)

}


