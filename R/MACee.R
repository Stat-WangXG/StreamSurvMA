
#####cee method
#' @title  Online updating CEE model averaging method
#' @description  The function provides the results of online updating CEE model averaging method in stream datasets.
#' @param par_ori  The intitial parameter value
#' @param data_train  The training dataset
#' @param data_test  The test dataset
#' @param p  The dimension of covarites used to construct the submodels
#' @param e  The dimension of important covarites contained in every submodel
#' @param n_block  The number of individuals at every stream data batch
#'
#' @return The CEE parameter estimators of submodels, the CEE model averaging weights and the predictive hazard ratio

MA_Cee <- function(par_ori,data_train,data_test,p,e,n_block){

  n_b<-n_block
  g_cee<-c(2:(p+1))

  updatebeta_cee <- list()
  updatecuma_cee <- list()
  updatecumu_cee <- list()
  updateprevinvvar_cee <- list()

  updateweight_cee <- matrix(0, nrow = (max(data_train$block)), ncol = p)
  updatenew11_cee <- list()
  pupdatenew11_cee <- list()

  for (i_2 in 1:(max(data_train$block))){

    blockDat_cee <- data_train[data_train$block == i_2, ]

    MAbeta_cee <- list()
    MAcuma_cee <- list()
    MAcumu_cee <- list()
    MAprevinvvar_cee <- list()
    new11_cee <- matrix(0, nrow = n_block, ncol = p)
    prenew11_cee <- matrix(0, nrow = dim(data_test)[1], ncol = p)

    for ( m_2 in 1:p){

      dims_cee <- m_2+e

      formula_m_cee=paste("Surv(time,status)~",paste(colnames(blockDat_cee)[1:dims_cee],sep="",collapse = "+"),sep = "")
      formula_m_cee=as.formula(formula_m_cee)

      if(i_2==1){

        cuma_cee <- previnvvar_cee <- matrix(0, dims_cee, dims_cee)
        prevbeta_cee <- cumu_cee <- matrix(0, dims_cee, 1)

        bczphFit_cee <- bczph.cee(formula_m_cee,
                                  cumA = cuma_cee, cumU = cumu_cee, prevBeta = prevbeta_cee,
                                  prevInvVar = previnvvar_cee, newdata = blockDat_cee,
                                  transform = "km")

        cumA_cee <- bczphFit_cee$cumA
        cumU_cee <- bczphFit_cee$cumU
        beta_cee <- bczphFit_cee$beta
        invvar_cee <- bczphFit_cee$invVar

      }else{

        cuma_cee <- updatecuma_cee[[i_2-1]][[m_2]]
        cumu_cee <- updatecumu_cee[[i_2-1]][[m_2]]
        prevbeta_cee <- updatebeta_cee[[i_2-1]][[m_2]]
        previnvvar_cee <- updateprevinvvar_cee[[i_2-1]][[m_2]]

        bczphFit_cee <- bczph.cee(formula_m_cee,
                                  cumA = cuma_cee, cumU = cumu_cee, prevBeta = prevbeta_cee,
                                  prevInvVar = previnvvar_cee, newdata = blockDat_cee,
                                  transform = "km")

        cumA_cee <- bczphFit_cee$cumA
        cumU_cee <- bczphFit_cee$cumU
        beta_cee <- bczphFit_cee$beta
        invvar_cee <- bczphFit_cee$invVar

      }

      MAbeta_cee[[m_2]] <- beta_cee
      MAcuma_cee[[m_2]] <- cumA_cee
      MAcumu_cee[[m_2]] <- cumU_cee
      MAprevinvvar_cee[[m_2]] <- invvar_cee

      new11_cee[,m_2]<-as.matrix(blockDat_cee[,1:dims_cee])%*%as.vector(beta_cee)

      prc_m_cee=as.matrix(cbind(data_test[,1:dims_cee]))
      prenew11_cee[,m_2]<-prc_m_cee%*%beta_cee
    }

    colnames(new11_cee) <- colnames(new11_cee, do.NULL = FALSE, prefix = "T.")
    data_cee <- as.data.frame(cbind(new11_cee=new11_cee, time=blockDat_cee$time, status=blockDat_cee$status))

    wformula_m_cee=paste("Surv(time,status)~",paste(colnames(data_cee)[1:p],sep="",collapse = "+"),sep = "")
    wformula_m_cee=as.formula(wformula_m_cee)

    if(i_2==1){

      res_cee <- stats::optim(par=par_ori,fn=transcox_loglik,
                              control=list(maxit=1e4,fnscale=-1),method = 'BFGS',
                              new11=new11_cee,time=blockDat_cee$time,status=blockDat_cee$status)

      w_cee <- res_cee$par

      Wbeta_cee <- append(x = exp(w_cee)/(sum(exp(w_cee))+1), 1/(sum(exp(w_cee))+1), after = p-1)

      WFit_cee <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cee,
            data = data_cee,
            init = Wbeta_cee,
            control = coxph.control(iter.max = 100)
          )
        )

      oriWinvVar_cee <- solve(WFit_cee$var)    #/n_b
      oriscore_cee <- as.matrix(do.call(scorefunc,list( batchData=data_cee, beta=Wbeta_cee, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron_cee<-kronecker(t(oriscore_cee),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron_cee<-kronecker(t(g_cee),I)   ### the part of 2 order deveriate of phi(mu) at the penalty

      WinvVar_cee <- t(Demu(w_cee)[[1]])%*%oriWinvVar_cee%*%Demu(w_cee)[[1]]-Scokron_cee%*%t(Demu(w_cee)[[2]])+(1/2)*log(n_b)*(Pekron_cee%*%t(Demu(w_cee)[[2]]))

      Infest<-WinvVar_cee%*%w_cee

    }else{

      WprevInvVar_cee <- WinvVar_cee
      WprevBeta_cee <- w_cee
      WprevInfest <- Infest

      res_cee2<-stats::optim(par=par_ori,fn=transcox_loglik,
                             control=list(maxit=1e4,fnscale=-1),method = 'BFGS',
                             new11=new11_cee,time=blockDat_cee$time,status=blockDat_cee$status)

      w_cee1 <- res_cee2$par
      Wbeta_cee1 <- append(x = exp(w_cee1)/(sum(exp(w_cee1))+1), 1/(sum(exp(w_cee1))+1), after = p-1)

      WFit_cee1 <-
        do.call(
          survival::coxph,
          list(
            formula = wformula_m_cee,
            data = data_cee,
            init = Wbeta_cee1,
            control = coxph.control(iter.max = 100)
          )
        )

      WInfor_cee1 <- solve(WFit_cee1$var)    #/n_b
      oriscore_cee1 <- as.matrix(do.call(scorefunc,list( batchData=data_cee, beta=Wbeta_cee1, dims=p)))    #/n_b ##### score function

      I<-diag(p-1)  ###identify matrix
      Scokron1_cee<-kronecker(t(oriscore_cee1),I)  ### the part of 2 order deveriate of phi(mu)

      Pekron1_cee<-kronecker(t(g_cee),I)

      WInfor_cee <- t(Demu(w_cee1)[[1]])%*%WInfor_cee1%*%Demu(w_cee1)[[1]]-Scokron1_cee%*%t(Demu(w_cee1)[[2]])+(1/2)*log(n_b)*(Pekron1_cee%*%t(Demu(w_cee1)[[2]]))

      WcumInvVar_cee <- WprevInvVar_cee + WInfor_cee
      WcumVar_cee <- solve(WcumInvVar_cee)   #####solve
      w_cee <- WcumVar_cee %*% (WprevInfest + WInfor_cee%*%w_cee1)

      Wbeta_cee <- append(exp(w_cee)/(sum(exp(w_cee))+1), 1/(sum(exp(w_cee))+1), after = p+1)

      WinvVar_cee = WcumInvVar_cee
      Infest <- WprevInfest + WInfor_cee%*%w_cee1
    }

    updateweight_cee[i_2,] <- t(Wbeta_cee)

    updatebeta_cee[[i_2]] <- MAbeta_cee
    updatecuma_cee[[i_2]] <- MAcuma_cee
    updatecumu_cee[[i_2]] <- MAcumu_cee
    updateprevinvvar_cee[[i_2]] <- MAprevinvvar_cee

    updatenew11_cee[[i_2]] <- new11_cee
    pupdatenew11_cee[[i_2]] <- prenew11_cee

  }
  ### output
  out <- list(Weight_Cee=updateweight_cee,
              Para_Cee=updatebeta_cee,
              Trainvarible_Cee=updatenew11_cee,
              Predictvarible_Cee=pupdatenew11_cee

  )
  return(out)

}

