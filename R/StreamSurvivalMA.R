#' @title             The survival prediction of renewable model averaging and its extension
#' @description       The function provides the survival prediction results of
#'      the cumulative estimating equation model averaging method (CEE),
#'      and the cumulative updated estimating equation model averaging method (CUEE),
#'      the renewable model averaging method and its extension for stream survival data.
#'      The function achieves the survival prediction of individuals
#'      and the online-updating weights under stream datasets and calculates the running time.
#'      The parameter estimators of every submodel can be obtained in it.
#' @param t           Time ponit
#' @param par_ori     The intitial parameter value
#' @param data_train  The training dataset
#' @param data_test   The test dataset
#' @param p           The dimension of covariates used to construct the submodels
#' @param e           The dimension of important covarites contained in each submodel
#' @param n_block     The number of individuals at every stream data batch
#'
#' @return            The predictive survival probability, online updating parameters of submodels, online updating model averaging weights and computation time
#' @export
#' @examples   library(survival); data_train <- packagedata("data_Train.csv")[,-1];
#'     data_test <- packagedata("data_Test.csv")[,-1];
#'     t=median(data_test$time); p=6; e=1; par_ori<-c(rep(1/(p-1),p-1));
#'     n_block=200; res<-StreamsurMA(t,par_ori,data_train,data_test,p,e,n_block)

StreamsurMA <- function(t,par_ori,data_train,data_test,p,e,n_block){

  ### renew method
  t.start <- proc.time()
  res_Renew<-MA_Renew(par_ori,data_train,data_test,p,e,n_block)  ### The renewable model averaging method
  t.use.Renew <- as.numeric(proc.time()[3]-t.start[3]) ### the running time

  ##### renew method with Taylor expansion of score function (The extension of renew method)
  t.start2 <- proc.time()
  res_ERenew<-MA_ERenew(par_ori,data_train,data_test,p,e,n_block)
  t.use.ERenew <- as.numeric(proc.time()[3]-t.start2[3])

  ### CUEE
  t.start3 <- proc.time()
  res_Cuee<-MA_Cuee(par_ori,data_train,data_test,p,e,n_block)
  t.use.Cuee <- as.numeric(proc.time()[3]-t.start3[3])

  ### Cee
  t.start4 <- proc.time()
  res_Cee<-MA_Cee(par_ori,data_train,data_test,p,e,n_block)
  t.use.Cee <- as.numeric(proc.time()[3]-t.start4[3])

  Presurp_Renew<-Stx(t,res_Renew$Predictvarible[[(max(data_train$block))]],res_Renew$Weight_renew[max(data_train$block),],data_test$time,data_test$status)
  Presurp_ERenew<-Stx(t,res_ERenew$Predictvarible_Erenew[[max(data_train$block)]],res_ERenew$Weight_Erenew[max(data_train$block),],data_test$time,data_test$status)
  Presurp_Cuee<-Stx(t,res_Cuee$Predictvarible_Cuee[[(max(data_train$block))]],res_Cuee$Weight_Cuee[max(data_train$block),],data_test$time,data_test$status)
  Presurp_Cee<-Stx(t,res_Cee$Predictvarible_Cee[[max(data_train$block)]],res_Cee$Weight_Cee[max(data_train$block),],data_test$time,data_test$status)

  ### output
  res_out <- list(PredictSurP_Renew=Presurp_Renew,
                  PredictSurP_ERenew=Presurp_ERenew,
                  PredictSurP_Cuee=Presurp_Cuee,
                  PredictSurP_Cee=Presurp_Cee,
                  MAWeight_Renew=res_Renew$Weight_renew[max(data_train$block),],
                  MAWeight_ERenew=res_ERenew$Weight_Erenew[max(data_train$block),],
                  MAWeight_Cuee=res_Cuee$Weight_Cuee[max(data_train$block),],
                  MAWeight_Cee=res_Cee$Weight_Cee[max(data_train$block),],
                  Para_Renew=res_Renew$Para_renew,
                  Para_ERenew=res_ERenew$Para_Erenew,
                  Para_Cuee=res_Cuee$Para_Cuee,
                  Para_Cee=res_Cee$Para_Cee,
                  Time_Renew=t.use.Renew,
                  Time_ERenew=t.use.ERenew,
                  Time_Cuee=t.use.Cuee,
                  Time_Cee=t.use.Cee

  )
  return(res_out)
}

# dir.create("inst/extdata", recursive = T)
# write.csv(data_train1[,1], file = "inst/extdata/trainset.csv")

# #dir.create("inst/extdata", recursive = T)
# write.csv(data_test, file = "inst/extdata/data_Test.csv")
