####For R package, LRT two-way Crossed ANOVA, Simultaneous ordered alternatives

################################################################################
##function to evaluate the test statistic value, critical value and the P-value of the test Psi1
#' To compute the test statistic value, critical value and the P-value of the test Psi1
#'
#' @param a number of levels of the row factor A
#' @param b number of levels of the column factor B
#' @param nominalsize the level of significance
#' @param cellsize matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param cellmean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param cellvar matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return Numeric vector consisting of observed test statistic, critical and P-values of the test Psi1 and the decision of the test
#' @export
#'
#' @examples
#' Psi_1(3,2,0.05,rbind(c(50,24),c(55,50),c(50,55)),rbind(c(69.8700,75.9833),c(73.2727,77.9140),c(77.4160,79.7090)),rbind(c(35.5678,10.6067),c(9.9028,12.9277),c(19.2650,15.0549)))
Psi_1<-function(a,b,nominalsize,cellsize,cellmean,cellvar){
  library(Iso)
  a<-a #number of levels of the row factor
  b<-b #number of levels of the column factor
  n<-cellsize #matrix of order a times b consisting of cell frequencies
  mu<-mean(cellmean)
  alpha<-apply(cellmean,1,mean)-mu
  beta<-apply(cellmean,2,mean)-mu
  nu<-cellmean-matrix(c(apply(cellmean,1,mean)),nrow=a,ncol=b,byrow = FALSE)-matrix(c(apply(cellmean,2,mean)),nrow=a,ncol=b,byrow=TRUE)+mu*matrix(1,a,b)
  var<-cellvar
  ####
  sam_mean<-cellmean #matrix of order a times b consisting of cell sample means
  sam_var<-cellvar   #matrix of order a times b consisting of cell sample variances
  mu_hat<-mu
  alpha_hat<-alpha
  beta_hat<-beta
  nu_hat<-nu

  ##MLEs under the null parameter space
  #initialization of the parameters under the null parameter space
  mu_null_0<-mu_hat
  var_null_0<-sam_var
  repeat
  {
    mu_null_1<-sum(n*sam_mean/var_null_0)/sum(n/var_null_0)
    var_null_1<-(((n-1)*sam_var)/n)+(sam_mean-mu_null_1)*(sam_mean-mu_null_1)
    if(max(abs(mu_null_1-mu_null_0),max(abs(var_null_1-var_null_0)))<0.00001)
    {
      break
    }
    mu_null_0<-mu_null_1
    var_null_0<-var_null_1
  }

  ####MLEs under the full parameter space
  ##initialization of the parameters under the full space
  var_full_0<-sam_var
  eta_full_0<-mu_hat*matrix(1,a,b)+matrix(c(beta_hat),a,b,byrow=TRUE)+matrix(c(alpha_hat),a,b,byrow=FALSE)+nu_hat
  repeat
  {
    w_matrix<-n/var_full_0
    eta_full_biviso<-biviso(y = sam_mean, w = w_matrix, eps = NULL, eps2 = 1e-9, ncycle = 50000)
    eta_full_1<-matrix(eta_full_biviso,a,b,byrow=FALSE)
    var_full_1<-(((n-1)*sam_var)/n)+(sam_mean-eta_full_1)*(sam_mean-eta_full_1)
    if(max(c(max(abs(eta_full_1-eta_full_0)),max(abs(var_full_1-var_full_0))))<0.00001)
    {
      break
    }
    var_full_0<-var_full_1
    eta_full_0<-eta_full_1
  }
  value_LRT<-(prod((var_full_0/var_null_0)^(n/2)))*exp(((1/2)*sum((1/var_full_1)*((n-1)*sam_var+n*(sam_mean-eta_full_1)*(sam_mean-eta_full_1))))-((1/2)*sum((1/var_null_1)*((n-1)*sam_var+n*(sam_mean-mu_null_1)*(sam_mean-mu_null_1)))))

  ##function to evaluate LRT value
  fun_LRT<-function(a,b,cellsize,cellmean,cellvar){
    a<-a
    b<-b
    n<-cellsize
    mu<-mean(cellmean)
    alpha<-apply(cellmean,1,mean)-mu
    beta<-apply(cellmean,2,mean)-mu
    nu<-cellmean-matrix(c(apply(cellmean,1,mean)),nrow=a,ncol=b,byrow = FALSE)-matrix(c(apply(cellmean,2,mean)),nrow=a,ncol=b,byrow=TRUE)+mu*matrix(1,a,b)
    var<-cellvar
    ##############################
    sam_mean<-array(NA,dim=c(a,b))
    sam_var<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data<-rnorm(n[i,j],mu+alpha[i]+beta[j]+nu[i,j],sqrt(var[i,j]))
        sam_mean[i,j]<-mean(data)
        sam_var[i,j]<-var(data)
      }
    }
    mu_hat<-sum(n*sam_mean)/sum(n)
    alpha_hat<-apply(sam_mean,1,mean)-mu_hat
    beta_hat<-apply(sam_mean,2,mean)-mu_hat
    nu_hat<-sam_mean-matrix(c(apply(sam_mean,1,mean)),nrow=a,ncol=b,byrow = FALSE)-matrix(c(apply(sam_mean,2,mean)),nrow=a,ncol=b,byrow=TRUE)+mu_hat*matrix(1,a,b)

    ##MLEs under the null parameter space
    #initialization of the parameters under the null space
    mu_null_0<-mu_hat
    var_null_0<-sam_var
    repeat
    {
      mu_null_1<-sum(n*sam_mean/var_null_0)/sum(n/var_null_0)
      var_null_1<-(((n-1)*sam_var)/n)+(sam_mean-mu_null_1)*(sam_mean-mu_null_1)
      if(max(abs(mu_null_1-mu_null_0),max(abs(var_null_1-var_null_0)))<0.00001)
      {
        break
      }
      mu_null_0<-mu_null_1
      var_null_0<-var_null_1
    }
    ##MLEs under the full parameter space
    #initialization of the parameters under the full space
    var_full_0<-sam_var
    eta_full_0<-mu_hat*matrix(1,a,b)+matrix(c(beta_hat),a,b,byrow=TRUE)+matrix(c(alpha_hat),a,b,byrow=FALSE)+nu_hat
    repeat
    {
      w_matrix<-n/var_full_0
      eta_full_biviso<-biviso(y = sam_mean, w = w_matrix, eps = NULL, eps2 = 1e-9, ncycle = 50000)
      eta_full_1<-matrix(eta_full_biviso,a,b,byrow=FALSE)
      var_full_1<-(((n-1)*sam_var)/n)+(sam_mean-eta_full_1)*(sam_mean-eta_full_1)
      if(max(c(max(abs(eta_full_1-eta_full_0)),max(abs(var_full_1-var_full_0))))<0.00001)
      {
        break
      }
      var_full_0<-var_full_1
      eta_full_0<-eta_full_1
    }
    LRT_value<-(prod((var_full_0/var_null_0)^(n/2)))*exp(((1/2)*sum((1/var_full_1)*((n-1)*sam_var+n*(sam_mean-eta_full_1)*(sam_mean-eta_full_1))))-((1/2)*sum((1/var_null_1)*((n-1)*sam_var+n*(sam_mean-mu_null_1)*(sam_mean-mu_null_1)))))
    return(LRT_value)
  }

  ####Parametric bootstrap step
  B<-5000 #number of bootstrap samples
  value_LRT_boot<-rep(NA,B)
  for(r in 1:B)
  {
    set.seed(17*r)
    value_LRT_boot[r]<-fun_LRT(a,b,n,array(0,c(a,b)),cellvar)
  }
  critical_value<-quantile(value_LRT_boot,nominalsize,names=F)
  P_value<-mean(value_LRT_boot<value_LRT)
  cat("The observed, critical and P-values of Psi_1 are respectively given by \n", value_LRT, critical_value, P_value)
  if(value_LRT<critical_value)
  {
    cat("\n Reject the null hypothesis at", nominalsize*100,"% level of significance")
  }
  else
  {
    cat("\n Do not reject the null hypothesis at", nominalsize*100,"% level of significance")
  }
}



################################################################################
##function to evaluate the test statistic value, critical value and the P-value of the test Psi2 and Psi3
#' To compute the test statistic values, critical values and the P-values of the test Psi2 and Psi3
#'
#' @param a number of levels of the row factor A
#' @param b number of levels of the column factor B
#' @param nominalsize the level of significance
#' @param cellsize matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param cellmean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param cellvar matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return Numeric vector consisting of observed test statistic, critical and P-values of the tests Psi2 and Psi3 and the decisions of these tests
#' @export
#'
#' @examples
#' Psi_2_Psi_3(3,2,0.05,rbind(c(50,24),c(55,50),c(50,55)),rbind(c(69.8700,75.9833),c(73.2727,77.9140),c(77.4160,79.7090)),rbind(c(35.5678,10.6067),c(9.9028,12.9277),c(19.2650,15.0549)))
Psi_2_Psi_3<-function(a,b,nominalsize,cellsize,cellmean,cellvar){
  library(Iso)
  a<-a
  b<-b
  n<-cellsize
  mu<-mean(cellmean)
  alpha<-apply(cellmean,1,mean)-mu
  beta<-apply(cellmean,2,mean)-mu
  nu<-cellmean-matrix(c(apply(cellmean,1,mean)),nrow=a,ncol=b,byrow = FALSE)-matrix(c(apply(cellmean,2,mean)),nrow=a,ncol=b,byrow=TRUE)+mu*matrix(1,a,b)
  var<-cellvar
  ####
  sam_mean<-cellmean
  sam_var<-cellvar

  ####
  T1<-array(NA,c((a-1),b))
  T2<-array(NA,c(a,(b-1)))
  for(i in 1:(a-1))
  {
    for(j in 1:b)
    {
      T1[i,j]<-(sam_mean[i+1,j]-sam_mean[i,j])/sqrt((sam_var[i+1,j]/n[i+1,j]+sam_var[i,j]/n[i,j]))
    }
  }
  for(i in 1:a)
  {
    for(j in 1:(b-1))
    {
      T2[i,j]<-(sam_mean[i,j+1]-sam_mean[i,j])/sqrt((sam_var[i,j+1]/n[i,j+1]+sam_var[i,j]/n[i,j]))
    }
  }
  value_1<-max(c(max(T1),max(T2)))
  value_2<-max(c(max(apply(T1,2,min)),max(apply(T2,1,min))))
  fun_simul<-function(a,b,cellsize,cellmean,cellvar){
    a<-a
    b<-b
    n<-cellsize
    mu<-mean(cellmean)
    alpha<-apply(cellmean,1,mean)-mu
    beta<-apply(cellmean,2,mean)-mu
    nu<-cellmean-matrix(c(apply(cellmean,1,mean)),nrow=a,ncol=b,byrow = FALSE)-matrix(c(apply(cellmean,2,mean)),nrow=a,ncol=b,byrow=TRUE)+mu*matrix(1,a,b)
    var<-cellvar
    ##############################
    sam_mean<-array(NA,dim=c(a,b))
    sam_var<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data<-rnorm(n[i,j],mu+alpha[i]+beta[j]+nu[i,j],sqrt(var[i,j]))
        sam_mean[i,j]<-mean(data)
        sam_var[i,j]<-var(data)
      }
    }

    ####
    T1<-array(NA,c((a-1),b))
    T2<-array(NA,c(a,(b-1)))
    for(i in 1:(a-1))
    {
      for(j in 1:b)
      {
        T1[i,j]<-(sam_mean[i+1,j]-sam_mean[i,j])/sqrt((sam_var[i+1,j]/n[i+1,j]+sam_var[i,j]/n[i,j]))
      }
    }
    for(i in 1:a)
    {
      for(j in 1:(b-1))
      {
        T2[i,j]<-(sam_mean[i,j+1]-sam_mean[i,j])/sqrt((sam_var[i,j+1]/n[i,j+1]+sam_var[i,j]/n[i,j]))
      }
    }
    value_1<-max(c(max(T1),max(T2)))
    value_2<-max(c(max(apply(T1,2,min)),max(apply(T2,1,min))))
    return(c(value_1,value_2))
  }

  ####Parametric bootstrap step
  B<-5000 #number of bootstrap samples
  value_simul_test1_boot<-rep(NA,B)
  value_simul_test2_boot<-rep(NA,B)
  for(r in 1:B)
  {
    set.seed(17*r)
    values_boot<-fun_simul(a,b,n,array(0,c(a,b)),cellvar)
    value_simul_test1_boot[r]<-values_boot[1]
    value_simul_test2_boot[r]<-values_boot[2]
  }
  critical_value_1<-quantile(value_simul_test1_boot,1-nominalsize,names=F)
  critical_value_2<-quantile(value_simul_test2_boot,1-nominalsize,names=F)
  P_value_1<-mean(value_simul_test1_boot>value_1)
  P_value_2<-mean(value_simul_test2_boot>value_2)

  cat("The observed, critical and P-values of Psi_2 are respectively given by \n", value_1, critical_value_1, P_value_1)
  if(value_1>critical_value_1)
  {
    cat("\n Reject the null hypothesis at", nominalsize*100,"% level of significance")
  }
  else
  {
    cat("\n Do not reject the null hypothesis at", nominalsize*100,"% level of significance")
  }

  cat("\n The observed, critical and P-values of Psi_3 are respectively given by \n", value_2, critical_value_2, P_value_2)
  if(value_2>critical_value_2)
  {
    cat("\n Reject the null hypothesis at", nominalsize*100,"% level of significance")
  }
  else
  {
    cat("\n Do not reject the null hypothesis at", nominalsize*100,"% level of significance")
  }
}
################################################################################



################################################################################
##function to evaluate one sided simultaneous CIs of the successive differences of the simple row and column effects
#' Title
#'
#' @param a number of levels of the row factor A
#' @param b number of levels of the column factor B
#' @param confidencelevel confidence level of the SCIs
#' @param cellsize matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param cellmean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param cellvar matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return Numerical values consisting of the lower limits of the one sided SCIs
#' @export
#'
#' @examples
#' SCI_Crossed(3,2,0.95,rbind(c(50,24),c(55,50),c(50,55)),rbind(c(69.8700,75.9833),c(73.2727,77.9140),c(77.4160,79.7090)),rbind(c(35.5678,10.6067),c(9.9028,12.9277),c(19.2650,15.0549)))
SCI_Crossed<-function(a,b,confidencelevel,cellsize,cellmean,cellvar){
  a<-a
  b<-b
  n<-cellsize
  mean_sam<-cellmean
  var_sam<-cellvar
  B<-5000 #number of bootstrap samples
  boot_value<-rep(NA,B)
  for(r in 1:B)
  {
    set.seed(17*r)
    sam_mean_boot<-array(NA,dim=c(a,b))
    sam_var_boot<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data_boot<-rnorm(n[i,j],mean_sam[i,j],sqrt(var_sam[i,j]))
        sam_mean_boot[i,j]<-mean(data_boot)
        sam_var_boot[i,j]<-var(data_boot)
      }
    }
    R_boot<-array(NA,dim=c(a-1,b))
    Q_boot<-array(NA,dim=c(a,b-1))
    for(i in 1:(a-1))
    {
      for(j in 1:b)
      {
        R_boot[i,j]<-((sam_mean_boot[i+1,j]-mean_sam[i+1,j])-(sam_mean_boot[i,j]-mean_sam[i,j]))/sqrt((sam_var_boot[i+1,j])/n[i+1,j]+(sam_var_boot[i,j])/n[i,j])
      }
    }
    for(i in 1:a)
    {
      for(j in 1:(b-1))
      {
        Q_boot[i,j]<-((sam_mean_boot[i,j+1]-mean_sam[i,j+1])-(sam_mean_boot[i,j]-mean_sam[i,j]))/sqrt((sam_var_boot[i,j+1])/n[i,j+1]+(sam_var_boot[i,j])/n[i,j])
      }
    }
    boot_value[r]<-max(c(max(R_boot),max(Q_boot)))
  }
  tau_approx<-quantile(boot_value,names = FALSE,probs = confidencelevel)
  ####
  CI_row_lower<-array(NA,dim=c(a-1,b))
  CI_col_lower<-array(NA,dim=c(a,b-1))
  for(i in 1:(a-1))
  {
    for(j in 1:b)
    {
      CI_row_lower[i,j]<-(mean_sam[i+1,j]-mean_sam[i,j])-tau_approx*sqrt((var_sam[i+1,j]/n[i+1,j])+(var_sam[i,j]/n[i,j]))
    }
  }
  for(i in 1:a)
  {
    for(j in 1:(b-1))
    {
      CI_col_lower[i,j]<-(mean_sam[i,j+1]-mean_sam[i,j])-tau_approx*sqrt((var_sam[i,j+1]/n[i,j+1])+(var_sam[i,j]/n[i,j]))
    }
  }
  cat("The lower limits of the one sided SCIs of Xi_2,1-Xi_1,1,...,Xi_a,1-Xi_a-1,1,...,Xi_2,b-Xi_1,b,...,Xi_a,b,Xi_a-1,b and
         Zeta_1,2-Zeta_1,1,...,Zeta_1,b-Zeta_1,b-1,...,Zeta_a,2-Zeta_a,1,...,Zeta_a,b-Zeta_a,b-1 \n at",
      confidencelevel*100,"% confidence level are respectively given by\n", c(CI_row_lower),"and \n",c(t(CI_col_lower)))
}
################################################################################


