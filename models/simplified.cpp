// This is a simple TMB model designed to match the original ADMB model from Stan's 2013 paper.

// Cole Monnahan; 9/2018

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA_INTEGER(ntows);
  // DATA_INTEGER(nlayers);
  DATA_VECTOR(BSA);
  DATA_VECTOR(BD);
  DATA_VECTOR(sum_SA1);
  DATA_VECTOR(sum_SA2);
  int ntows=BSA.rows();
  vector<Type> eta(ntows);
  vector<Type> eta_c(ntows);
  vector<Type> BSA_hat(ntows);
  
  //
  PARAMETER(log_q);
  PARAMETER(log_a);
  PARAMETER(d2);
  PARAMETER(b_BD);
  PARAMETER(log_c);
  PARAMETER(logSigma);
  PARAMETER(cb_BD);
  Type q=exp(log_q);
  Type a=exp(log_a);
  Type c=exp(log_c);

  Type nll=0.0;	   

  eta = d2+b_BD*BD;
  eta_c = cb_BD*BD;
  BSA_hat=1/(1/
   (q*sum_SA1 +exp(eta)*sum_SA2 + exp(eta_c)*c)+1/a);
  Type sigmasq=exp(2*logSigma);
  nll=0.5*(ntows*log(2*PI*sigmasq)+
	   pow(log(BSA)-log(BSA_hat),Type(2)).sum()/sigmasq);
  ADREPORT(a);
  ADREPORT(q);
  ADREPORT(c);
  vector<Type> d1=exp(eta)*sum_SA2 + exp(eta_c)*c;
  REPORT(d1);
  REPORT(BSA_hat);
  return(nll);
}

