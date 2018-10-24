// This is the same as factor_analysis_pois.cpp but includes a spatial effect

#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX(Y_sp);       	// Responses
  DATA_INTEGER(n_f);        // Number of factors
  DATA_INTEGER(n_x); // Number of vertices in mesh
  DATA_IVECTOR(x_s);	// Convert site s to vertex x
  DATA_MATRIX(X_sj);		// Covariate design matrix
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER_MATRIX(beta_jp);
  PARAMETER_VECTOR(Loadings_vec);

  // Random effect
  PARAMETER_ARRAY(Omega_sf);
  PARAMETER(logsigma); // obs variance
  PARAMETER(logweight); // log weight for Poisson-link
  PARAMETER(log_kappa);


  Type sigma=exp(logsigma);
  Type sig22;
  sig22=pow(sigma,2)/2;
  
  using namespace density;
  int n_p = Y_sp.row(0).size();    // n_p is number of species
  int n_s = Y_sp.col(0).size();	   // num of observations
  int n_j = X_sj.row(0).size();	   // num of covariates
  Type jnll = 0;

  // Unpack loadings matrix
  matrix<Type> Loadings_pf(n_p, n_f);
  int Count = 0;
  for(int f=0; f<n_f; f++){
    for(int p=0; p<n_p; p++){
      if(p>=f){
        Loadings_pf(p,f) = Loadings_vec(Count);
        Count++;
      }else{
        Loadings_pf(p,f) = 0.0;
      }
    }
  }

  // Spatial variables
  Type log_tau = log( 1 / (exp(log_kappa) * sqrt(4*M_PI)) );  // Ensures that MargSD = 1
  Type Range = sqrt(8) / exp( log_kappa );
  Eigen::SparseMatrix<Type> Q = exp(log_kappa*4)*M0 + Type(2.0)*exp(log_kappa*2)*M1 + M2;
  for(int f=0; f<n_f; f++){
    jnll += SCALE( GMRF(Q), 1/exp(log_tau))( Omega_xf.col(f) );
  }
  
  // Calculate the numbers-density and weight in each strata
  matrix<Type> logweights(n_s, n_p);
  matrix<Type> lognumbers(n_s, n_p);
  matrix<Type> encounter( n_s, n_p);
  matrix<Type> catchrate( n_s, n_p);
  matrix<Type> logdensity(n_s, n_p); // predicted density in each strata
  vector<Type> BT_hat(n_s); // predicted BT catch
  matrix<Type> eta_sp( n_s, n_p );
  eta_sp = X_sj * beta_jp;
  for(int s=0; s<n_s; s++){
    for(int p=0; p<n_p; p++){
      // I assume the covariates and correlation between strata effects the
      // lognumbers not the weight
      lognumbers(s,p) = eta_sp(s,p);
      for(int f=0; f<n_f; f++){
	// This is the log-density in each of the three strata (p)
	lognumbers(s,p) += Omega_sf(s,f) * Loadings_pf(p,f);
      }
      logweights(s,p)=logweight; // assume constant weight
      encounter(s,p)=Type(1.0)-exp(-1*exp(lognumbers(s,p)));
      catchrate(s,p)=exp(lognumbers(s,p))*exp(logweights(s,p))/encounter(s,p);
      logdensity(s,p)=lognumbers(s,p)+logweights(s,p);
      // Now calculate the likelihood
      // Likelihood for the BT data (first column) I need to sum across the
      // first two strata. I guess this is where the correlation is induced.
      BT_hat(s)=log(exp(logdensity(s,0))+exp(logdensity(s,1)));
      if(Y_sp(s,p)==0){
	jnll-=log(1-encounter(s,p));
      } else {
	jnll-=log(encounter(s,p));
	// jnll-=dlnorm(catches(s,p), log(catchrate(s,p))-pow(sigma,2)/2, sigma, true);
	if(p==0)
	  jnll-=dlnorm(Y_sp(s,0), BT_hat(s) -sig22, sigma, true);
	if(p==1)
	  // The AT1 is just the middle strata
	  jnll-=dlnorm(Y_sp(s,1), logdensity(s,1)-sig22, sigma, true);
	if(p==2)
	  // and AT2 is the last strata
	  jnll-=dlnorm(Y_sp(s,2), logdensity(s,2)-sig22, sigma, true);
      }
    }
  }

  // Reporting
  REPORT(logdensity);
  REPORT( log_kappa );
  REPORT( log_tau );
  REPORT( Range );
  REPORT(BT_hat);
  REPORT( Omega_sf );
  REPORT( Loadings_pf );
  REPORT( jnll );
  ADREPORT(sigma);
  return jnll;
}
