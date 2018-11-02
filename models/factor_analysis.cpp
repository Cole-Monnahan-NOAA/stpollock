/// A factor analysis model modified from Thorson's spatial_factor_analysis
// TMB model. Here I'm using it to try and closely match Stan's old method. 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX(Y_sp);       	// Responses
  DATA_INTEGER(n_f);        // Number of factors
  // DATA_INTEGER(n_x); // Number of vertices in mesh
  //  DATA_IVECTOR(x_s);	// Convert site s to vertex x
  DATA_MATRIX(X_sj);		// Covariate design matrix

  // Parameters
  PARAMETER_MATRIX(beta_jp);
  PARAMETER_VECTOR(Loadings_vec);

  // Random effect
  PARAMETER_ARRAY(Omega_sf);
  PARAMETER(logsigma); // obs variance
  Type Sigma=exp(logsigma);
  
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
  // Type log_tau = log( 1 / (exp(log_kappa) * sqrt(4*M_PI)) );  // Ensures that MargSD = 1
  // Type Range = sqrt(8) / exp( log_kappa );
  // Eigen::SparseMatrix<Type> Q = exp(log_kappa*4)*M0 + Type(2.0)*exp(log_kappa*2)*M1 + M2;
  // for(int f=0; f<n_f; f++){
  //   jnll += SCALE( GMRF(Q), 1/exp(log_tau))( Omega_xf.col(f) );
  // }
  for(int ii=0; ii<n_s; ii++){
    for(int jj=0; jj<n_f; jj++){
      jnll-=dnorm(Omega_sf(ii,jj), Type(0),  Type(1), true);
    }
  }
  // Likelihood contribution from observations
  matrix<Type> logdensity_sp( n_s, n_p );
  matrix<Type> eta_sp( n_s, n_p );
  eta_sp = X_sj * beta_jp;
  for(int s=0; s<n_s; s++){
  for(int p=0; p<n_p; p++){
    // Predictor
    logdensity_sp(s,p) = eta_sp(s,p);
    for(int f=0; f<n_f; f++){
      // This is the log-density in each of the three strata (p)
      logdensity_sp(s,p) += Omega_sf(s,f) * Loadings_pf(p,f);
    }
    // jnll -= dpois( Y_sp(s,p), exp(logdensity_sp(s,p)), true );
  }}

  vector<Type> BT_hat(n_s);
  // The likelihood 
  for(int s=0; s<n_s; s++){
    // Likelihood for the BT data (first column) I need to sum across the
    // first two strata. I guess this is where the correlation is induced.
    BT_hat(s)=log(exp(logdensity_sp(s,0))+exp(logdensity_sp(s,1)));
    // BT_hat(s)=logdensity_sp(s,0)+logdensity_sp(s,1);
    jnll-=dnorm(Y_sp(s,0), BT_hat(s), Sigma, true);
    // The AT1 is just the middle strata
    jnll-=dnorm(Y_sp(s,1), logdensity_sp(s,1), Sigma, true);
    // and AT2 is the last strata
    jnll-=dnorm(Y_sp(s,2), logdensity_sp(s,2), Sigma, true);
  }
  
  // Reporting
  // REPORT( log_kappa );
  //  REPORT( log_tau );
  // REPORT( Range );
  REPORT(BT_hat);
  REPORT( Omega_sf );
  REPORT( Loadings_pf );
  REPORT( jnll );
  REPORT( logdensity_sp );
  ADREPORT(Sigma);
  ADREPORT(Omega_sf(0,0));
  return jnll;
}
