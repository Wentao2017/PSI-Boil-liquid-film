    /* turbulent surface tension */
    for_m(m){ 
      real rhol = liquid.rho()->value();
      real sigma = mixed.sigma()->value();
//      real kappa = conc.kappa()->value();
      for_vmijk(xyz,m,i,j,k){
	      
    	  const real c_0_hmp = 0.5;
    	  const real c_sigma_hmp = 10.0;     /*f=kappa*(sigma+sigma*((1+(kappa*dl/c_0)^beta)^(gamma/beta)-1)*c_0*c_sigma)*gradphi*/
          const real beta_hmp = 5.0;
	  const real gamma_hmp = 1.0;

 	  real dl = pow(xyz.dV(m,i,j,k),1.0/3.0);
     	  //real sigma_t = C_hmp * rhol * vmag2 * dl;
	  real kappa_limiter;
	  if(kappa[i][j][k] > 1e+23){
	      kappa_limiter = 0;
          }else{
              kappa_limiter = kappa[i][j][k];
          }  		  
  	  real sigma_t = sigma * (pow((1 + pow((fabs(kappa_limiter) * dl / c_0_hmp), beta_hmp)), (gamma_hmp / beta_hmp)) - 1) * c_0_hmp * c_sigma_hmp;                
          sigma_t = std::min(100.0*sigma,sigma_t);
          xyz[m][i][j][k] *= (sigma+sigma_t)/sigma;
          sig_t[m][i][j][k] = sigma_t;  // for visualization
      }
    }
