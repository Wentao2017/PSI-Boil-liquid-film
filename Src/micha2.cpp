    /* turbulent surface tension */
    for_m(m){ 
      real rhol = liquid.rho()->value();
      real sigma = mixed.sigma()->value();
//      real kappa = conc.kappa()->value();
      for_vmijk(xyz,m,i,j,k){
    	  const real C_hmp2 = 1.0;
    	  const real a_hmp2 = 1.0;     /*f=kappa*(sigma+sigma*C_hmp2*(kappa*dl)^a_hmp2)*gradphi*/
 	  real dl = pow(xyz.dV(m,i,j,k),1.0/3.0);
     	  //real sigma_t = C_hmp * rhol * vmag2 * dl;
	  real kappa_limiter;
	  if(kappa[i][j][k] > 1e+23){
	      kappa_limiter = 0;
          }else{
              kappa_limiter = kappa[i][j][k];
          }  		  
  	  real sigma_t = sigma * C_hmp2 * pow((kappa_limiter * dl),a_hmp2);                        
          sigma_t = std::min(100.0*sigma,sigma_t);
          xyz[m][i][j][k] *= (sigma+sigma_t)/sigma;
          sig_t[m][i][j][k] = sigma_t;  // for visualization
      }
    }
