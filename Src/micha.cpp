    /* turbulent surface tension */
    for_m(m){ 
      real rhol = liquid.rho()->value();
      real sigma = mixed.sigma()->value();
      const real C_hmp = 2.0;
      for_vmijk(xyz,m,i,j,k){
        real utmp = uvw[Comp::u()][i][j][k];
        real vtmp = uvw[Comp::v()][i][j][k];
        real wtmp = uvw[Comp::w()][i][j][k];
        if (m==Comp::u()) {
          vtmp = (uvw[Comp::v()][i-1][j][k] + uvw[Comp::v()][i-1][j+1][k]
                 +uvw[Comp::v()][i  ][j][k] + uvw[Comp::v()][i  ][j+1][k])*0.25;
          wtmp = (uvw[Comp::w()][i-1][j][k] + uvw[Comp::w()][i-1][j][k+1]
                 +uvw[Comp::w()][i  ][j][k] + uvw[Comp::w()][i  ][j][k+1])*0.25;
        }
        if (m==Comp::v()) {
          utmp = (uvw[Comp::u()][i][j-1][k] + uvw[Comp::u()][i+1][j-1][k]
                 +uvw[Comp::u()][i][j  ][k] + uvw[Comp::u()][i+1][j  ][k])*0.25;
          wtmp = (uvw[Comp::w()][i][j-1][k] + uvw[Comp::w()][i][j-1][k+1]
                 +uvw[Comp::w()][i][j  ][k] + uvw[Comp::w()][i][j  ][k+1])*0.25;
        }
        if (m==Comp::w()) {
          utmp = (uvw[Comp::u()][i][j][k-1] + uvw[Comp::u()][i+1][j][k-1]
                 +uvw[Comp::u()][i  ][j][k] + uvw[Comp::u()][i+1][j][k  ])*0.25;
          vtmp = (uvw[Comp::v()][i][j][k-1] + uvw[Comp::v()][i][j+1][k-1]
                 +uvw[Comp::v()][i][j  ][k] + uvw[Comp::v()][i][j+1][k  ])*0.25;
        } 
        real vmag2 = utmp*utmp + vtmp*vtmp + wtmp*wtmp;
        real dl = pow(xyz.dV(m,i,j,k),1.0/3.0);
        //real sigma_t = C_hmp * rhol * vmag2 * dl;
        real sigma_t = C_hmp * mixed.rho(m,i,j,k) * vmag2 * dl;
        sigma_t = std::min(100.0*sigma,sigma_t);
        xyz[m][i][j][k] *= (sigma+sigma_t)/sigma;
        sig_t[m][i][j][k] = sigma_t;  // for visualization
      }
    }
