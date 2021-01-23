#include "Include/psi-boil.h"
#include "exchange_all.cpp"

/* domain size */
const real LX = 0.01;  // domain length in x (axial)
const real LY = 0.00125;  // domain length in y
const real LZ = 0.01;     // domain length in z (wall-wall)

/* grid size */
const int gLevel = 4;
const int NX  = 32*gLevel;  // uniform grid
const int NY  =  4*gLevel;
const int NZ  = 32*gLevel;

/* initial and boundary condition */
const real dpdx = 2.0*1.49e+4;     // pressure drop given as boundary condition
const real filmThick_init=2.5e-4;  // used for initial condition
const real superfv_l_init=0.5;     // used for initial condition
const real superfv_g_init=25.0;    // used for initial condition

/* constants */
const real gravity = 9.8;

/******************************************************************************/
main(int argc, char ** argv) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx ( Range<real>(0.0,LX), NX, Periodic::yes() );
  Grid1D gy ( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::yes() );
  Grid1D gz ( Range<real>(0.0,LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);           // vel
  Vector sig_t(d);           // output for turbulent_surface_tension
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d); // concentration
  Scalar mu_t(d), wd(d), ws(d);    // turbulence model

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  kappa = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 1.0 ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 1.0 ) );

  mu_t.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 0.0 ) );
  mu_t.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );
  wd = mu_t.shape();

  press.bc().add( BndCnd( Dir::imin(), BndType::periodic(),  dpdx*LX) );
  press.bc().add( BndCnd( Dir::imax(), BndType::periodic(), -dpdx*LX) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  exchange_all(press,dpdx*LX);

  /* copy b.c. from p */
  f = p.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter liquid(d),gas(d);
  liquid.rho   (997.6);
  liquid.mu    (8.887e-4);
  gas.rho      (1.184);
  gas.mu       (1.0447e-5);
  Matter mixed (liquid, gas, & c);
  mixed.sigma  (0.072);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 5000000;
  real tint = 0.005;       // time-interval for entire field data [s]
  real tint2 = 1.0/10000;  // time-interval for sliced planes [s], 10k Hz
  const int nint = 100000; // interval of step for bck files
  const real dxmin = d.dxyz_min();  // minimum grid spacing
  const real dtmax  = 1.0 * pow(gas.rho()->value()*pow(dxmin,3.0)
                    / (2.0*boil::pi*mixed.sigma()->value()),0.5);
  const real dt_init  = 5.0e-7;  // will be modified by initial velocity
  Times time(ndt, dt_init);
  boil::oout<<"dtmax= "<<dtmax<<"\n";
  /* set iint and iint2 */
  int iint = 0;   // counter for tint
  int iint2 = 0;  // counter for tint2

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.convection_set(TimeScheme::forward_euler());
  Pressure pr(p, f, uvw, time, solver, &mixed);
  VOF conc  (c,  g, kappa, uvw, time, solver);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /*--------+
  |  model  |
  +--------*/
  Model tm;
  Distance di(wd, ws, uvw, time, solver);
  di.compute();

  /*-------------------------------+
  |  initial condition or restart  |
  +-------------------------------*/
  #include "initCond.cpp"

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /*------------------+
    |  reset body force |
    +------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* surface tension */
    conc.tension(&xyz, mixed, c);

    /* turbulent surface tension */
#if 1
    #include "micha2.cpp"
#endif

    /* gravity force */
    Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k){
      xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    }

    /*-------------------+
    |  turbulence model  |
    +-------------------*/
    tm.smagorinsky( & ns, & mu_t, 0.173 );
    tm.tau_wall( & ns, wd, & xyz);

    /* essential for moving front */
    ns.discretize( & mu_t );
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-14));
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(pmin>press[i][j][k]) pmin=press[i][j][k];
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      press[i][j][k] -= pmin;
    }
    press.bnd_update();
    exchange_all(press,dpdx*LX);
#if 0
    /* limit velocity */
    for_m(m) {
      for_vmijk(uvw,m,i,j,k){
        uvw[m][i][j][k]=std::max(-100.0,std::min(100.0,uvw[m][i][j][k]));
      }
    }
    uvw.exchange_all();
#endif

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc.advance();
    conc.totalvol();

    /*-------------+
    |  dt control  |
    +-------------*/
    time.control_dt( ns.cfl_max(),0.25,dtmax);

    /*-----------------------------------+
    |  statistics: superfitial velocity  |
    +-----------------------------------*/
    #include "superfv.cpp"

    /*--------------------------------------------+
    |  output free surface point & sliced planes  |
    +--------------------------------------------*/
  //  #include "optfs.cpp"

    /*---------------------------+
    |  output entire field data  |
    +---------------------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      uvw.exchange_all();
      boil::plot->plot(uvw,c,press,mu_t,"uvw-c-press-mu_t",iint);
      boil::plot->plot(sig_t,c,"sig_t-c",iint);
      boil::plot->plot(xyz,c,"xyz-c",iint);
      boil::plot->plot(kappa,c,"kappa-c",iint);
      iint++;
    }

    /*---------------------+
    |  output backup data  |
    +---------------------*/
    #include "optbck.cpp"
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}

/*-----------------------------------------------------------------------------+
 '$Id: main-dam.cpp,v 1.12 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
