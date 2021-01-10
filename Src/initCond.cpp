  /*-------------------+
  |  check if restart  |
  +-------------------*/
  std::fstream input;
  int irun = 0;
  if(boil::cart.iam()==0){
    input.open("run.txt", std::ios::in);
    if( !input.fail() ) {
      input >> irun;
      std::cout<<"read irun.  irun= "<<irun<<"\n";
    }
    input.close();
  }
  boil::cart.sum_int(&irun);
  if (irun==1){
    boil::oout<<"exit job due to irun=1"<<"\n";
    exit(0);
  }

  if(boil::cart.iam()==0){
    std::fstream output;
    output.open("run.txt", std::ios::out);
    output << 1 << boil::endl;
    output.close();
  }

  int ts=0;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    /* set time */
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
    /* load data */
    uvw.load("uvw", ts);
    press.load("press", ts);
    c.load("c", ts);
    uvw.exchange_all();
    press.exchange_all();
    c.exchange_all();
    /* set iint and iint2 */
    iint = int(time.current_time()/tint) + 1;
    iint2 = int(time.current_time()/tint2) + 1;
    boil::oout<<"iint= "<<iint<<" "<<iint2<<"\n";

  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "######################" << boil::endl;
    /*--------------------+
    |  initial condition  |
    +--------------------*/
    /* initial volume fraction */
    for_vijk(c,i,j,k)
      c[i][j][k] = 0.0;

    for_vijk(c,i,j,k){
      if(c.zn(k+1)<=filmThick_init){
        c[i][j][k]=1.0;
      } else if((c.zn(k)<filmThick_init)&&(filmThick_init<c.zn(k+1))){
        c[i][j][k]=(filmThick_init-c.zn(k))/c.dzc(k);
      } else if(LZ-filmThick_init<=c.zn(k)){
        c[i][j][k]=1.0;
      } else if((c.zn(k)<LZ-filmThick_init)&&(LZ-filmThick_init<c.zn(k+1))) {
        c[i][j][k]=(c.zn(k+1)-(LZ-filmThick_init))/c.dzc(k);
      }
    }
    c.exchange_all();
    conc.init();

    /* initial velocity */
    Comp m = Comp::u();
    real vf_l_ave = 2.0*filmThick_init/LZ;
    real vf_g_ave = 1.0 - vf_l_ave;
    boil::oout<<"# vf_l_ave, vf_g_ave= "<<vf_l_ave<<" "<<vf_g_ave<<"  "
              <<superfv_l_init/vf_l_ave<<" "<<superfv_g_init/vf_g_ave<<"\n";
    for_vmijk(uvw,m,i,j,k){
      if(c[i][j][k]>0.5){
        uvw[m][i][j][k] = superfv_l_init/vf_l_ave;
      } else {
        uvw[m][i][j][k] = superfv_g_init/vf_g_ave;
      }
    }
    uvw.exchange_all();

    /* initial time increment based on velocity */
    time.set_dt(0.25*LX/real(NX)/(superfv_g_init/vf_g_ave));
    boil::oout<<"# set dt= "<<0.25*LX/real(NX)/(superfv_g_init/vf_g_ave)<<"\n";

    /* output entire field */
    boil::plot->plot(uvw,c,press,mu_t,"uvw-c-press-mu_t",0);
    iint++;

    /* output free surface point data and sliced planes */
    #include "optfs.cpp"
  }
