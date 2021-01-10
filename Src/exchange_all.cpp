/******************************************************************************/
void exchange_all(Scalar & sca, real v_ex) {

  int dir = -1;
  const int s_x=sca.si();
  const int e_x=sca.ei();
  const int s_y=sca.sj();
  const int e_y=sca.ej();
  const int s_z=sca.sk();
  const int e_z=sca.ek();
  const int o_x=sca.ox();
  const int o_y=sca.oy();
  const int o_z=sca.oz();

  if( sca.bc().count() == 0 ) {
    OMS(Warning: exchanging a variable without boundary conditions);
  }

  /*-------------------------------+
  |                                |
  |  buffers for parallel version  |
  |                                |
  +-------------------------------*/
  real * sbuff_s, * sbuff_e, * rbuff_s, * rbuff_e;

  par_request req_s1,req_s2,req_r1,req_r2;
  par_status  ps;

  /* allocate memory for buffers */
  const int n = boil::maxi(sca.ni(), sca.nj(), sca.nk())+1;
  assert(n > 0);

  sbuff_s = new real [n*n];
  sbuff_e = new real [n*n];
  rbuff_s = new real [n*n];
  rbuff_e = new real [n*n];

  /*----------------+
  |  I - direction  |
  +----------------*/
  if(dir == -1 || dir == 0) {

    /* not decomposed in I direction */
    if( sca.domain()->dim(Comp::i()) == 1 ) {
      if( sca.bc().type(Dir::imin(), BndType::periodic()) && 
          sca.bc().type(Dir::imax(), BndType::periodic()) ){
#if 0
       std::cout<<"Gotcha!! "<<dir<<"\n";
       std::cout<<"s_x, e_x= "<<s_x<<" "<<e_x<<" "<<o_x<<"\n";
       std::cout<<e_x+1<<" "<<s_x + o_x<<" "<<v_ex<<"\n";
#endif
        real v_min =  v_ex;
        real v_max = -v_ex;
        for_avjk(sca,j,k) {
          sca[e_x+1][j][k] = sca[s_x + o_x][j][k] + v_max;
          sca[s_x-1][j][k] = sca[e_x - o_x][j][k] + v_min;
        }
      }
      //std::cout<<"I not decomposed\n";
      //exit(0);
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_avjk(sca,j,k) {
        int l = k*sca.nj()+j;
        sbuff_e[l] = sca[e_x - o_x][j][k];   // buffer i end
        sbuff_s[l] = sca[s_x + o_x][j][k];   // buffer i start
        rbuff_e[l] = sca[e_x +  1 ][j][k];
        rbuff_s[l] = sca[s_x -  1 ][j][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], sca.nj()*sca.nk(),
                            par_real, sca.domain()->neighbour(Dir::imax()),
                                      sca.domain()->neighbour(Dir::imin()), Tag(0));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], sca.nj()*sca.nk(),
                            par_real, sca.domain()->neighbour(Dir::imin()),
                                      sca.domain()->neighbour(Dir::imax()), Tag(1));
#else
      for_avjk(sca,j,k) {
        int l = k*sca.nj()+j;
        rbuff_e[l] = sca[e_x +  1 ][j][k];
        rbuff_s[l] = sca[s_x -  1 ][j][k];
      }

      if( sca.domain()->neighbour(Dir::imin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], sca.nj()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::imin()), Tag(0), & req_r1 );
      }
      if( sca.domain()->neighbour(Dir::imax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], sca.nj()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::imax()), Tag(1), & req_r2 );
      }

      for_avjk(sca,j,k) {
        int l = k*sca.nj()+j;
        sbuff_e[l] = sca[e_x - o_x][j][k];   // buffer i end
        sbuff_s[l] = sca[s_x + o_x][j][k];   // buffer i start
      }

      if( sca.domain()->neighbour(Dir::imax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], sca.nj()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::imax()), Tag(0), & req_s1 );
      }
      if( sca.domain()->neighbour(Dir::imin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], sca.nj()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::imin()), Tag(1), & req_s2 );
      }

      if( sca.domain()->neighbour(Dir::imax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( sca.domain()->neighbour(Dir::imin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( sca.domain()->neighbour(Dir::imin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( sca.domain()->neighbour(Dir::imax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
#if 1
      real v_min =  v_ex;
      real v_max = -v_ex;
      if(sca.bc().type_decomp(Dir::imin())) v_min=0.0;
      if(sca.bc().type_decomp(Dir::imax())) v_max=0.0;
      for_avjk(sca,j,k) {
        int l = k*sca.nj()+j;
        sca[e_x+1][j][k] = rbuff_e[l] + v_max;
        sca[s_x-1][j][k] = rbuff_s[l] + v_min;
      }
#else
      for_avjk(sca,j,k) {
        int l = k*sca.nj()+j;
        sca[e_x+1][j][k] = rbuff_e[l];   // buffer i end
        sca[s_x-1][j][k] = rbuff_s[l];   // buffer i start
      }
#endif
    }
  }
  
  /*----------------+
  |  J - direction  |
  +----------------*/
  if(dir == -1 || dir == 1) {

    /* not decomposed in J direction */
    if( sca.domain()->dim(Comp::j()) == 1 ) {
      if( sca.bc().type(Dir::jmin(), BndType::periodic()) && 
          sca.bc().type(Dir::jmax(), BndType::periodic()) )
        for_avik(sca,i,k) {
          sca[i][e_y+1][k] = sca[i][s_y + o_y][k];
          sca[i][s_y-1][k] = sca[i][e_y - o_y][k];
        }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_avik(sca,i,k) {
        int l = k*sca.ni()+i;
        sbuff_e[l] = sca[i][e_y - o_y][k];   // buffer j end
        sbuff_s[l] = sca[i][s_y + o_y][k];   // buffer j start
        rbuff_e[l] = sca[i][e_y +  1 ][k];
        rbuff_s[l] = sca[i][s_y -  1 ][k];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], sca.ni()*sca.nk(),
                            par_real, sca.domain()->neighbour(Dir::jmax()),
                                      sca.domain()->neighbour(Dir::jmin()), Tag(2));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], sca.ni()*sca.nk(),
                            par_real, sca.domain()->neighbour(Dir::jmin()),
                                      sca.domain()->neighbour(Dir::jmax()), Tag(3));
#else
      for_avik(sca,i,k) {
        int l = k*sca.ni()+i;
        rbuff_e[l] = sca[i][e_y +  1 ][k];
        rbuff_s[l] = sca[i][s_y -  1 ][k];
      }

      if( sca.domain()->neighbour(Dir::jmin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], sca.ni()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::jmin()), Tag(2), & req_r1 );
      }
      if( sca.domain()->neighbour(Dir::jmax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], sca.ni()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::jmax()), Tag(3), & req_r2 );
      }

      for_avik(sca,i,k) {
        int l = k*sca.ni()+i;
        sbuff_e[l] = sca[i][e_y - o_y][k];   // buffer j end
        sbuff_s[l] = sca[i][s_y + o_y][k];   // buffer j start
      }
  
      if( sca.domain()->neighbour(Dir::jmax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], sca.ni()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::jmax()), Tag(2), & req_s1 );
      }
      if( sca.domain()->neighbour(Dir::jmin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], sca.ni()*sca.nk(), par_real, 
                           sca.domain()->neighbour(Dir::jmin()), Tag(3), & req_s2 );
      }

      if( sca.domain()->neighbour(Dir::jmax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( sca.domain()->neighbour(Dir::jmin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( sca.domain()->neighbour(Dir::jmin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( sca.domain()->neighbour(Dir::jmax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
      for_avik(sca,i,k) {
        int l = k*sca.ni()+i;
        sca[i][e_y+1][k] = rbuff_e[l];   // buffer j end
        sca[i][s_y-1][k] = rbuff_s[l];   // buffer j start
      }
    }
  }
  
  /*----------------+
  |  K - direction  |
  +----------------*/
  if(dir == -1 || dir == 2) {

    /* not decomposed */
    if( sca.domain()->dim(Comp::k()) == 1 ) {
      if( sca.bc().type(Dir::kmin(), BndType::periodic()) && 
          sca.bc().type(Dir::kmax(), BndType::periodic()) )
        for_avij(sca,i,j) {
          sca[i][j][e_z+1] = sca[i][j][s_z + o_z];
          sca[i][j][s_z-1] = sca[i][j][e_z - o_z];
        }
    } 
    /* decomposed */
    else {
#ifdef SENDRECV
      for_avij(sca,i,j) {
        int l = j*sca.ni()+i;
        sbuff_e[l] = sca[i][j][e_z - o_z];   // buffer k end
        sbuff_s[l] = sca[i][j][s_z + o_z];   // buffer k start
        rbuff_e[l] = sca[i][j][e_z +  1 ];
        rbuff_s[l] = sca[i][j][s_z -  1 ];
      }
        /* send last and receive first */
        boil::cart.sendrecv(&sbuff_e[0], &rbuff_s[0], sca.ni()*sca.nj(),
                            par_real, sca.domain()->neighbour(Dir::kmax()),
                                      sca.domain()->neighbour(Dir::kmin()), Tag(4));

        /* send first and receive last */
        boil::cart.sendrecv(&sbuff_s[0], &rbuff_e[0], sca.ni()*sca.nj(),
                            par_real, sca.domain()->neighbour(Dir::kmin()),
                                      sca.domain()->neighbour(Dir::kmax()), Tag(5));
#else
      for_avij(sca,i,j) {
        int l = j*sca.ni()+i;
        rbuff_e[l] = sca[i][j][e_z +  1 ];
        rbuff_s[l] = sca[i][j][s_z -  1 ];
      }

      if( sca.domain()->neighbour(Dir::kmin()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_s[0], sca.ni()*sca.nj(), par_real,
                           sca.domain()->neighbour(Dir::kmin()), Tag(4), & req_r1 );
      }
      if( sca.domain()->neighbour(Dir::kmax()) != par_proc_null ) {
        boil::cart.irecv( &rbuff_e[0], sca.ni()*sca.nj(), par_real,
                           sca.domain()->neighbour(Dir::kmax()), Tag(5), & req_r2 );
      }

      for_avij(sca,i,j) {
        int l = j*sca.ni()+i;
        sbuff_e[l] = sca[i][j][e_z - o_z];   // buffer k end
        sbuff_s[l] = sca[i][j][s_z + o_z];   // buffer k start
      }

      if( sca.domain()->neighbour(Dir::kmax()) != par_proc_null ) {
        boil::cart.isend( &sbuff_e[0], sca.ni()*sca.nj(), par_real, 
                           sca.domain()->neighbour(Dir::kmax()), Tag(4), & req_s1 );
      }  
      if( sca.domain()->neighbour(Dir::kmin()) != par_proc_null ) {
        boil::cart.isend( &sbuff_s[0], sca.ni()*sca.nj(), par_real, 
                           sca.domain()->neighbour(Dir::kmin()), Tag(5), & req_s2 );
      }

      if( sca.domain()->neighbour(Dir::kmax()) != par_proc_null )
        boil::cart.wait( & req_s1, & ps );
      if( sca.domain()->neighbour(Dir::kmin()) != par_proc_null )
        boil::cart.wait( & req_s2, & ps );
      if( sca.domain()->neighbour(Dir::kmin()) != par_proc_null )
        boil::cart.wait( & req_r1, & ps );
      if( sca.domain()->neighbour(Dir::kmax()) != par_proc_null )
        boil::cart.wait( & req_r2, & ps );
#endif
      for_avij(sca,i,j) {
        int l = j*sca.ni()+i;
        sca[i][j][e_z+1] = rbuff_e[l];   // buffer k end
        sca[i][j][s_z-1] = rbuff_s[l];   // buffer k start
      }
    }
  }
  
  delete [] sbuff_s;
  delete [] sbuff_e;
  delete [] rbuff_s;
  delete [] rbuff_e;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_exchange_all.cpp,v 1.16 2011/03/28 07:44:59 sato Exp $'/
+-----------------------------------------------------------------------------*/
