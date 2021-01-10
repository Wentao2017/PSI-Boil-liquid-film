    /*---------------------------------+
    |  output free surface point data  |
    +---------------------------------*/
    const int n_planes = 4;
    const real dx = LX/real(NX);
    const real dx_planes = LX/real(n_planes);
    real x_fs[n_planes],z_fs[n_planes],z_fs2[n_planes];
    for(int i_plane = 0; i_plane<n_planes; i_plane++) {
      x_fs[i_plane]=-boil::exa;
      z_fs[i_plane]=-boil::exa;
      z_fs2[i_plane]=-boil::exa;
      for_vi(c,i){
        if((dx_planes*real(i_plane)<=c.xc(i))&&
           (c.xc(i)<dx_planes*real(i_plane)+dx)) {
          x_fs[i_plane]=c.xc(i);
          if(c.zn(c.sk())==0.0&&approx(c.yn(c.sj()),-0.5*LY)) {
#if 1
            for_vk(c,k) {
              if(c[i][c.sj()][k]<0.99){
                z_fs[i_plane]=c.zn(k)+c[i][c.sj()][k]*c.dzc(k);
                break;
              }
            }
#endif
#if 1
            real v0=c[i][c.sj()][c.sk()];
            if(v0<0.99){
              z_fs2[i_plane]=v0*c.dzc(c.sk());
              break;
            }
            for_vk(c,k) {
              real v1=c[i][c.sj()][k];
              if(v1<0.5){
                real r0=(0.5-v0)/(v1-v0);
                z_fs2[i_plane]=c.zc(k)*r0+c.zc(k-1)*(1.0-r0);
                break;
              }
              v0=v1;
            }
#endif
          }
        }
      }
    }
    boil::oout<<"z_fs= "<<time.current_time()<<" ";
    for(int i_plane = 0; i_plane<n_planes; i_plane++) {
      boil::cart.max_real(&x_fs[i_plane]);
      boil::cart.max_real(&z_fs[i_plane]);
      boil::oout<<x_fs[i_plane]<<" "<<z_fs[i_plane]<<" ";
    }
    boil::oout<<"\n";

    boil::oout<<"z_fs2= "<<time.current_time()<<" ";
    for(int i_plane = 0; i_plane<n_planes; i_plane++) {
      boil::cart.max_real(&x_fs[i_plane]);
      boil::cart.max_real(&z_fs2[i_plane]);
      boil::oout<<x_fs[i_plane]<<" "<<z_fs2[i_plane]<<" ";
    }
    boil::oout<<"\n";



    /*-----------------------+
    |  output sliced planes  |
    +-----------------------*/
    if((time.current_time())/(tint2) >= real(iint2) ) {
      for(int i_plane = 0; i_plane<n_planes; i_plane++) {
        int i=i_plane*(NX/n_planes);
        std::string fname;
        char* ch;
        fname = "c" + std::to_string(i_plane);
        ch = &fname[0];
        c.save_range(c.sI()+i, Range<int>(c.sJ(),c.eJ()),
                             Range<int>(c.sK(),c.eK()),
                             ch,iint2);
        m=Comp::u();
        fname = "u" + std::to_string(i_plane);
        ch = &fname[0];
        uvw.save_range(uvw.sI(m)+i,
                       Range<int>(uvw.sJ(m),uvw.eJ(m)),
                       Range<int>(uvw.sK(m),uvw.eK(m)),
                       m,ch,iint2);
        m=Comp::v();
        fname = "v" + std::to_string(i_plane);
        ch = &fname[0];
        uvw.save_range(uvw.sI(m)+i,
                       Range<int>(uvw.sJ(m),uvw.eJ(m)),
                       Range<int>(uvw.sK(m),uvw.eK(m)),
                       m,ch,iint2);
        m=Comp::w();
        fname = "w" + std::to_string(i_plane);
        ch = &fname[0];
        uvw.save_range(uvw.sI(m)+i,
                       Range<int>(uvw.sJ(m),uvw.eJ(m)),
                       Range<int>(uvw.sK(m),uvw.eK(m)),
                       m,ch,iint2);
      }
      iint2++;
    }
