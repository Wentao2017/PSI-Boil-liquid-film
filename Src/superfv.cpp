    /* calculate superficial velocity at cell center */
    real area_l=0.0;
    real area_v=0.0;
    real area_sum=0.0;
    real vel_l=0.0;
    real vel_v=0.0;

    for_vijk(c,i,j,k){
        real a_l = c.dSx(i,j,k) * c[i][j][k];
        real a_v = c.dSx(i,j,k) * (1.0- c[i][j][k]);
        real utmp = 0.5*(uvw[Comp::u()][i  ][j][k]
                       + uvw[Comp::u()][i+1][j][k]);
        vel_l  += a_l * utmp;
        vel_v  += a_v * utmp;
        area_l += a_l;
        area_v += a_v;
        area_sum += c.dSx(i,j,k);
    }

    boil::cart.sum_real(&area_l);
    boil::cart.sum_real(&area_v);
    boil::cart.sum_real(&area_sum);
    boil::cart.sum_real(&vel_l);
    boil::cart.sum_real(&vel_v);

    vel_l /= area_sum;
    vel_v /= area_sum;

    /* calculate superficial velocity at face center */

    real area_lf=0.0;
    real area_vf=0.0;
    real area_sumf=0.0;
    real vel_lf=0.0;
    real vel_vf=0.0;

    for_vijk(c,i,j,k){
        real a_lf = c.dSx(i,j,k) * 0.5 * (c[i-1][j][k]+c[i][j][k]);
        real a_vf = c.dSx(i,j,k) * 0.5 * (2.0-c[i-1][j][k]-c[i][j][k]);
        real utmp = uvw[Comp::u()][i][j][k];
        vel_lf  += a_lf * utmp;
        vel_vf  += a_vf * utmp;
        area_lf += a_lf;
        area_vf += a_vf;
        area_sumf += c.dSx(i,j,k);
    }

    boil::cart.sum_real(&area_lf);
    boil::cart.sum_real(&area_vf);
    boil::cart.sum_real(&area_sumf);
    boil::cart.sum_real(&vel_lf);
    boil::cart.sum_real(&vel_vf);

    vel_lf /= area_sumf;
    vel_vf /= area_sumf;

    boil::oout<<"superficialVel= "<<time.current_time()<<" "
              <<vel_l<<" "<<vel_v<<"  "
              <<vel_lf<<" "<<vel_vf<<"\n";

