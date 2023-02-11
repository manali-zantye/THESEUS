*Battery design
mod_bat.lo = 0;

mod_bat.up = max_mod;

nc_p.up = nc_po*max_mod*0.001;

nc_e.up = nc_eo*max_mod*0.001;

P_bat.lo(t) = -nc_po*max_mod*0.001;

P_bat.up(t) = nc_po*max_mod*0.001;




*Battery outputs
y_cell.lo(t, 'SOD') = 0.4;

y_cell.up(t, 'SOD') = 0.85;

y_cell.lo(t, 'V') = 1.61289;

y_cell.up(t, 'V') = 2.23266;


fix_curr.lo = 5.79358;

fix_curr.up = 6.92459;


*Added bounds

ecell.lo(t) = 0;

ecell.up(t) = 110;



flag_1.lo(t) = 0;

flag_1.up(t) = 1;

flag_2.lo(t) = 0;

flag_2.up(t) = 1;

flag_c.lo(t) = 0;

flag_c.up(t) = 1;

flag_d.lo(t) = 0;

flag_d.up(t) = 1;


df_cell.lo(t, 'V', st_eq) = -3;

df_cell.up(t, 'V', st_eq) = 3;

df1_cell.lo(t, 'V', st_eq) = -3;

df1_cell.up(t, 'V', st_eq) = 3;

df2_cell.lo(t, 'V', st_eq) = -3;

df2_cell.up(t, 'V', st_eq) = 3;

df3_cell.lo(t, 'V', st_eq) = 0;

df3_cell.up(t, 'V', st_eq) = 3;


*rc_ng.up(t) = 645;

dem_um.up(t) = 645;

png_waste.up(t) = 645;

P_grid.up(t) = 645;




yc.up(t) = nc_po*1000/cell_mod;

yd.up(t) = nc_po*1000/cell_mod;


y_ng.lo(t,'p_ng') = minlf_ng*ynom_ng('p_ng');

y_ng.up(t,'p_ng') = 645;

*fix_yng.lo = minlf_ng*ynom_ng('p_ng');

*fix_yng.up = 645;



