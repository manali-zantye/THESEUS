*Battery design
nc_p.lo = 0;
nc_p.up = Pub;
nc_e.lo = 0;
nc_e.up = Eub;

*Battery outputs
y_cell.lo(t, 'SOC') = s_lb;
y_cell.up(t, 'SOC') = s_ub;
y_cell.lo(t, 'V') = V_lb;
y_cell.up(t, 'V') = V_ub;

flag_c.lo(t) = 0;
flag_c.up(t) = 1;
flag_d.lo(t) = 0;
flag_d.up(t) = 1;

dem_um.up(t) = ynom_ng('p_ng') + 5;
png_waste.up(t) = ynom_ng('p_ng') + 5;


yc.up(t) = Inom*V_ub; 
yd.up(t) = Inom*V_ub;

y_ng.lo(t,'p_ng') = minlf_ng*ynom_ng('p_ng');
y_ng.up(t,'p_ng') = ynom_ng('p_ng');


*y_cell.fx('1', 'V') = Vnom;
*y_cell.fx('2', 'V') = Vnom;




