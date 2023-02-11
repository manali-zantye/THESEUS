*Battery design
nc_p.lo = 0;
nc_p.up = 200;
nc_e.lo = 0;
nc_e.up = 800;



*Battery outputs
y_cell.lo(t, 'SOC') = 0.4;
y_cell.up(t, 'SOC') = 0.85;
y_cell.lo(t, 'V') = 3;
y_cell.up(t, 'V') = 4.25;


flag_c.lo(t) = 0;
flag_c.up(t) = 1;
flag_d.lo(t) = 0;
flag_d.up(t) = 1;


dem_um.up(t) = ynom_ng('p_ng') + 5;
png_waste.up(t) = ynom_ng('p_ng') + 5;




yc.up(t) = 25.5;
yd.up(t) = 25.5;


y_ng.lo(t,'p_ng') = minlf_ng*ynom_ng('p_ng');
y_ng.up(t,'p_ng') = ynom_ng('p_ng') + 5;




