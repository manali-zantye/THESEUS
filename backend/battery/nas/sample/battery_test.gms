option minlp = baron;
option optcr = 1e-3;
option optca = 1e-3;

Sets
t "decision stages: stage 1 is initialization, not decision" /1*%gams.user1%/

st_eq "cell state for the 4 equations" /2d, 1d, 1c, 2c/

st "actual cell states we need"  /c, d/

var "cell output variables"
/
SOD, V
/

inp "cell input variables: current(A)" /I_cell/

;

*NGCC sets
Sets
out_ng "NGCC output vars"
/p_ng/;
;



*master k,z,j,a,b,h,g sets for cell ROMs
Sets
$include "input_master_sets.gms"
;

*to determine limit on z,j,a,b,h,g
Sets
p1_map(var,st_eq,z)
n1_map(var,st_eq,j)
p2_map(var,st_eq,a)
n2_map(var,st_eq,b)
q1_map(var,st_eq,h)
r1_map(var,st_eq,g)
;

*actual limits on z,j,a,b,h,g (p1,n1,p2,n2,r1,q1)
Parameter
p1(var,st_eq)/
$include "p1_input_togams.gms"
/
n1(var,st_eq)/
$include "n1_input_togams.gms"
/
p2(var,st_eq)/
$include "p2_input_togams.gms"
/
n2(var,st_eq)/
$include "n2_input_togams.gms"
/
q1(var,st_eq)/
$include "q1_input_togams.gms"
/
r1(var,st_eq)/
$include "r1_input_togams.gms"
/
;

p1_map(var,st_eq,z) = no;
n1_map(var,st_eq,j) = no;
p2_map(var,st_eq,a) = no;
n2_map(var,st_eq,b) = no;
q1_map(var,st_eq,h) = no;
r1_map(var,st_eq,g) = no;

p1_map(var,st_eq,z)$(ord(z) le p1(var,st_eq)) = yes;
n1_map(var,st_eq,j)$(ord(j) le n1(var,st_eq)) = yes;
p2_map(var,st_eq,a)$(ord(a) le p2(var,st_eq)) = yes;
n2_map(var,st_eq,b)$(ord(b) le n2(var,st_eq)) = yes;
q1_map(var,st_eq,h)$(ord(h) le q1(var,st_eq)) = yes;
r1_map(var,st_eq,g)$(ord(g) le r1(var,st_eq)) = yes;


Set
k_map1(var,st_eq,k,z,j)
k_map2(var,st_eq,k,a,b)
k_map3(var,st_eq,k,h,g)
;


k_map1(var,st_eq,k,z,j) = no;
k_map1(var,st_eq,k,z,j)$(p1_map(var,st_eq,z) and n1_map(var,st_eq,j) and (ord(k) eq ( (ord(z) - 1)*n1(var,st_eq) + ord(j)))) = yes;


k_map2(var,st_eq,k,a,b) = no;
k_map2(var,st_eq,k,a,b)$(p2_map(var,st_eq,a) and n2_map(var,st_eq,b) and (ord(k) eq (n1(var,st_eq)*p1(var,st_eq) + (ord(a) - 1)*n2(var,st_eq) + ord(b)))) = yes;


k_map3(var,st_eq,k,h,g) = no;
k_map3(var,st_eq,k,h,g)$(q1_map(var,st_eq,h) and r1_map(var,st_eq,g) and (ord(k) eq (n1(var,st_eq)*p1(var,st_eq) + n2(var,st_eq)*p2(var,st_eq) + (ord(h)-1)*r1(var,st_eq) + ord(g)))) = yes;

Parameter theta(var,st_eq,k) "parameter vector for cell ROM"/
$include "input_h_vector.gms"
/;



*NGCC pars
Parameter

ynom_ng(out_ng)
/
$include pnom_ng.gms
/

P_dem(t) "power demand at time t for grid in MWe"/
$include ngcc_demand.gms
/

price(t) "electricity price in $/MWh"/
$include df_price.gms
/

;

Scalar
delta_T "time resolution in hrs" /%gams.user2%/

time_hor "time horizon in consideration in hrs" /%gams.user3%/

emin_cell "min energy cap of cell in Wh" /21.901/
emax_cell "max energy cap of cell in Wh" /107/

cell_mod "no. of cells in one module: 50kW/15.46W" /3234/

max_mod "max no. of battery modules in parallel" /4500/

minlf_ng "min load factor for ngcc" /0.4/

nc_eo "nominal energy cap of single battery module in kWh"

nc_eo_cell "nominal energy cap of single battery module in Wh" /107/

nc_po "nominal power o/p of single battery module in kW" /50/

eeta_acdc "a/c d/c conv efficiency, both ways" /0.95/

epsi "small value for power constraint" /0.01/

ro_ng "unit ramp rate of NGCC in %/min.MW" /5/

rr_ng "NGCC ramp rate in MW"

crf "capital recovery factor"

disc_fac "discounting factor"

r_disc "annual discount rate" /0.1/

t_lf "lifetime of plant operation in years" /15/

unders_pen "undersupply penalty in $/MWh" /1000/

cyctol "Tolerance for meeting cycling condition" /0.1/
;

nc_eo = nc_eo_cell*cell_mod*0.001;

disc_fac = (1/r_disc)*(1 - ((1/(1+r_disc))**t_lf));

crf = 1/disc_fac;

rr_ng = ro_ng*delta_t*60*ynom_ng('p_ng')/100;




Scalar

cc_ic/
$include "cc_ic.gms"
/

omfo_bat/
$include "omfo_bat.gms"
/

omvo_bat/
$include "omvo_bat.gms"
/

omvo_ng/
$include "omvo_ng.gms"
/

sc_ng/
$include "sc_ng.gms"
/

overs_pen/
$include "overs_pen.gms"
/ 

;



*Cell vars
Positive variable

u_actual(t,inp)

ecell(t) "energy cap of cell at time t in Wh"

nc_e "nominal energy capacity of battery in MWh (AC)"

nc_p "nominal power capacity of battery in MW (AC)"

yc(t) "Power charged to single cell at time t in W"
yd(t) "Power discharged from single cell at time t in W"


*NGCC vars
y_ng(t,out_ng) "Actual value of NGCC plant outputs"

png_waste(t) "power plant energy wasted at time t in MW"

dem_um(t) "unmet net demand at time t"

lic_bat "levelized inv cost of battery in $"

omf_bat "fixed O&M cost of battery in $"

omv_bat(t) "var O&M cost of battery in $/hr at time t"

omv_battot "total var costs of battery in $"

omv_ng(t) "var NGCC O&M cost in $/hr at time t"

omv_ngtot

overs_cost(t) "oversupply cost at time t in $"

tot_overc "total oversupply cost in time horizon of operation in $"

unders_cost(t) "undersupply cost at time t in $"

tot_underc "total undersupply cost in time horizon of operation in $"

rc_ng

rc_ngtot

P_grid(t)

fix_yng

fix_curr
 
;


Variable
y_cell(t,var)

df_cell(t,var,st_eq) "output var exp for single cell at time t and state st"

df1_cell(t,var,st_eq) "main summation has been split into 3 terms"
df2_cell(t,var,st_eq)
df3_cell(t,var,st_eq)

flag_c(t)
flag_d(t)
flag_1(t)
flag_2(t)

P_bat(t) "Battery power charged/discharged in MW"
tot_cost

cycpen "Penalty for cycling constraint violation"

price_cost(t)
price_costtot
;




Integer variable
mod_bat "no.of battery modules in parallel"
;


Binary variable 
z_idle(t)
;





*splitting df in 3 terms

Equation df_cell_eq_Volt; df_cell_eq_Volt(t,var,st_eq)$( ord(var) eq 2 ).. df_cell(t,var,st_eq) =e= df1_cell(t,var,st_eq) + df2_cell(t,var,st_eq) + df3_cell(t,var,st_eq);



Equation f1_cell_eq_Volt_1c; f1_cell_eq_Volt_1c(t,var)$((ord(var) eq 2) ).. df1_cell(t,var,'1c') =e= sum(z$p1_map(var,'1c',z),sum((j,k)$(n1_map(var,'1c',j) and k_map1(var,'1c',k,z,j)),(power((u_actual(t - (ord(j) - 1),'I_cell')*(-1) ),ord(z)))*theta(var,'1c',k) ));

Equation f1_cell_eq_Volt_2c; f1_cell_eq_Volt_2c(t,var)$((ord(var) eq 2) ).. df1_cell(t,var,'2c') =e= sum(z$p1_map(var,'2c',z),sum((j,k)$(n1_map(var,'2c',j) and k_map1(var,'2c',k,z,j)),(power((u_actual(t - (ord(j) - 1),'I_cell')*(-1) ),ord(z)))*theta(var,'2c',k) ));

Equation f1_cell_eq_Volt_1d; f1_cell_eq_Volt_1d(t,var)$((ord(var) eq 2) ).. df1_cell(t,var,'1d') =e= sum(z$p1_map(var,'1d',z),sum((j,k)$(n1_map(var,'1d',j) and k_map1(var,'1d',k,z,j)),(power((u_actual(t - (ord(j) - 1),'I_cell') ),ord(z)))*theta(var,'1d',k) ));

Equation f1_cell_eq_Volt_2d; f1_cell_eq_Volt_2d(t,var)$((ord(var) eq 2) ).. df1_cell(t,var,'2d') =e= sum(z$p1_map(var,'2d',z),sum((j,k)$(n1_map(var,'2d',j) and k_map1(var,'2d',k,z,j)),(power((u_actual(t - (ord(j) - 1),'I_cell') ),ord(z)))*theta(var,'2d',k) ));



Equation f2_cell_eq_Volt; f2_cell_eq_Volt(t,var,st_eq)$( (ord(var) eq 2) ).. df2_cell(t,var,st_eq) =e= sum(a$p2_map(var,st_eq,a),sum((b,k)$(n2_map(var,st_eq,b) and k_map2(var,st_eq,k,a,b)),(power((y_cell(t - (ord(b) - 1),'SOD')),ord(a)))*theta(var,st_eq,k) ));



Equation f3_cell_eq_Volt; f3_cell_eq_Volt(t,var,st_eq)$( ord(var) eq 2 ).. df3_cell(t,var,st_eq) =e= sum(h$q1_map(var,st_eq,h),sum((g,k)$(r1_map(var,st_eq,g) and k_map3(var,st_eq,k,h,g)),(power((y_cell(t-ord(g), var)),ord(h)))*theta(var,st_eq,k) ));



Equation flag_charge_eq; flag_charge_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_c(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) - yd(t) + yc(t)  ) / 2;

Equation flag_disc_eq; flag_disc_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_d(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) + yd(t) - yc(t)  ) / 2;

Equation flag_single_eq; flag_single_eq(t).. sqrt(power( (y_cell(t,'SOD') - 0.558) , 2) )*flag_1(t) =e= ( sqrt(power( (y_cell(t,'SOD') - 0.558) , 2) ) + y_cell(t,'SOD') - 0.558 ) / 2;

Equation flag_two_eq; flag_two_eq(t).. sqrt(power( (y_cell(t,'SOD') - 0.558) , 2) )*flag_2(t) =e= ( sqrt(power( (y_cell(t,'SOD') - 0.558) , 2) ) - y_cell(t,'SOD') + 0.558 ) / 2;



Equation combo_outeq; combo_outeq(t,var)$( (ord(t) ne 1) and (ord(var) eq 2) ).. y_cell(t,var) =e=  ( ( flag_1(t)*flag_c(t)*df_cell(t,var,'1c') ) + ( flag_1(t)*flag_d(t)*df_cell(t,var,'1d') ) + ( flag_2(t)*flag_c(t)*df_cell(t,var,'2c') ) + ( flag_2(t)*flag_d(t)*df_cell(t,var,'2d') ) )*(1 - z_idle(t) )  + y_cell(t-1,var)*z_idle(t);


Equation new_powereq; new_powereq(t).. sqrt( power( ( yc(t) - yd(t) ), 2 ) ) =e= y_cell(t,'V')*u_actual(t,'I_cell')*(1 - z_idle(t) );


*Equation ecellexp; ecellexp(t).. ecell(t) =e=  emax_cell + 59.501 - 162.13*y_cell(t,'SOD');

Equation ecellexp; ecellexp(t).. ecell(t) =e=  144.6 - 162.13*y_cell(t,'SOD');



Equation ebalcell; ebalcell(t)$(ord(t) ne card(t)).. ecell(t+1) =e= ecell(t) + ( yc(t) - yd(t) )*delta_T;


Equation pconn; pconn(t).. P_bat(t) =e= cell_mod*mod_bat*power(10,-6)*( yd(t) - yc(t) );


Equation energy_cap; energy_cap.. nc_e =e= nc_eo*mod_bat*0.001;

Equation pwr_cap; pwr_cap.. nc_p =e= nc_po*mod_bat*0.001;



Equation ramp_ng1; ramp_ng1(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =g= -rr_ng;

Equation ramp_ng2; ramp_ng2(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =l= rr_ng;



Equation balint; balint(t)$(ord(t) ne card(t)).. P_dem(t) =e= P_bat(t) + y_ng(t,'p_ng') + P_grid(t) - png_waste(t) + dem_um(t);



Equation batinv; batinv.. lic_bat =e= cc_ic*(nc_e**0.95)*1000*crf*time_hor/8760;


Equation batfix; batfix.. omf_bat =e= omfo_bat*nc_p*1000*time_hor/8760;


Equation batvar1; batvar1(t).. omv_bat(t) =e= omvo_bat*( yc(t) + yd(t) )*cell_mod*mod_bat*power(10,-6)*delta_T;


Equation batvar3; batvar3.. omv_battot =e= sum(t$(ord(t) ne card(t) ),omv_bat(t));


Equation ngvar; ngvar(t).. omv_ng(t) =e= omvo_ng*y_ng(t,'p_ng')*delta_T/(0.761 + 0.239*( y_ng(t,'p_ng')/ynom_ng('p_ng') ));

Equation ngvartot; ngvartot.. omv_ngtot =e= sum(t$(ord(t) ne card(t) ), omv_ng(t));




Equation ngcyc1; ngcyc1(t)$(ord(t) ne card(t)).. rc_ng(t) =g= sc_ng*(y_ng(t+1,'p_ng') - y_ng(t,'p_ng'));

Equation ngcyc2; ngcyc2(t)$(ord(t) ne card(t)).. rc_ng(t) =g= -sc_ng*(y_ng(t+1,'p_ng') - y_ng(t,'p_ng'));

Equation ngcyctot; ngcyctot.. rc_ngtot =e= sum(t$(ord(t) ne card(t)), rc_ng(t));


Equation ngovereq; ngovereq(t).. overs_cost(t) =e= png_waste(t)*overs_pen*delta_T;

Equation sumovereq; sumovereq.. tot_overc =e= sum(t$(ord(t) ne card(t) ), overs_cost(t));


Equation ngundereq; ngundereq(t).. unders_cost(t) =e= dem_um(t)*unders_pen*delta_T;

Equation sumundereq; sumundereq.. tot_underc =e= sum(t$(ord(t) ne card(t) ), unders_cost(t));


Equation objeq; objeq.. tot_cost =e= (lic_bat + omf_bat + omv_battot + omv_ngtot + tot_underc + tot_overc + rc_ngtot + price_costtot)/(1e3);



Equation ecellcy1; ecellcy1(t)$(ord(t) eq card(t)).. ( ecell(t) - ecell('1') )  =l= cyctol*ecell('1'); 

Equation ecellcy2; ecellcy2(t)$(ord(t) eq card(t)).. ( ecell(t) - ecell('1') )  =g= -1*cyctol*ecell('1'); 




Equation pricecosteq; pricecosteq(t)$(ord(t) ne card(t)).. price_cost(t) =e= P_grid(t)*price(t)*delta_T;

Equation pricecosttoteq; pricecosttoteq.. price_costtot =e= sum(t$(ord(t) ne card(t)), price_cost(t) );


Equation fix_curreq; fix_curreq(t).. u_actual(t, 'I_cell') =e= fix_curr;

Equation buylim; buylim(t).. P_grid(t) =l= (-1)*P_bat(t)*flag_c(t);

*Equation fix_yngeq; fix_yngeq(t).. y_ng(t,'p_ng') =e= fix_yng;



$include "bounds.gms";

$include "fixyc.gms";

$include "fixyd.gms";

fix_curr.fx = 6;

P_grid.fx(t) = 0;


 
mod_bat.fx = 0;

z_idle.fx(t) = 1;

model storage /all/;

option reslim = 3600;

solve storage using MINLP minimizing tot_cost;


mod_bat.lo = 0;

mod_bat.up = max_mod;

z_idle.lo(t) = 0;

z_idle.up(t) = 1;

$include "specify_reslim.gms";

solve storage using MINLP minimizing tot_cost;



*****Post processing******
execute_unload "results.gdx";
execute "gdx2sqlite -i results.gdx -o results.db";
