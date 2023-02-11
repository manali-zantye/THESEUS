option minlp = baron;
option optcr = 1e-3;
option optca = 1e-3;

Sets
t "decision stages: stage 1 is initialization, not decision" /1*%gams.user1%/
var "cell output variables"/SOC, V/
inp "cell input variables: current(A)" /I_cell/
out_ng "NGCC output vars" /p_ng/
;

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
emax_cell "max energy cap of cell in Wh" /87.6/
minlf_ng "min load factor for ngcc" /0.4/
ro_ng "unit ramp rate of NGCC in %/min.MW" /5/
rr_ng "NGCC ramp rate in MW"
crf "capital recovery factor"
disc_fac "discounting factor"
r_disc "annual discount rate" /0.1/
t_lf "lifetime of plant operation in years" /6/
unders_pen "undersupply penalty in $/MWh" /1000/
cyctol "Tolerance for meeting cycling condition" /0.1/
Rcell_lib "Single cell internal resistance (Ohms)" /0.0082/
a1_lib /0.82/
a2_lib /-0.72/
a3_lib /0.69/
a4_lib /3.45/
ncpocell_lib "Power output of single cell (W): will update based on current value" /21.9/
ncho_lib "Max storage discharge hours (hr) for LiB" /4/
eta_rte "Round-trip efficiency of technology i (yrs)"/0.95/
;

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


Positive variable
u_actual(t,inp)
yc(t) "Power charged to single cell at time t in W"
yd(t) "Power discharged from single cell at time t in W"
OCV_lib(t) "Open circuit voltage in V"
y_cell(t,var)
ncell_lib
nc_e "nominal energy capacity of battery in MWh (AC)"
nc_p "nominal power capacity of battery in MW (AC)"
E_bat(t)
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
fix_curr
E_cell(t) "Energy capacity of each cell (Wh)"
;

Variable
flag_c(t)
flag_d(t)
P_bat(t) "Battery power charged/discharged in MW"
tot_cost
;

Binary variable
z_idle(t)
;

Equation flag_charge_eq; flag_charge_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_c(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) - yd(t) + yc(t)  ) / 2;
Equation flag_disc_eq; flag_disc_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_d(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) + yd(t) - yc(t)  ) / 2;

Equation combo_outeq; combo_outeq(t,var)$( ord(var) eq 2 ).. y_cell(t,var) =e=  OCV_lib(t) - u_actual(t,'I_cell')*Rcell_lib*( -flag_c(t) + flag_d(t) )*(1 - z_idle(t) );
Equation ocveq_lib; ocveq_lib(t).. OCV_lib(t) =e= a1_lib*(y_cell(t,'SOC')**3) + a2_lib*(y_cell(t,'SOC')**2) + a3_lib*y_cell(t,'SOC') + a4_lib;
Equation new_powereq; new_powereq(t).. sqrt( power( ( yc(t) - yd(t) ), 2 ) ) =e= y_cell(t,'V')*u_actual(t,'I_cell')*(1 - z_idle(t) );
Equation pconn; pconn(t).. P_bat(t) =e= ncell_lib*power(10,-6)*( yd(t) - yc(t) );
Equation ncelleq_lib; ncelleq_lib.. ncell_lib =e= nc_p*(10**6)/25.5;

Equation epratio_lib; epratio_lib.. nc_e =e= ncho_lib*nc_p;
Equation lib_f2; lib_f2(t).. E_bat(t) =e= nc_e*y_cell(t,'SOC');
Equation ebalcell; ebalcell(t)$(ord(t) ne card(t)).. E_bat(t+1) =e= E_bat(t) - ( eta_rte*flag_c(t) + flag_d(t) )*P_bat(t)*delta_T;
Equation ecellcy1; ecellcy1(t)$(ord(t) eq card(t)).. ( E_bat(t) - E_bat('1') )  =l= cyctol*E_bat('1'); 
Equation ecellcy2; ecellcy2(t)$(ord(t) eq card(t)).. ( E_bat(t) - E_bat('1') )  =g= -1*cyctol*E_bat('1'); 
Equation fix_curreq; fix_curreq(t).. u_actual(t, 'I_cell') =e= fix_curr;

Equation maxcapeq; maxcapeq(t).. E_bat(t) =l= nc_e;
Equation pbndd1; pbndd1(t).. P_bat(t) =l= nc_p;
Equation pbndd2; pbndd2(t).. (-1)*P_bat(t) =l= nc_p;


Equation ramp_ng1; ramp_ng1(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =g= -rr_ng;
Equation ramp_ng2; ramp_ng2(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =l= rr_ng;
Equation balint; balint(t)$(ord(t) ne card(t)).. P_dem(t) =e= P_bat(t) + y_ng(t,'p_ng') - png_waste(t) + dem_um(t);
Equation batinv; batinv.. lic_bat =e= cc_ic*(nc_e**0.95)*1000*crf*time_hor/8760;
Equation batfix; batfix.. omf_bat =e= omfo_bat*nc_p*1000*time_hor/8760;
Equation batvar1; batvar1(t).. omv_bat(t) =e= omvo_bat*( yc(t) + yd(t) )*ncell_lib*power(10,-6)*delta_T;
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
Equation objeq; objeq.. tot_cost =e= (lic_bat + omf_bat + omv_battot + omv_ngtot + tot_underc + tot_overc + rc_ngtot)/(1e3);

$include "bounds.gms";
$include "fixyc.gms";
$include "fixyd.gms";
fix_curr.fx = 6;

nc_p.fx = 0;
z_idle.fx(t) = 1;
model storage /all/;
option reslim = 3600;
solve storage using MINLP minimizing tot_cost;

nc_p.lo = 0;
nc_p.up = 200;
z_idle.lo(t) = 0;
z_idle.up(t) = 1;
$include "specify_reslim.gms";
solve storage using MINLP minimizing tot_cost;

*****Post processing******
execute_unload "results.gdx";
execute "gdx2sqlite -i results.gdx -o results.db";