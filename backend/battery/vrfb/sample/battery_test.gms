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
;

Scalar
delta_T "time resolution in hrs" /%gams.user2%/
time_hor "time horizon in consideration in hrs" /%gams.user3%/
minlf_ng "min load factor for ngcc" /0.4/
ro_ng "unit ramp rate of NGCC in %/min.MW" /5/
rr_ng "NGCC ramp rate in MW"
crf "capital recovery factor"
disc_fac "discounting factor"
m1 "Parameter in linear model of power plant efficiency ratio"/0.761/
m2 "Parameter in linear model of power plant efficiency ratio"/0.239/
r_disc "annual discount rate" /0.10/
unders_pen "undersupply penalty in $/MWh" /1000/
cyctol "Tolerance for meeting cycling condition" /0.1/
t_lf "lifetime of storage operation in years" /15/
eta_rte "Round-trip efficiency of technology i (yrs)"/0.7/
Eub "Upper bound on storage capacity of technology i (MWh)"/1000/
Pub "Upper bound on discharge power capacity of technology i (MW)"/100/
s_lb "Lower bounds for state variable: here SOC"/0.3/
s_ub "Upper bounds for state variable"/1/
V_lb "lower bound on cell voltage (V)" /0.1/
V_ub "upper bound on cell voltage (V)" /1.8/
Vnom "nominal cell voltage (V)" /1.6/
Inom "nominal cell current (A)" /0.4/

c1_vrfb /-0.037322528321698/
c2_vrfb /-0.002757770570363/
c3_vrfb /0.000003656406576/
c4_vrfb /0.000003162999015/
c5_vrfb /0.775209348664613/
c6_vrfb /0.225021392375511/
ncpocell_vrfb "Power output of single cell (W)" 
Etank_vrfb "Energy of battery per unit volume of electrolyte at fully charged condition (Wh/L)" /25/
aspect_vrfb "Aspect ratio for tank sizing: h/D" /5/
rho_vrfb "average electrolyte density (g/cm3)" /1.35/
;

disc_fac = (1/r_disc)*(1 - ((1/(1+r_disc))**t_lf));
crf = 1/disc_fac;
rr_ng = ro_ng*delta_t*60*ynom_ng('p_ng')/100;
ncpocell_vrfb = Vnom*Inom;

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
Fcell_vrfb(t) "Flow for single cell (mL/hr)"
Wpump_vrfb(t) "Pump power (MW)"
Vtank_vrfb "Volume of the electrolyte tank (L)"
Dtank_vrfb "Diameter of the electrolyte tank (cm)"
htank_vrfb "Height of the electrolyte tank (cm)"
dP_vrfb "Pump pressure increase (Pa)" 
Paux(t) "total power consumed by auxiliary process equipment at time t"
lcos "storage lcos ($/MWh)"
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

**Storage energy balance equations
Equation flag_charge_eq; flag_charge_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_c(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) - yd(t) + yc(t)  ) / 2;
Equation flag_disc_eq; flag_disc_eq(t).. sqrt( power(( yd(t) - yc(t) ),2) )*flag_d(t) =e= ( sqrt( power(( yd(t) - yc(t) ),2) ) + yd(t) - yc(t)  ) / 2;
Equation maxcapeq; maxcapeq(t).. E_bat(t) =l= nc_e;
Equation eubeq; eubeq.. nc_e =l= Eub;
*Equation pbndd1; pbndd1(t).. P_bat(t) =l= nc_p;
*Equation pbndd2; pbndd2(t).. (-1)*P_bat(t) =l= nc_p;
Equation pbndd3; pbndd3.. nc_p =l= Pub;
Equation ebalcell; ebalcell(t)$(ord(t) ne card(t)).. E_bat(t+1) =e= E_bat(t) - ( eta_rte*flag_c(t) + flag_d(t) )*P_bat(t)*delta_T;
Equation ecellcy1; ecellcy1(t)$(ord(t) eq card(t)).. ( E_bat(t) - E_bat('1') )  =l= cyctol*E_bat('1'); 
Equation ecellcy2; ecellcy2(t)$(ord(t) eq card(t)).. ( E_bat(t) - E_bat('1') )  =g= -1*cyctol*E_bat('1'); 
Equation fix_curreq; fix_curreq(t).. u_actual(t, 'I_cell') =e= Inom*(1 - z_idle(t));
Equation fcellvrfbeq; fcellvrfbeq(t).. Fcell_vrfb(t) =e= 1800*(1 - z_idle(t));

**f1
Equation combo_outeq; combo_outeq(t,var)$( (ord(var) eq 2) and (ord(t) ne 1) ).. y_cell(t,var) =e=  ( c1_vrfb*( -flag_c(t) + flag_d(t) )*u_actual(t,'I_cell') + c2_vrfb*( -flag_c(t-1) + flag_d(t-1) )*u_actual(t-1,'I_cell') + c3_vrfb*Fcell_vrfb(t) + c4_vrfb*Fcell_vrfb(t-1) + c5_vrfb*y_cell(t-1,var) + c6_vrfb*y_cell(t-2,var) )*(1 - z_idle(t)) + y_cell(t-1,var)*z_idle(t);
*Equation combo_outeq; combo_outeq(t,var)$( (ord(var) eq 2) and (ord(t) ne 1) and (ord(t) ne 2) ).. y_cell(t,var) =e= c1_vrfb*( -flag_c(t) + flag_d(t) )*u_actual(t,'I_cell') + c2_vrfb*( -flag_c(t-1) + flag_d(t-1) )*u_actual(t-1,'I_cell') + c3_vrfb*Fcell_vrfb(t) + c4_vrfb*Fcell_vrfb(t-1) + c5_vrfb*y_cell(t-1,var) + c6_vrfb*y_cell(t-2,var) ;
Equation combo_outeq1; combo_outeq1.. y_cell('1','V') =e= y_cell('2','V');
Equation new_powereq; new_powereq(t).. sqrt( power( ( yc(t) - yd(t) ), 2 ) ) =e= y_cell(t,'V')*u_actual(t,'I_cell')*(1 - z_idle(t) );
Equation pconn; pconn(t).. P_bat(t) =e= ( yd(t) - yc(t) )*nc_p/ncpocell_vrfb;

Equation vrfb_vol; vrfb_vol.. Vtank_vrfb =e= nc_e*(10**6)/Etank_vrfb;
Equation vrfb_dia; vrfb_dia.. Dtank_vrfb =e= 10*( (Vtank_vrfb*4/(pi*aspect_vrfb))**(1/3) );
Equation vrfb_ht; vrfb_ht.. htank_vrfb =e= Dtank_vrfb*aspect_vrfb;
Equation vrfb_pr; vrfb_pr.. dP_vrfb =e= 10000 + (rho_vrfb*980.991*htank_vrfb/10);
Equation vrfb_pump; vrfb_pump(t).. Wpump_vrfb(t) =e= ( dP_vrfb*Fcell_vrfb(t)*(1 - z_idle(t))*(nc_p*(10**6)/ncpocell_vrfb) )/(3600*0.8*(10**12));
Equation pauxeq; pauxeq(t).. Paux(t) =e= Wpump_vrfb(t);

**f2
Equation vrfb_f2; vrfb_f2(t).. E_bat(t) =e= nc_e*y_cell(t,'SOC');

***Storage cost equations
Equation batinv; batinv.. lic_bat =e= cc_ic*(nc_e**0.95)*1000*crf*time_hor/8760;
Equation batfix; batfix.. omf_bat =e= omfo_bat*nc_p*1000*time_hor/8760;
Equation batvar1; batvar1(t).. omv_bat(t) =e= omvo_bat*( yc(t) + yd(t) )*(nc_p/ncpocell_vrfb)*delta_T;
Equation batvar3; batvar3.. omv_battot =e= sum(t$(ord(t) ne card(t) ),omv_bat(t));

***Power plant equations
Equation ramp_ng1; ramp_ng1(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =g= -rr_ng;
Equation ramp_ng2; ramp_ng2(t)$( ord(t) ne card(t) ).. y_ng(t+1,'p_ng') - y_ng(t,'p_ng') =l= rr_ng;

***Power plant cost equations
Equation ngvar; ngvar(t).. omv_ng(t) =e= omvo_ng*y_ng(t,'p_ng')*delta_T/(m1 + m2*( y_ng(t,'p_ng')/ynom_ng('p_ng') ));
Equation ngvartot; ngvartot.. omv_ngtot =e= sum(t$(ord(t) ne card(t) ), omv_ng(t));
Equation ngcyc1; ngcyc1(t)$(ord(t) ne card(t)).. rc_ng(t) =g= sc_ng*(y_ng(t+1,'p_ng') - y_ng(t,'p_ng'));
Equation ngcyc2; ngcyc2(t)$(ord(t) ne card(t)).. rc_ng(t) =g= -sc_ng*(y_ng(t+1,'p_ng') - y_ng(t,'p_ng'));
Equation ngcyctot; ngcyctot.. rc_ngtot =e= sum(t$(ord(t) ne card(t)), rc_ng(t));

***Overall energy balance and cost eqs
Equation balint; balint(t)$(ord(t) ne card(t)).. P_dem(t) =e= P_bat(t) + y_ng(t,'p_ng') - png_waste(t) + dem_um(t) - Paux(t);
Equation ngovereq; ngovereq(t).. overs_cost(t) =e= png_waste(t)*overs_pen*delta_T;
Equation sumovereq; sumovereq.. tot_overc =e= sum(t$(ord(t) ne card(t) ), overs_cost(t));
Equation ngundereq; ngundereq(t).. unders_cost(t) =e= dem_um(t)*unders_pen*delta_T;
Equation sumundereq; sumundereq.. tot_underc =e= sum(t$(ord(t) ne card(t) ), unders_cost(t));
Equation objeq; objeq.. tot_cost =e= (lic_bat + omf_bat + omv_battot + omv_ngtot + tot_underc + tot_overc + rc_ngtot)/(100);


$include "bounds.gms";
$include "fixyc.gms";
$include "fixyd.gms";

nc_p.fx = 0;
z_idle.fx(t) = 1;
model storage /all/;
option reslim = 3600;
solve storage using MINLP minimizing tot_cost;


nc_p.lo = 0;
nc_p.up = Pub;
z_idle.lo(t) = 0;
z_idle.up(t) = 1;
z_idle.fx('1') = 1;
z_idle.fx('2') = 1;
$include "voltage_iniguess.gms";
$include "specify_reslim.gms";
solve storage using MINLP minimizing tot_cost;

*****Post processing******
execute_unload "results.gdx";
execute "gdx2sqlite -i results.gdx -o results.db";