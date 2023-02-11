option minlp = baron;
option optcr = 1e-4;
option optca = 1e-4;

Sets
t "Time stages" /1*%gams.user1%/
i "storage technologies" /ces, htts, caes, h2, phs, pcm, nas, vrfb, lib/
b "storage state of operation" /c, d, idle/
;

Scalars
delta_t "time resolution in hrs" /%gams.user2%/
time_hor "time horizon in consideration in hrs" /24/
r_disc "annual discount rate" /0.1/
overs_pen "Oversupply penalty in $/MWh"/155/ 
unders_pen "Undersupply penalty in $/MWh"/1000/ 
rr_fp "fossil plant ramp rate in MW"
cyctol "Tolerance for meeting storage cycling condition" /0/
epsi_lb "Lower bound on storage size is selected (MWh)"/0.1/
beta_scal "Scaling factor for economies of numbers"/0.85/
;

$include "pars/fp_pars.gms";
rr_fp = ro_fp*delta_t*60*Pnom_fp/100;
$include "pars/stor_pars.gms";

Parameters
disc_fac(i) "discounting factor for technology i"
crf(i) "capital recovery factor for technology i"
P_dem(t) "power demand at time t for grid in MW"/
$include pars/ngcc_demand.gms
/
;


*for co2 capture 
Scalars
gamma_co2 "co2 capture efficiency" /0.9/
co_co2 "co2 capture spec cap cost ($/MW)" /810000/
eco2_ref "reference emission intensity (ton/MWh)" /0.76/
i17 "cepci" /257.5/
i02 "cepci" /128.7/
mu_abs "energy penalty of absorption" /0.02/
mu_des "energy penalty of desorption" /0.04/
mu_comp "energy penalty of compression" /0.02/
ra_max "max absorption rate" /1/
rd_max "max desoprtion rate" /1/
cap_rate "regulatory capture rate" /0/
tlf_co2 "lifetime of co2 capture plant in years" /20/
df_co2 "discounting factor for co2 capture"
crf_co2 "capital recovery factor for co2 capture"
ro_ra "max ramp rate of CO2 absorption system (1/h)" /1/
ro_rd "max ramp rate of CO2 desorption system (1/h)" /1/
;


disc_fac(i) = (1/r_disc)*( 1 - ( (1/(1+r_disc))**tlf(i) ) );
crf(i) = 1/disc_fac(i);
df_co2 = (1/r_disc)*( 1 - ( (1/(1+r_disc))**tlf_co2 ) );
crf_co2 = 1/df_co2;

Positive variables
*dummy
Pd_max(i) "Decision: Max discharging power output of technology i in MW"
x(i) "Decision:Storage size of technology i in MWh"
fl(i, t) "Decision: Flow variable for technology i at time t"
P_fp(t) "Decision: Power op of fossil plant at time t (MW)"
Pc_max(i) "Decision: Max charging power output of technology i in MW"
st(i, t) "State of technology i at time t"
E(i, t) "Energy capacity of technology i at time t (MWh)"
eff_fp(t) "Efficiency of fossil plant at time t"
cov_fp(t) "Fossil plant variable o&m cost ($)"
crc_fp(t) "Fossil plant cycling cost ($)"
cov_stor(i,t) "Storage variable o&m cost ($)"
cov_fptot "Total Fossil plant variable o&m cost ($)"
crc_fptot "Total Fossil plant cycling cost ($)"
cov_stortot(i) "Total Storage variable o&m cost ($)"
civ_stor(i) "Storage investment cost scaled over scheduling horizon ($)"
cof_stor(i) "Storage fixed o&m cost scaled over scheduling horizon ($)"
pfp_waste(t) "Oversupply of electricity to grid at time t in MW"
dem_um(t) "undersupply of electricity at time t"
overs_cost(t) "Oversupply cost at time t"
tot_overc "Total oversupply cost over scheduling horizon"
unders_cost(t) "Undersupply cost at time t"
tot_underc "Total undersupply cost over scheduling horizon"
fp_effratio(t) "Ratio of actual efficiency to nominal efficiency for fossil plant"
lcos "Levelized cost of storage"

*co2 capture variables
mco2_flue(t) "co2 flowrate in flue gas (ton/hr)" 
mco2_capt(t) "co2 capture flowrate (ton/hr)"
mco2_em(t) "co2 flowrate in emissions to environment (ton/hr)"
tax_co2(t) "co2 emission tax paid ($)"
tax_co2tot "Total co2 emission tax over scheduling horizon ($)"
rev_co2(t) "revenue obtained from co2 sale ($)"
rev_co2tot "Total revenue obtained over scheduling horizon ($)"
P_co2(t) "Power consumed by co2 capture system (MW)"
cc_co2 "capital cost of co2 capture system ($)"
ra_co2(t) "co2 capture absorption rate"
rd_co2(t) "co2 capture desorption rate"
;

Binary variable
y(i) "Decision: Selection of technology i"
z_op(i,t,b) "Decision: State of operation of technology i at time t"

*CO2 capture binary variables
co2_sel "for co2 capture selection"
;

Integer variable
train(i) "No. of parallel trains of technology i";

Variable 
*obj
tot_cost "Total system cost over operating horizon (1000 $)"
P(i,t) "Power o/p of technology i at time t (MW): positive for discharge, negative for charge"
f1(i, t) "Set of eqns for power capacity of technology i at time t: function of state st(i,t), flow fl(i,t), operation z_op(i,t,b) and size x(i)"
f2(i, t) "Set of eqns for energy capacity of technology i at time t: function of state st(i,t), and size x(i)"

*CO2 capture
co2cost_tot "Total co2 capture cost including co2 tax, co2 revenue and capital cost ($)"
;

*Equation dumeq; dumeq.. obj =e= dummy;

***Storage energy balance equations
Equation sumsteq; sumsteq(i,t).. sum(b, z_op(i,t,b) ) =e= 1;
Equation eubeq; eubeq(i).. x(i) =l= Eub(i)*y(i);
Equation elbeq; elbeq(i).. x(i) =g= epsi_lb*y(i);
Equation maxcapeq; maxcapeq(i, t).. E(i, t) =l= x(i);
Equation energybal; energybal(i, t)$(ord(t) ne card(t)).. E(i, t + 1) =e= E(i, t) -  (eta_rte(i)*(st('caes',t)**exp_rte(i))*z_op(i,t,'c') + z_op(i,t,'d'))*P(i, t)*delta_t;
Equation pbndc1; pbndc1(i,t).. (-1)*z_op(i,t,'c')*P(i,t) =l= Pc_max(i);
Equation pbndc2; pbndc2(i).. Pc_max(i) =l= Pc_ub(i)*y(i);
Equation pbndd1; pbndd1(i,t).. z_op(i,t,'d')*P(i,t) =l= Pd_max(i);
Equation pbndd2; pbndd2(i).. Pd_max(i) =l= Pd_ub(i)*y(i);
Equation ecellcy1; ecellcy1(i, t)$(ord(t) eq card(t)).. ( E(i, t) - E(i, '1') )  =l= cyctol*E(i, '1'); 
Equation ecellcy2; ecellcy2(i, t)$(ord(t) eq card(t)).. ( E(i, t) - E(i, '1') )  =g= -1*cyctol*E(i, '1'); 
Equation storramp1; storramp1(i,t)$( ord(t) ne card(t) ).. ( z_op(i,t + 1,'d') - z_op(i,t + 1,'c') )*P(i, t + 1) - ( z_op(i,t,'d') - z_op(i,t,'c') )*P(i, t) =l= ro_stor(i)*delta_T*60;
Equation storramp2; storramp2(i,t)$( ord(t) ne card(t) ).. ( z_op(i,t + 1,'d') - z_op(i,t + 1,'c') )*P(i, t + 1) - ( z_op(i,t,'d') - z_op(i,t,'c') )*P(i, t) =g= -ro_stor(i)*delta_T*60;
Equation train_bnd1; train_bnd1(i).. y(i)*1 =l= train(i);
Equation train_bnd2; train_bnd2(i).. train(i) =l= y(i)*5;

**Storage Technology specific models
Equation pspec; pspec(i, t).. P(i, t) =e= f1(i, t);
Equation espec; espec(i, t).. E(i, t) =e= f2(i, t);
$include "pars/f1.gms";
$include "pars/f2.gms";

***Storage cost equations
Equation costeq1; costeq1(i).. civ_stor(i) =e= ( c11(i)*( x(i)**a11(i) ) + c12(i)*( Pd_max(i)**a12(i) ) + c13(i)*( Pc_max(i)**a13(i) ) )*crf(i)*time_hor/8760;
Equation costeq2; costeq2(i).. cof_stor(i) =e= ( c21(i)*( x(i)**a21(i) ) + c22(i)*( Pd_max(i)**a22(i) ) + c23(i)*( Pc_max(i)**a23(i) ) )*time_hor/8760;
Equation costeq31; costeq31(i, t).. cov_stor(i,t) =e= (c31(i)*z_op(i,t,'d')*P(i,t) - c32(i)*z_op(i,t,'c')*P(i,t))*delta_t; 
Equation costeq32; costeq32(i).. cov_stortot(i) =e= sum(t$(ord(t) ne card(t) ), cov_stor(i,t) );

***Power plant equations
Equation ramp_fp1; ramp_fp1(t)$( ord(t) ne card(t) ).. P_fp(t+1) - P_fp(t) =g= -rr_fp;
Equation ramp_fp2; ramp_fp2(t)$( ord(t) ne card(t) ).. P_fp(t+1) - P_fp(t) =l= rr_fp;
Equation fpeffeq; fpeffeq(t).. fp_effratio(t) =e= m1 + m2*( P_fp(t) / Pnom_fp);

***Power plant cost equations
Equation costeq41; costeq41(t).. cov_fp(t)*fp_effratio(t) =e= c4*P_fp(t)*delta_t;
Equation costeq42; costeq42.. cov_fptot =e= sum(t$(ord(t) ne card(t) ), cov_fp(t) );
Equation costeq51; costeq51(t)$( ord(t) ne card(t) ).. crc_fp(t) =g= c5*(P_fp(t+1) - P_fp(t));
Equation costeq52; costeq52(t)$(ord(t) ne card(t)).. crc_fp(t) =g= -c5*(P_fp(t+1) - P_fp(t));
Equation costeq53; costeq53.. crc_fptot =e= sum( t$(ord(t) ne card(t)), crc_fp(t) );

***Overall energy balance and cost eqs
Equation balint; balint(t)$(ord(t) ne card(t)).. P_dem(t) =e= sum(i, P(i, t)*train(i) ) + P_fp(t) - pfp_waste(t) + dem_um(t) - P_co2(t);
Equation overeq; overeq(t).. overs_cost(t) =e= pfp_waste(t)*overs_pen*delta_T;
Equation sumovereq; sumovereq.. tot_overc =e= sum(t$(ord(t) ne card(t) ), overs_cost(t));
Equation undereq; undereq(t).. unders_cost(t) =e= dem_um(t)*unders_pen*delta_T;
Equation sumundereq; sumundereq.. tot_underc =e= sum(t$(ord(t) ne card(t) ), unders_cost(t));

***co2 system equations
***Mass balance
Equation mco2flueeq; mco2flueeq(t).. mco2_flue(t) =e= eco2*P_fp(t);
Equation mco2capteq; mco2capteq(t).. mco2_capt(t) =e= rd_co2(t)*gamma_co2*mco2_flue(t);
Equation mco2emeq; mco2emeq(t).. mco2_em(t) =e= mco2_flue(t)*(1 - ra_co2(t)*gamma_co2);
Equation mco2bal; mco2bal.. sum(t$(ord(t) ne card(t) ), mco2_flue(t) ) =e= sum(t$(ord(t) ne card(t) ), mco2_capt(t) ) + sum(t$(ord(t) ne card(t) ), mco2_em(t) );
Equation rabnd; rabnd(t).. ra_co2(t) =l= co2_sel*ra_max;
Equation rdbnd; rdbnd(t).. rd_co2(t) =l= co2_sel*rd_max;
Equation ramp_ra1; ramp_ra1(t)$( ord(t) ne card(t) ).. ra_co2(t+1) - ra_co2(t) =g= -ro_ra*delta_T;
Equation ramp_ra2; ramp_ra2(t)$( ord(t) ne card(t) ).. ra_co2(t+1) - ra_co2(t) =l= ro_ra*delta_T;
Equation ramp_rd1; ramp_rd1(t)$( ord(t) ne card(t) ).. rd_co2(t+1) - rd_co2(t) =g= -ro_rd*delta_T;
Equation ramp_rd2; ramp_rd2(t)$( ord(t) ne card(t) ).. rd_co2(t+1) - rd_co2(t) =l= ro_rd*delta_T;
Equation caprateeq; caprateeq.. sum(t$(ord(t) ne card(t) ), mco2_capt(t) ) =g= sum(t$(ord(t) ne card(t) ), mco2_flue(t) )*cap_rate;

**Energy balance
Equation powerco2eq; powerco2eq(t).. P_co2(t) =e= mu_abs*(Pnom_fp/eff_fpmax)*ra_co2(t) + (mu_des + mu_comp)*(Pnom_fp/eff_fpmax)*rd_co2(t);

**Cost equations
Equation co2taxeq1; co2taxeq1(t).. tax_co2(t) =e= ctax*mco2_em(t)*delta_T;
Equation co2taxeq2; co2taxeq2.. tax_co2tot =e= sum(t$(ord(t) ne card(t)), tax_co2(t) );
Equation co2price1; co2price1(t).. rev_co2(t) =e= cprice*mco2_capt(t)*delta_T;
Equation co2price2; co2price2.. rev_co2tot =e= sum(t$(ord(t) ne card(t)), rev_co2(t) );
Equation capcostco2eq; capcostco2eq.. cc_co2 =e= (co_co2*(i17/i02)*Pnom_fp*eco2*co2_sel/eco2_ref)*crf_co2*time_hor/8760;
Equation totco2costeq; totco2costeq.. co2cost_tot =e= tax_co2tot - rev_co2tot + cc_co2;

****Objective and LCOS
Equation objeq1; objeq1.. tot_cost =e= ( sum(i, civ_stor(i)*(train(i)**beta_scal) ) + sum(i, cof_stor(i)*train(i) ) + sum(i, cov_stortot(i)*train(i) ) + cov_fptot + tot_underc + tot_overc + crc_fptot + co2cost_tot)/(1e5);
Equation objeq2; objeq2(i).. lcos(i)*( sum(t$(ord(t) ne card(t) ),  z_op(i,t,'d')*P(i,t)*delta_t )  )*train(i) =e= ( civ_stor(i)*(train(i)**beta_scal) + cof_stor(i)*train(i) + cov_stortot(i)*train(i) +  c4*sum( t$(ord(t) ne card(t)),  (-1)*z_op(i,t,'c')*P(i,t)*delta_t )*train(i) );
Equation lcos_ub; lcos_ub(i).. lcos(i) =l= 10000*y(i);

co2_sel.fx = 0;

$include "pars/fixzc.gms";
$include "pars/fixzd.gms";
$include "bounds.gms";
y.fx(i) = 0;
z_op.fx(i,t,'idle') = 1;
model storage /all/;
option reslim = 3600;
solve storage using MINLP minimizing tot_cost;

y.lo(i) = 0;
y.up(i) = 1;
z_op.lo(i,t,'idle') = 0;
z_op.up(i,t,'idle') = 1;
$include "tech_fx.gms";
y.fx('h2') = 0;
option reslim = 400;
solve storage using MINLP minimizing tot_cost;

y.lo('h2') = 0;
y.up('h2') = 1;
$include "tech_fx.gms";
$include "specify_reslim.gms";
solve storage using MINLP minimizing tot_cost;


*****Post processing******
execute_unload "results.gdx";
execute "gdx2sqlite -i results.gdx -o results.db";