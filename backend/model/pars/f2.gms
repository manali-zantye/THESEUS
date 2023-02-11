Equation ces_f2; ces_f2(t).. f2('ces', t) =e= 1.138*0.0001*st('ces',t);

Equation htts_f2; htts_f2(t)..f2('htts', t) =e= st('htts', t)*W_d*R_f/3600;

Equation caes_f2; caes_f2(t).. f2('caes', t) =e= 71.0429*0.85*(0.0624*(log(st('caes', t)))*st('caes', t) - 0.0196*st('caes', t) - 0.2676) - 0.00115172;

Scalar nch_phs "max storage duration for phs (hrs)" /8/;
Equation nce_ncpphs; nce_ncpphs.. x('phs') =l= nch_phs*Pd_max('phs');
Equation phs_f2; phs_f2(t).. f2('phs', t) =e= st('phs', t);

Equation h2_f2; h2_f2(t).. f2('h2', t) =e= st('h2',t)*0.001*(nc_dis_h2 + (nc_pocomp_h2/1.95));

Equation pcm_f2; pcm_f2(t).. f2('pcm', t) =e= st('pcm',t);
