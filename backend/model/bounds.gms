st.lo(i, t) = s_lb(i);
st.up(i, t) = s_ub(i);

fl.lo(i, t) = l_lb(i);
fl.up(i,t) = l_ub(i);

P_fp.lo(t) = minlf_fp*Pnom_fp;
P_fp.up(t) = Pnom_fp;

train.lo(i) = 0;
train.up(i) = 5;
