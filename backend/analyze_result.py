import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import os
import math
path = os.getcwd()

tech_arr = ['ces', 'htts', 'caes', 'h2', 'phs', 'pcm', 'nas', 'vrfb', 'lib']


def tech_select(stor_inp, tech_arr): # Function to specify which technologies to fix as 0
    f = open("model/tech_fx.gms", "w+")
    for i in range(len(stor_inp)):
        if stor_inp[i] == 0:
            f.write("y.fx('%s') = 0;\n" %(tech_arr[i]))
    f.close()
    return

def fixycyd(): ## Function to fix yc and yd as zero depending on demand
    ngcc_dem = pd.read_csv('model/pars/ngcc_demand.gms', header=None, delimiter=' ', index_col=0)
    ngcc_dem.index.rename('time', inplace=True)
    ngcc_dem.rename(columns={1: 'Demand'}, inplace=True)

    for i in ngcc_dem.index:
        if ngcc_dem.loc[i, 'Demand'] <= ngcc_dem.mean()['Demand']:
            ngcc_dem.loc[i, 'flag'] = 0
        else:
            ngcc_dem.loc[i, 'flag'] = 1

    for i in ngcc_dem.index:
        if ngcc_dem.loc[i, 'flag'] == 0:
            ngcc_dem.loc[i, 'yd'] = 0
        else:
            ngcc_dem.loc[i, 'yc'] = 0

    f = open('model/pars/fixzd.gms', 'w+') 
    for j in ngcc_dem['yd'].dropna().index:    
        f.write("z_op.fx(i, '%d', 'd') = 0;\n" % (j) )
    f.close()

    f = open('model/pars/fixzc.gms', 'w+') 
    for k in ngcc_dem['yc'].dropna().index:    
        f.write("z_op.fx(i, '%d', 'c') = 0;\n" % (k) )
    f.close()

    return

def gmscall(): ## Function to call gams based on input file
    ngcc_dem = pd.read_csv('model/pars/ngcc_demand.gms', header=None, delimiter=' ', index_col=0)
    inp1 = len(ngcc_dem) + 1
    inp2 = 24/(len(ngcc_dem))
    os.system('cd model;gams main.gms u1 = ' + str(inp1) + ' u2 = ' + str(inp2) + ';cd ..;')
    return

def plot_intop(): # # Function to plot integrated system operation
    
    #read output
    results_file = 'results.db'
    database = 'model/' + results_file
    dat = sqlite3.connect(database)
    
    scal_var_gams = pd.read_sql_query("SELECT * FROM scalarvariables", dat).set_index('name').level
    scalars = pd.read_sql_query("SELECT * FROM scalars", dat).set_index('name')

    train = pd.read_sql_query("SELECT * FROM train", dat).set_index('i', drop=True).level

    P_stor_all = pd.DataFrame(pd.read_sql_query("SELECT * FROM P", dat)[['i','t','level']])
    for tech_ind in train.index:
        for k in P_stor_all.index:
            if P_stor_all.loc[k,'i'] == tech_ind:
                P_stor_all.loc[k,'level'] = P_stor_all.loc[k,'level']*train.loc[tech_ind]

    P_stor_all.set_index(['i','t'], drop=True, inplace=True)
    P_stor_sum = P_stor_all.groupby(by='t').sum()
    P_stor_sum.index = P_stor_sum.index.astype('int64')
    P_stor_sum.sort_index(ascending=True, inplace=True)
    P_stor = P_stor_sum[:-1].level

    P_fp = pd.read_sql_query("SELECT * FROM P_fp", dat).set_index('t', drop=True).level[:-1]
    P_fp.index = P_fp.index.astype('int64')

    P_co2 = pd.read_sql_query("SELECT * FROM P_co2", dat).set_index('t', drop=True).level[:-1]
    P_co2.index = P_co2.index.astype('int64')

    pfp_waste = pd.read_sql_query("SELECT * FROM pfp_waste", dat).set_index('t', drop=True).level[:-1]
    pfp_waste.index = pfp_waste.index.astype('int64')

    dem_um = pd.read_sql_query("SELECT * FROM dem_um", dat).set_index('t', drop=True).level[:-1]
    dem_um.index = dem_um.index.astype('int64')

    net_dem = pd.read_sql_query("SELECT * FROM P_dem", dat).set_index('t', drop=True).value
    net_dem.index = net_dem.index.astype('int64')

    yop = pd.DataFrame({'Pstor': P_stor, 'Pfp': P_fp, 'Pco2': (-1)*P_co2,'Pos': pfp_waste, 'Pus': dem_um, 'Pdem':net_dem})

    yop['Pd'] = yop['Pstor']
    yop.loc[yop.loc[yop['Pd'] < 0].index, 'Pd'] = 0

    yop['Pc'] = yop['Pstor']
    yop.loc[yop.loc[yop['Pc'] > 0].index, 'Pc'] = 0
    
    #Plot overall system operation
    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 25})
    p_index = np.array(yop.index).copy()
    width = 0.8

    plt.bar(p_index, 
        yop.loc[:,'Pfp'], 
        width, 
        label='NGCC net output',
        color='tab:blue')

    plt.bar(p_index, 
            yop.loc[:,'Pd'], 
            width,
            bottom = yop.loc[:,'Pfp'], 
            label='Storage discharge',
            color='tab:orange')

    plt.bar(p_index, 
            yop.loc[:,'Pc'], 
            width,
            bottom = yop.loc[:,'Pfp'], 
            label='Storage charge',
            color='tab:purple')

    plt.bar(p_index, 
        yop.loc[:,'Pco2'], 
        width, 
            
        bottom = yop.loc[:,'Pfp'] + yop.loc[:,'Pc'], 
        label='CO$_2$ capture consumption',
        color='tab:green')

    plt.bar(p_index, 
            yop.loc[:,'Pus'], 
            width, 
            bottom = yop.loc[:,'Pfp'] + yop.loc[:,'Pd'], 
            label='Undersupply',
            color='tab:red')

    plt.plot(p_index, 
             yop.loc[:,'Pdem'], 
             label='Net demand', 
             linewidth=3, 
             color='xkcd:dark green')

    plt.plot(p_index, 
             yop.loc[:,'Pfp'], 
             label='NGCC gross output', 
             linewidth=3, 
             color='xkcd:dark red')

    plt.xlabel('Time (hrs)')
    plt.legend(ncol=2, prop={'size':15}, columnspacing=0.3)
    plt.ylabel('Power (MW)')
    plt.ylim([0, net_dem.max()+50])
    plt.margins(x=0)
    plt.xticks( np.arange(0, len(np.array(yop.index)), math.ceil(len(np.array(yop.index)) / 24) ), np.arange(0, 24, 1) ); 
    plt.savefig("results/op_profile.png", dpi=300, bbox_inches='tight')
        
    return

def storage_design(): # Function to get optimal storage size

    results_file = 'results.db'
    database = 'model/' + results_file
    dat = sqlite3.connect(database)
    scal_var_gams = pd.read_sql_query("SELECT * FROM scalarvariables", dat).set_index('name').level
    y = pd.read_sql_query("SELECT * FROM y", dat).set_index('i', drop=True).level
    x = pd.read_sql_query("SELECT * FROM x", dat).set_index('i', drop=True).level
    Pd_max = pd.read_sql_query("SELECT * FROM Pd_max", dat).set_index('i', drop=True).level
    train = pd.read_sql_query("SELECT * FROM train", dat).set_index('i', drop=True).level
    lcos = (pd.read_sql_query("SELECT * FROM lcos", dat).set_index('lcos', drop=True).level)
    stor_des = pd.DataFrame({'y': y, 'x': x, 'Pd_max': Pd_max, 'train': train, 'lcos': lcos})
    stor_des.to_csv('results/stor_des.csv')

    mco2_flue = pd.read_sql_query("SELECT * FROM mco2_flue", dat).set_index('t', drop=True).level[:-1]
    mco2_flue.index = mco2_flue.index.astype('int64')
    mco2_capt = pd.read_sql_query("SELECT * FROM mco2_capt", dat).set_index('t', drop=True).level[:-1]
    mco2_capt.index = mco2_capt.index.astype('int64')
    co2_desgn = pd.DataFrame({'co2sel': scal_var_gams.loc['co2_sel'], 'co2cost_tot': scal_var_gams.loc['co2cost_tot'], 
                        'cc_co2': scal_var_gams.loc['cc_co2'], 'tax_co2tot': scal_var_gams.loc['tax_co2tot'], 
                        'rev_co2tot': scal_var_gams.loc['rev_co2tot'], 'per_capt': (mco2_capt.sum())*100/(mco2_flue.sum()) }, index=[0])
    co2_desgn.to_csv('results/co2_des.csv')
    
    return

def stor_op(tech_arr): # Function to plot storage system operation for one train

    results_file = 'results.db'
    database = 'model/' + results_file
    dat = sqlite3.connect(database)
    scalars = pd.read_sql_query("SELECT * FROM scalars", dat).set_index('name')
    l = pd.DataFrame(pd.read_sql_query("SELECT * FROM fl", dat).set_index(['i','t'], drop=True).level)
    s = pd.DataFrame(pd.read_sql_query("SELECT * FROM st", dat).set_index(['i','t'], drop=True).level)
    z_op = pd.DataFrame(pd.read_sql_query("SELECT * FROM z_op", dat).set_index(['i','t'], drop=True))
    z_c = pd.DataFrame(z_op.loc[z_op['b'] == 'c'].level)
    z_d = pd.DataFrame(z_op.loc[z_op['b'] == 'd'].level)
    z_idle = pd.DataFrame(z_op.loc[z_op['b'] == 'idle'].level)
    df_pltdet = pd.DataFrame(index=tech_arr, columns=['smult','lmult','slab', 'llab'])
    df_pltdet.loc['ces'] = [1/907.185, 1/3600, 'Mass of air stored (ton)', 'Air mass flowrate (kg/s)']
    df_pltdet.loc['htts'] = [1/907.185, 1, 'Mass of molten salt in hot tank (ton)', 'Molten salt flowrate (kg/s)']
    df_pltdet.loc['pcm'] = [1, 1, 'Energy stored (MWh)', 'Heat transfer fluid flowrate (kg/s)']
    df_pltdet.loc['caes'] = [1, 1, 'Storage pressure (bar)', 'Air mass flowrate (kg/s)']
    df_pltdet.loc['phs'] = [1, 1, 'Storage energy (MWh)', 'Storage charge/discharge rate (MW)']
    df_pltdet.loc['h2'] = [1, 1, 'Mass of H$_2$ stored (kg)', 'H$_2$ mass flowrate (kg/hr)']

    for i in tech_arr:
        l_new = l.loc[i][:-1]
        l_new.index = l_new.index.astype('int64')

        s_new = s.loc[i][:-1]
        s_new.index = s_new.index.astype('int64')

        z_op_new = z_op.loc[i][:-1]
        z_op_new.index = z_op_new.index.astype('int64')

        zc_new = z_c.loc[i][:-1]
        zc_new.index = zc_new.index.astype('int64')

        zd_new = z_d.loc[i][:-1]
        zd_new.index = zd_new.index.astype('int64')

        zidle_new = z_idle.loc[i][:-1]
        zidle_new.index = zidle_new.index.astype('int64')

        stor_op = pd.DataFrame({'l': l_new.level, 's': s_new.level, 'z_c': zc_new.level, 'z_d': zd_new.level, 
                                'z_idle': zidle_new.level})


        plt.figure(figsize=(15, 10))
        plt.rcParams.update({'font.size': 25})
        p_index = np.array(stor_op.index).copy()
        width = 0.8

        plt.bar(p_index, 
                stor_op['s']*df_pltdet.loc[i,'smult'], 
                width, 
                color='tab:green')

        plt.xlabel('Time (hrs)')
        plt.ylabel(df_pltdet.loc[i,'slab'])
        plt.twinx()
        plt.plot(p_index, (stor_op['l']*(1 - stor_op['z_idle'])*df_pltdet.loc[i,'lmult'] ), linewidth=2, color='red', label='Charge')

        plt.yticks(color='red')
        plt.ylabel(df_pltdet.loc[i,'llab'], color='red')
        plt.margins(x=0)
        plt.xticks( np.arange(0, len(np.array(stor_op.index)), math.ceil(len(np.array(stor_op.index)) / 24) ) , np.arange(0, 24, 1) ); 
        plt.savefig("results/" + str(i) + "_stor_op.png", dpi=300, bbox_inches='tight')

    return


def co2_op():

    results_file = 'results.db'
    database = 'model/' + results_file
    dat = sqlite3.connect(database)
    scalars = pd.read_sql_query("SELECT * FROM scalars", dat).set_index('name')
    mco2_flue = pd.read_sql_query("SELECT * FROM mco2_flue", dat).set_index('t', drop=True).level[:-1]
    mco2_flue.index = mco2_flue.index.astype('int64')
    mco2_capt = pd.read_sql_query("SELECT * FROM mco2_capt", dat).set_index('t', drop=True).level[:-1]
    mco2_capt.index = mco2_capt.index.astype('int64')
    mco2_em = pd.read_sql_query("SELECT * FROM mco2_em", dat).set_index('t', drop=True).level[:-1]
    mco2_em.index = mco2_em.index.astype('int64')
    ra_co2 = pd.read_sql_query("SELECT * FROM ra_co2", dat).set_index('t', drop=True).level[:-1]
    ra_co2.index = ra_co2.index.astype('int64')
    rd_co2 = pd.read_sql_query("SELECT * FROM rd_co2", dat).set_index('t', drop=True).level[:-1]
    rd_co2.index = rd_co2.index.astype('int64')
    yop_co2 = pd.DataFrame({'mflue': mco2_flue, 'mcapt': mco2_capt, 'mem': mco2_em})

    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 25})
    p_index = np.array(yop_co2.index).copy()
    width = 0.8

    plt.bar(p_index, 
            yop_co2.loc[:,'mem'], 
            width, 
            label='CO$_2$ emissions',
            color='tab:olive')

    plt.bar(p_index, 
            yop_co2.loc[:,'mcapt'], 
            width,
            bottom = yop_co2.loc[:,'mem'], 
            label='CO$_2$ captured',
            color='tab:red')

    plt.plot(p_index, 
             yop_co2.loc[:,'mflue'], 
             label='CO$_2$ in flue gas', 
             linewidth=3, 
             color='xkcd:dark green')

    plt.xlabel('Time (hrs)')
    plt.legend(ncol=1, prop={'size':25})
    plt.ylabel('CO$_2$ flowrate (ton/hr)')
    plt.ylim([0, (yop_co2.loc[:,'mcapt'] + yop_co2.loc[:,'mem']).max() + 50])
    plt.margins(x=0)
    plt.xticks( np.arange(0, len(np.array(yop_co2.index)), math.ceil(len(np.array(yop_co2.index)) / 24) ), np.arange(0, 24, 1) ); 
    plt.savefig("results/" + "co2_op.png", dpi=300, bbox_inches='tight')

    #plot rates of co2 abs/des

    # plt.figure(figsize=(15, 10))
    # plt.rcParams.update({'font.size': 25})
    # p_index = np.array(yop_co2.index).copy()
    # ra_co2.plot(label='ra')
    # # rd_co2.plot(label='rd')
    # plt.legend()
    # plt.savefig("results/" + "rard_co2.png", dpi=300, bbox_inches='tight')

    return


def combined_storop(tech_arr):

    results_file = 'results.db'
    database = 'model/' + results_file
    dat = sqlite3.connect(database)

    train = pd.read_sql_query("SELECT * FROM train", dat).set_index('i', drop=True).level
    P_stor_all = pd.DataFrame(pd.read_sql_query("SELECT * FROM P", dat)[['i','t','level']])
    for tech_ind in train.index:
        for k in P_stor_all.index:
            if P_stor_all.loc[k,'i'] == tech_ind:
                P_stor_all.loc[k,'level'] = P_stor_all.loc[k,'level']*train.loc[tech_ind]

    P_stor_all.set_index(['i','t'], drop=True, inplace=True)
    Pind_stor = pd.DataFrame(index=P_stor_all.loc['ces', :][:-1].index.astype('int64'), columns=tech_arr)

    for i in tech_arr:
        temp_arr = P_stor_all.loc[i, :].level
        temp_arr.index = temp_arr.index.astype('int64')
        Pind_stor[i] = temp_arr

    Pind_stor.columns = ['CES', 'HTTS', 'CAES', 'H$_2$', 'PHS', 'PCM', 'NAS', 'VRFB', 'LIB']

    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 25})
    p_index = np.array(Pind_stor.index).copy()
    width = 1

    cumval=0
    for j in range(len(Pind_stor.columns)):
        plt.bar(p_index, 
        Pind_stor.loc[:, Pind_stor.columns[j]], 
        width,
        bottom = cumval, 
        label= Pind_stor.columns[j])
        cumval = cumval + Pind_stor[Pind_stor.columns[j]]

    plt.legend(ncol=4, prop={'size':20})
    plt.margins(x=0)
    plt.xticks( np.arange(0, len(np.array(Pind_stor.index)), math.ceil(len(np.array(Pind_stor.index)) / 24) ), np.arange(0, 24, 1) ); 
    plt.xlabel('Time (hrs)')
    plt.ylabel('Power flow to/from storage (MW)')
    plt.ylim([cumval.min() - 50, cumval.max() + 50]);
    plt.savefig("results/" + "combined_stor_op.png", dpi=300, bbox_inches='tight')

    return


