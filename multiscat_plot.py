# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 16:42:27 2023

This is a suggested method for plotting the output results of multiscat.

@author: SamLambrick
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Default theme
sns.set_theme()

def import_multiscat(fname):    
    """Impot standrad multiscat output into a pandas data frame."""
    d = pd.read_csv(fname, skiprows=7, delim_whitespace=True, 
                    header=None, names=['#','n1','n2','I'])
    d.drop(columns=['#'], inplace=True)
    return(d)


d_tmp = import_multiscat('lif_30deg_45deg.out')


def diffraction_plot_square(d):
    d2 = d.pivot_table(index='n1', columns='n2', values='I')
    plt.figure()
    ax = sns.heatmap(d2, cmap='viridis', cbar_kws={'label' : '$P(n_1,n_2)$'})
    ax.set_aspect('equal')
    ax.set_xlabel('$n_2$')
    ax.set_ylabel('$n_1$')
    ax.set_title('$\\theta_i = 30^\\circ$')


def import_multiscat_phi(fname):
    with open(fname) as f:
        lines = f.readlines()
    phis = []
    data = []
    data_set = []
    i = 1
    while i < len(lines):
        if lines[i][1] == 'R':
            phis.append(float(lines[i + 5].split()[6]))
            if i > 1:
                data.append(data_set)
                data_set = []
            i = i + 6
        else:
            data_set.append([float(s) for s in lines[i].split()[1:]])
            if i + 1 == len(lines):
                data.append(data_set)
            i = i + 1
    return(np.array(phis), data)

phis, data = import_multiscat_phi('lif_30deg.out')


d0 = data[np.where(phis == 0)[0][0]]
d0 = pd.DataFrame(np.array(d0), columns=['n1', 'n2', 'I'])

d45 = data[np.where(phis == 45)[0][0]]
d45 = pd.DataFrame(np.array(d45), columns=['n1', 'n2', 'I'])

d175 = data[np.where(phis == 17.5)[0][0]]
d175 = pd.DataFrame(np.array(d175), columns=['n1', 'n2', 'I'])

d275 = data[np.where(phis == 27.5)[0][0]]
d275 = pd.DataFrame(np.array(d275), columns=['n1', 'n2', 'I'])

d325 = data[np.where(phis == 32.5)[0][0]]
d325 = pd.DataFrame(np.array(d325), columns=['n1', 'n2', 'I'])


# 1,0 direction peaks we want
d0t = d0[(d0['n2'] == 0) & (d0['n1'] <= 0) & (d0['n1'] > -4)]
# 1,1 direction peaks we want
d45t = d45[(d45['n2'] == d45['n1']) & (d45['n2'] < 0) & (d45['n2'] > -4)]
# And the other 3 diffraction peaks that we want
d175t = d175[(d175['n2'] == -1) & (d175['n1'] == -3)]
d275t = d275[(d275['n2'] == -1) & (d275['n1'] == -2)]
d325t = d325[(d325['n2'] == -2) & (d325['n1'] == -3)]

d0t_2 = d0t[d0t['n1'] != 0]
d0t_2 = d0t_2.rename(columns = {'n1': 'n2', 'n2': 'n1'})
d175t_2 = d175t.rename(columns = {'n1': 'n2', 'n2': 'n1'})
d275t_2 = d275t.rename(columns = {'n1': 'n2', 'n2': 'n1'})
d325t_2 = d325t.rename(columns = {'n1': 'n2', 'n2': 'n1'})

d = pd.concat([d0t, d45t, d175t, d275t, d325t, d0t_2, d175t_2, d275t_2, d325t_2])
d['n1'] = d['n1'].astype(int)
d['n2'] = d['n2'].astype(int)

d_tmp = d[d['n1'] < 0]
d_tmp = d_tmp.apply(lambda x: -x if x.name == 'n1' else x)
d = pd.concat([d, d_tmp])

d_tmp = d[d['n2'] < 0]
d_tmp = d_tmp.apply(lambda x: -x if x.name == 'n2' else x)
d = pd.concat([d, d_tmp])

diffraction_plot_square(d)
