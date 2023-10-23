# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 16:42:27 2023

This is a suggested method for plotting the output results of multiscat.

@author: SamLambrick
"""

import numpy as np
import pandas as pd
import seaborn as sns

# Default theme
sns.set_theme()

def import_multiscat(fname):    
    """Impot standrad multiscat output into a pandas data frame."""
    d = pd.read_csv(fname, skiprows=7, delim_whitespace=True, 
                    header=None, names=['#','n1','n2','I'])
    d.drop(columns=['#'], inplace=True)
    return(d)

d = import_multiscat('diffrac10001.out')

d2 = d.pivot(index='n1', columns='n2', values='I')
ax = sns.heatmap(d2, cmap='viridis', cbar_kws={'label' : '$P(n_1,n_2)$'})
ax.set_aspect('equal')
ax.set_xlabel('$n_1$')
ax.set_ylabel('$n_2$')
