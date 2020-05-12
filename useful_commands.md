#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 13:31:55 2020

@author: transcend
"""

OVERLAPS
#Generate a dataframe from a list of Class objects
ovlpDF = pd.DataFrame([vars(s) for s in scan.ovlpClusterList])

#Make a DataFrame of lists into a one-dimensional list
flat_df = []
for sublist in ovlpDF['spacing']:
    for item in sublist:
        flat_df.append(item)
        
#Turn the flat list into a flat DF
flat_df = pd.DataFrame(flat_df)

#Count the occurences of values in a DataFrame
flat_df_count = pd.DataFrame(flat_df[0].value_counts())

#Sort the DF by value
order = [-8, -7, -6, -5, -4, -3, -2, -1] """(or np.arange(-8,0)).tolist()))"""
flat_df_count = flat_df_count.loc[order]

#Generating a bar plot from a flat dataframe
fig, ax = plt.subplots()
flat_df_count[0].plot(ax=ax, kind='bar')

-----------

#UNPROCESSED

unpDF = pd.DataFrame([vars(s) for s in scan.unpClusterList]) #Generate a dataframe from a list of Class objects

flat_df_unp = [] #Make a DataFrame of lists into a one-dimensional list
for sublist in unpDF['spacing']:
    for item in sublist:
        flat_df_unp.append(item)
        
flat_df_unp = pd.DataFrame(flat_df_unp) #Turn the flat list into a flat DF

flat_df_count_unp = pd.DataFrame(flat_df_unp[0].value_counts()) #Count the occurences of values in a DataFrame

order = np.arange(-8,36).tolist() #this needs to be changed, to account for loss of overlaps (-8 first, just add lists) 
flat_df_count_unp = flat_df_count_unp.loc[order] #Sort the DF by value

import matplotlib.pyplot as plt #Generating a bar plot from a flat dataframe
fig, ax = plt.subplots()
flat_df_count_unp[0].plot(ax=ax, kind='bar')