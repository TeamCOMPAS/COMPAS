# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

fname = 'COMPAS_Output/COMPAS_Output.h5'
myf = h5.File(fname, 'r')
myf.keys()

SPs = myf['BSE_System_Parameters']
SNe = myf['BSE_Supernovae']
SNe.keys()

#print(SPs['Stellar_Type(1)'][()])
#print(SPs['Stellar_Type@ZAMS(1)'][()])
printCompasDetails(SNe)

SLs = myf['BSE_Switch_Log']
SLs.keys()

for key in SLs.keys():
    print(SLs[key][()])


