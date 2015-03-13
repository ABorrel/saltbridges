'''
Created on 10 mars 2015

@author: borrel
'''
from os import listdir
from re import search



def ScriptConvertICBtoPDB (pr_init, p_filout):

    filout = open (p_filout, "w")
    l_file = listdir(pr_init)

    for f in l_file : 
        if search("\.icb", f):
            filout.write ("openFile \"" + str (pr_init + f) +"\" 0 yes no no no \" append\"\n")
            filout.write ("write pdb " + f[0:-4] + " \""  +  pr_init + f[0:-4] + ".pdb\"\n")
    filout.close ()
    
    return p_filout
            
            

#icm/def> openFile "/home/borrel/saltBridgesProject/GPCRDock2010/PDB_conserved/D3_0400_0001.icb" 0 yes no no no " append" 
#icm/D3_A_1> write pdb a_D3_A_1. "/home/borrel/icm-browser-pro-3.8-3/D3_A_1.pdb"


        
        
        