import os
import CPUUsage
from time import sleep

naccess = "/home/borrel/softwares/naccess/naccess"



def runNACESS (path_file_in, probe = "", verbose = 1):
    """
    Run NACCESS with file pdb 
    args: -> filin pdb
    return: path files .asa and .rsa
    """
    
    # check if file exist
    if os.path.isfile(path_file_in[0:-4] + ".asa") and os.path.isfile(path_file_in[0:-4] + ".rsa") : 
        print "continue"
        return [path_file_in[0:-4] + ".asa", path_file_in[0:-4] + ".rsa"]
    
#     
    time_cpu =  CPUUsage.CPUsage().compute()
    while float(time_cpu) > 70 :
        print "pause ", CPUUsage.CPUsage().compute()
        sleep (0.5)
        time_cpu =  CPUUsage.CPUsage().compute()
        continue
    
    cmd_run = naccess + " " + path_file_in + "&"
    os.system (cmd_run)
    
    return [path_file_in[0:-4] + ".asa", path_file_in[0:-4] + ".rsa"]

