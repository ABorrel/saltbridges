import os
import CPUUsage
from time import sleep
import psutil

naccess = "/home/borrel/softwares/naccess/naccess"



def runNACESS (path_file_in, pr_result, probe = "", verbose = 1, multi_run = 0):
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
    if multi_run == 1 : 
        

        #CPU  
        time_cpu =  CPUUsage.CPUsage().compute()
        
        #MEM
        proc = psutil.phymem_usage()
    
        memory = proc.percent
        print memory
    
        while float(time_cpu) > 70 or float(memory) > 30 :
            print "pause ", CPUUsage.CPUsage().compute(), memory
            sleep (0.5)
            time_cpu =  CPUUsage.CPUsage().compute()
            memory = psutil.phymem_usage().percent
            continue
        
        cmd_run = naccess + " " + path_file_in + "&"
        print cmd_run
        os.system (cmd_run)
        
        try : 
            os.system("mv *.asa " +  pr_result)
            os.system("mv *.rsa " +  pr_result)
            os.system("rm *.log")
        except : pass
        
        return [path_file_in[0:-4] + ".asa", path_file_in[0:-4] + ".rsa"]
    
    
    else : 
        
        cmd_run = naccess + " " + path_file_in
        print cmd_run
        os.system (cmd_run)
        try : 
            os.system("mv *.asa " +  pr_result)
            os.system("mv *.rsa " +  pr_result)
            os.system("rm *.log")
        except : pass
        
        return [path_file_in[0:-4] + ".asa", path_file_in[0:-4] + ".rsa"]
    
    
    return "ERROR"
    
      
        
        
