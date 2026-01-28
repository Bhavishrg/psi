import os
import sys
import re

# directory = 'Results_Emgraph_Init/2_PC/10000_Nodes/TestRun'
directory = sys.argv[1]


aggpath = os.path.join(directory, "agregate_stat.log")
agg = open(aggpath, "w")
agg.truncate(0) 


max_p = ""
max_time = 0

preproc_time = 0
online_time = 0

online_comm = 0
preproc_comm = 0

for filename in os.listdir(directory):
    if "party_0" in filename:
        
        filepath = os.path.join(directory, filename)
        f = open(filepath, "r")
        for x in f:
            if "preproc time:" in x:
                time = re.findall(r'\d+\.\d+', x)
                agg.write(filename + " preproc_time: " + time[0] + " ms" + "\n")
                preproc_time = float(time[0])
                
            if "preproc sent:" in x:
                comm = re.findall(r'\d+', x)
                preproc_comm = float(comm[0])
                agg.write(filename + " preproc_comm: " + comm[0] + " bits\n")
        f.close()
    
    else:
        
        filepath = os.path.join(directory, filename)
        f = open(filepath, "r")
        for x in f:
            if "online time:" in x:
                time = re.findall(r'\d+\.\d+', x)
                agg.write(filename + " " + time[0] + " ms" + "\n")
                t = float(time[0])
                
                if max_time < t:
                    max_time =  t
                    max_p = filename
                    
            if "online sent:" in x:
                comm = re.findall(r'\d+', x)
                online_comm += float(comm[0])
                agg.write(filename + " online_comm: " + comm[0] + " bits" "\n")
                
        f.close()

agg.write("\nOnline time: " + str(max_time) + "ms\n")
agg.write("Online comm: " + str(online_comm) + "bits\n")
agg.write("Max time party: " + max_p + "\n")

agg.write("\nPreproc time: " + str(preproc_time) + "ms\n")
agg.write("Preproc comm: " + str(preproc_comm) + "bits\n")
agg.close()






