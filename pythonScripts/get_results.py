import os
import sys
import re

# directory = 'Results_Emgraph_Init/2_PC/10000_Nodes/TestRun'
directory = sys.argv[1]


aggpath = os.path.join(directory, "agregate_stat.log")
agg = open(aggpath, "w")
agg.truncate(0) 

max_time = 0
max_p = ""

for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    f = open(filepath, "r")
    for x in f:
        if "Online Eval" in x:
            print(x)
        if "online test time:" in x:
            time = re.findall(r'\d+\.\d+', x)
            agg.write(filename + " " + time[0] + "\n")
            t = float(time[0])
            
            if max_time < t:
                max_time =  t
                max_p = filename
    f.close()

agg.write("\nMax time: " + str(max_time) + "\n")
agg.write("\nMax time party: " + max_p + "\n")
agg.close()






