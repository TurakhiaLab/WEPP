import os
import sys
import subprocess
import math
import random

cmd = "matUtils extract -i global_assignments.pb -m {} -v temp".format(str(sys.argv[1]))
os.system(cmd)
os.system("grep -v \"#\" temp")

cmd = "grep -v \"#\" temp | wc -l"
count = subprocess.check_output(cmd, shell=True)
frac = int(str(sys.argv[2])[2])
lim = math.ceil(float(frac) * float(count) / 10)

for i in range(1, int(lim)+1):
    print(random.randint(1,int(count)))