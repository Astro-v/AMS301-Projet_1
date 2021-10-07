import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

fil = open("jacobiSeq.txt","r")
line = fil.readlines()
fil.close()
line = [line[i][:-1] for i in range(0,len(line))]
for i in range(0,len(line)):
    line[i] = line[i].split()
    for j in range(0,len(line[i])):
        line[i][j]=float(line[i][j])

plt.figure()
plt.imshow(line)
plt.savefig("seq.png", dpi=150)

fil = open("jacobiPara.txt","r")
line = fil.readlines()
fil.close()
line = [line[i][:-1] for i in range(0,len(line))]
for i in range(0,len(line)):
    line[i] = line[i].split()
    for j in range(0,len(line[i])):
        line[i][j]=float(line[i][j])

plt.figure()
plt.imshow(line)
plt.savefig("par.png", dpi=150)

P = [1,2,5,10,20,40,80]
L100 = [3.49605,1.70169,0.710672,0.383819,0.215062,0.182067,0.190505]
L100e = [3.49605/p for p in P]
L1000 = [406.828,204.353,86.8638,35.2449,17.9862,10.8216,5.85413]
L1000e = [406.828/p for p in P]
plt.figure()
plt.loglog(P,L100)
plt.loglog(P,L100e)
plt.savefig("scalForte100.png", dpi=150)
plt.figure()
plt.loglog(P,L1000)
plt.loglog(P,L1000e)
plt.savefig("scalForte1000.png", dpi=150)

P = [1,2,5,10,20,40,80]
S1 = [100*np.sqrt(p) for p in P]
# [100, 141, 224, 316, 447, 632, 894]
L1 = [3.01389,3.13379,3.33244,3.44574,3.78009,4.45922, 5.03782] # proportionnalite selon x et y
S1 = [100*p for p in P]
L2 = [3.01389,3.21376,3.29145,3.38424,3.60109,4.41684,4.82369] # selon x seulement
S = [L1[0]/L1[i] for i in range(0,len(L1))]
plt.figure()
plt.loglog(P,S)
plt.loglog(P,[1 for i in range(len(S))])
plt.savefig("scalFaible.png", dpi=150)