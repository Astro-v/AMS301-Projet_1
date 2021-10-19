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
plt.colorbar()
plt.savefig("seqJ.png", dpi=150)

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
plt.colorbar()
plt.savefig("parJ.png", dpi=150)

fil = open("gaussSeidelSeq.txt","r")
line = fil.readlines()
fil.close()
line = [line[i][:-1] for i in range(0,len(line))]
for i in range(0,len(line)):
    line[i] = line[i].split()
    for j in range(0,len(line[i])):
        line[i][j]=float(line[i][j])

plt.figure()
plt.imshow(line)
plt.colorbar()
plt.savefig("seqGS.png", dpi=150)


fil = open("gaussSeidelPara.txt","r")
line = fil.readlines()
fil.close()
line = [line[i][:-1] for i in range(0,len(line))]
for i in range(0,len(line)):
    line[i] = line[i].split()
    for j in range(0,len(line[i])):
        line[i][j]=float(line[i][j])

plt.figure()
plt.imshow(line)
plt.colorbar()
plt.savefig("parGS.png", dpi=150)

P = [1,2,5,10,20,40,80]
L100 = [3.49605,1.70169,0.710672,0.383819,0.215062,0.182067,0.190505]
L100e = [3.49605/p for p in P]
L1000 = [406.828,204.353,86.8638,35.2449,17.9862,10.8216,5.85413]
L1000e = [406.828/p for p in P]
plt.figure()
plt.loglog(P,L100,label="Mesure")
plt.loglog(P,L100e,label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Temps de calcul")
plt.legend()
plt.savefig("scalForteJ100.png", dpi=150)
plt.figure()
plt.loglog(P,L1000,label="Mesure")
plt.loglog(P,L1000e,label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Temps de calcul")
plt.legend()
plt.savefig("scalForteJ1000.png", dpi=150)

P = [1,2,5,10,20,40,80]
S1 = [100*np.sqrt(p) for p in P]
# [100, 141, 224, 316, 447, 632, 894]
# [1000, 1414, 2236, 3162, 4472, 6325, 8944]
# L1 = [3.01389,3.13379,3.33244,3.44574,3.78009,4.45922, 5.03782] # proportionnalite selon x et y
L1 = [4.77241,4.85676,5.19079,5.39416,5.72231,7.25802]
S1 = [100*p for p in P]
L2 = [3.01389,3.21376,3.29145,3.38424,3.60109,4.41684,4.82369] # selon x seulement
S = [L1[0]/L1[i] for i in range(0,len(L1))]
plt.figure()
plt.loglog(P,S,label="Mesure")
plt.loglog(P,[1 for i in range(len(S))],label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Efficacite")
plt.legend()
plt.savefig("scalFaibleJ.png", dpi=150)

P = [1,2,5,10,20,40,80]
L100 = [3.32907,1.68622,0.729999,0.380173,0.230128,0.181992,0.195936]
L100e = [3.32907/p for p in P]
L1000 = [318.941,159.278,68.2452,34.9886,18.8786,12.0077,7.18206]
L1000e = [318.941/p for p in P]
plt.figure()
plt.loglog(P,L100,label="Mesure")
plt.loglog(P,L100e,label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Temps de calcul")
plt.legend()
plt.savefig("scalForteG100.png", dpi=150)
plt.figure()
plt.loglog(P,L1000,label="Mesure")
plt.loglog(P,L1000e,label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Temps de calcul")
plt.legend()
plt.savefig("scalForteG1000.png", dpi=150)

P = [1,2,5,10,20,40,80]
S1 = [1000*np.sqrt(p) for p in P]
# [100, 141, 224, 316, 447, 632, 894]
# [1000, 1414, 2236, 3162, 4472, 6325, 8944]
L1 = [3.32907,3.38073,3.42005,3.69482,4.09652,5.01967, 5.53075] # proportionnalite selon x et y
S1 = [100*p for p in P]
L2 = [3.01389,3.21376,3.29145,3.38424,3.60109,4.41684,4.82369] # selon x seulement
S = [L1[0]/L1[i] for i in range(0,len(L1))]
plt.figure()
plt.loglog(P,S,label="Mesure")
plt.loglog(P,[1 for i in range(len(S))],label="Ideal")
plt.xlabel("Nombre de processus")
plt.ylabel("Efficacite")
plt.legend()
plt.savefig("scalFaibleG.png", dpi=150)