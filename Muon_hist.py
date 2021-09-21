import numpy as np
import pandas as pn
import matplotlib.pyplot as plt

bi=100
file1="Muon_dR_inFile_All_muons.txt"
file2="Muon_dR_inFile_decay_muon.txt"
file3="Muon_dR_min_decay_muon.txt"                 

df=pn.read_csv(file1, sep='\t',names=["dR","fileNumber"]) # Data file
dg=pn.read_csv(file2, sep='\t',names=["dR","fileNumber"]) # Data file
dh=pn.read_csv(file3, sep='\t',names=["dR","fileNumber"]) # Data file

plt.hist(df["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='All muons cross-checked dR', color='b')
plt.hist(dg["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Decay muons cross-checked dR', color='r')
plt.hist(dh["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Each decay muon's min dR value", color='g')
plt.legend()
#plt.xlabel("Invariant mass, GeV")
plt.xlabel("dR, GeV")
plt.ylabel("Frequency")
plt.title("Muon private MC files Chi")
plt.xticks(np.arange(0,11,1))
plt.savefig("Result_Plot_Final.png")
plt.show()
