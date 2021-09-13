import numpy as np
import pandas as pn
import matplotlib.pyplot as plt

bi=100
file1="Muon_invM_DR_inFile.txt"                   

df=pn.read_csv(file1, sep='\t',names=["mass","invM_over_dR"]) # Data file

#plt.hist(df["mass"],bins=bi,range=[0,10], alpha=0.5, histtype=u'step', label='Di-muon inv.mass', color='b')
plt.hist(df["invM_over_dR"],bins=bi,range=[0,20], alpha=0.5, histtype=u'step', label='Di-muon inv.mass over dR', color='r')
plt.legend()
#plt.xlabel("Invariant mass, GeV")
plt.xlabel("Invariant mass over dR, GeV")
plt.ylabel("Frequency")
plt.title("Muon test MC")
#plt.xticks(np.arange(0,11,1))
plt.savefig("Result_Plot.png")
plt.show()
