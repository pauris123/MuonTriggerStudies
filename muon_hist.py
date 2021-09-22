import numpy as np
import pandas as pn
import matplotlib.pyplot as plt

bi=100
file1="Muon_dR_from_JPsi_genP.txt"
file2="Muon_dR_from_Matching_genP_L1T.txt"
file3="Di_Muon_dR_invM_invM_o_dR.txt"                 

df=pn.read_csv(file1, sep='\t',names=["dR","event","i","file"]) # Data file for di-muon dR from J-Psi in genParticles
dg=pn.read_csv(file2, sep='\t',names=["dR","event","igen","iL1T","file"]) # Data file for muon matching dR from genParticles to L1T muons
dh=pn.read_csv(file3, sep='\t',names=["dR","invM","invM_dR","file"]) # Data file for di-muon dR,inv.M,inv.M/dR from J/Psi in L1T

fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
n, bins, patches = ax1.hist(df["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Di-muon dR from J/Psi in genParticles', color='b')
#ax1.hist(df["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Di-muon dR from J/Psi in genParticles', color='b')
ax1.set_xlabel("dR")
ax1.set_ylabel("Frequency")
#ax1.xticks(np.arange(0,11,1))
#ax1.savefig("dR_JPsiMuons_genP.png")
ax1.set_title("Di-muon dR from J/Psi in genParticles") # Nezinu vai der, maybe setlabel

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
n, bins, patches = ax2.hist(dg["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Muon matching dR from genParticles to L1T', color='r')
#ax2.hist(dg["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Muon matching dR from genParticles to L1T', color='r')
ax2.set_xlabel("dR")
ax2.set_ylabel("Frequency")
#ax2.xticks(np.arange(0,11,1))
#ax2.savefig("Muon_matching_dR")
ax2.set_title("Muon matching dR from genParticles to L1T") # Nezinu vai der, maybe setlabel

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
n, bins, patches = ax3.hist(dh["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon dR from J/Psi in L1T", color='g')
#ax3.hist(dh["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon dR from J/Psi in L1T", color='g')
ax3.set_xlabel("dR")
ax3.set_ylabel("Frequency")
#ax3.xticks(np.arange(0,11,1))
#ax3.savefig("dR_JPsiMuons_L1T.png")
ax3.set_title("Di-muon dR from J/Psi in L1T") # Nezinu vai der, maybe setlabel

fig4 = plt.figure()
ax4 = fig4.add_subplot(1,1,1)
n, bins, patches = ax4.hist(dh["invM"],bins=bi,range=[0,15], alpha=0.8, histtype=u'step', label="Di-muon inv.M from J/Psi in L1T", color='m')
#ax4.hist(dh["invM"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon inv.M from J/Psi in L1T", color='m')
ax4.set_xlabel("invM, GeV")
ax4.set_ylabel("Frequency")
#ax4.xticks(np.arange(0,11,1))
#ax4.savefig("InvM_JPsiMuons_L1T.png")
ax4.set_title("Di-muon inv.M from J/Psi in L1T") # Nezinu vai der, maybe setlabel

fig5 = plt.figure()
ax5 = fig5.add_subplot(1,1,1)
n, bins, patches = ax5.hist(dh["invM_dR"],bins=bi,range=[0,15], alpha=0.8, histtype=u'step', label="Di-muon inv.M/dR from J/Psi in L1T", color='k')
#ax5.hist(dh["inM/dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon inv.M/dR from J/Psi in L1T", color='k')
ax5.set_xlabel("invM/dR")
ax5.set_ylabel("Frequency")
#ax5.xticks(np.arange(0,11,1))
#ax5.savefig("InvM_o_dR_JPsiMuons_L1T.png")
ax5.set_title("Di-muon inv.M/dR from J/Psi in L1T") # Nezinu vai der, maybe setlabel

#plt.hist(df["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Di-muon dR from J/Psi in genParticles', color='b')
#plt.hist(dg["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label='Muon matching dR from genParticles to L1T', color='r')
#plt.hist(dh["dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon dR from J/Psi in L1T", color='g')
#plt.hist(dh["invM"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon inv.M from J/Psi in L1T", color='m')
#plt.hist(dh["inM/dR"],bins=bi,range=[0,10], alpha=0.8, histtype=u'step', label="Di-muon inv.M/dR from J/Psi in L1T", color='k')
#plt.legend()
#plt.xlabel("Inv.M, GeV")
#plt.xlabel("Inv.M/dR, GeV")
#plt.xlabel("dR, GeV")
#plt.ylabel("Frequency")
#plt.title("Muon private MC files Chi")
#plt.xticks(np.arange(0,11,1))
#plt.savefig("Result_Plot_Final.png")
plt.show()


