# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

import time
bigBang=time.time()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pn

def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s) #lxplus nepatik d
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

l1Muons, l1MuonLabel = Handle("BXVector<l1t::Muon>"), "gmtStage2Digis:Muon:RECO"
genParticles, genParticlesLabel = Handle("vector<reco::GenParticle>"), "genParticles"

seenIt = {} # list of things we've seen (so that we dump them in full only once)
names=pn.read_csv("Step3_FileNames.txt", sep='\t',names=["nos"]) 
print("I have the libraries and names")

notik=0
tik=len(names["nos"])

for aiziet in range(notik,tik):
    print("\n {0}/{1}, {2}".format(aiziet+1,tik,wall_time(time.time()-bigBang)))

    start=time.time()       
    rootfile="/eos/cms/store/group/phys_bphys/chic_hlt/"+names["nos"][aiziet]
    events = Events("root://eoscms/"+rootfile)
    print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))
    
    muon_charge = []
    muon_pt = []
    muon_eta = []
    muon_phi = []
    
    muon_pdgId = []
    muon_pt_gen = []
    muon_eta_gen = []
    muon_phi_gen = []
    
    for iev,event in enumerate(events):
        event.getByLabel(l1MuonLabel, l1Muons)
        event.getByLabel(genParticlesLabel, genParticles)

        #print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())

        # l1Muons
        for i,l1Mu in enumerate(l1Muons.product()):
            #print("i {0}, pt {1:1f}, eta {2:2f}, phi {3:2f}, charge {4:1f}".format(i, l1Mu.pt(), l1Mu.eta(), l1Mu.phi(), l1Mu.charge()))
            muon_charge.append(l1Mu.charge())
            muon_pt.append(l1Mu.pt())
            muon_eta.append(l1Mu.eta())
            muon_phi.append(l1Mu.phi())    
    
         #genParticles
        for i,genP in enumerate(genParticles.product()):
            if (genP.pdgId() == 13 or genP.pdgId() == -13):
                #print("i {0}, pt {1:1f}, eta {2:2f}, phi {3:2f}, pdgId {4:1f}".format(i, genP.pt(), genP.eta(), genP.phi(), genP.pdgId())) 
                muon_pdgId.append(genP.pdgId())
                muon_pt_gen.append(genP.pt())
                muon_eta_gen.append(genP.eta())
                muon_phi_gen.append(genP.phi())
                
    muon_dR = []
    muon_dR_min = []
    
    for k in range(len(muon_charge)):
    
        for j in range(len(muon_pdgId)):
            if (muon_charge[k] == -1 and muon_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_pdgId[j] == -13):  
                
                muon_dR.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5)
        
            else:
                continue
                
    fo=open("Muon_dR_inFile.txt","a")
    for z in range(len(muon_dR)):
        fo.write(str(muon_dR[z])+"\t"+str(aiziet)+"\n")
    fo.close()                
                
    print("List of cross-checked muon dR is this long -> "+str(len(muon_dR)))            
    print("GenParticles muons found -> "+str(len(muon_pdgId))+" L1T muons found -> "+str(len(muon_charge)))                           
                
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
                              
  
