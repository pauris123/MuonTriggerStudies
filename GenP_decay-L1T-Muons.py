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
    #print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))
    
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
        for genP in genParticles.product():
            if (genP.pdgId() == 20443 and genP.numberOfDaughters() == 2) :
                # Prints the last Chi that decays into J/PSi and Photon
                #print("Mother: pt {0:1f}, pdgId {1:2f}, N of D {2}, Daughter 1 pdgId {3}, Daughter 2 pdgId {4}".format( genP.pt(), genP.pdgId(), genP.numberOfDaughters(),genP.daughter(0).pdgId(),genP.daughter(1).pdgId())) 
                # Prints J/Psi or Photon parameters
                #print("Daughter 1: pdgId {0}, pt {1}, eta {2}, phi {3}, N o D {4}".format(genP.daughter(0).pdgId(),genP.daughter(0).pt(),genP.daughter(0).eta(),genP.daughter(0).phi(),genP.daughter(0).numberOfDaughters()))
                
                if genP.daughter(0).numberOfDaughters() > 0:
                    
                    for y in range(genP.daughter(0).numberOfDaughters()):
                        # Prints J/Psi or photon decay products. In this case muons or same strange photons
                        #print("Grand-Daughter 1.{0}: pdgId {1}, pt {2}, eta {3}, phi {4}, N o D {5}".format(y+1,genP.daughter(0).daughter(y).pdgId(),genP.daughter(0).daughter(y).pt(),genP.daughter(0).daughter(y).eta(),genP.daughter(0).daughter(y).phi(),genP.daughter(0).daughter(y).numberOfDaughters()))
                        if (genP.daughter(0).daughter(y).pdgId() == 13 or genP.daughter(0).daughter(y).pdgId() == -13):
                            
                            muon_pdgId.append(genP.daughter(0).daughter(y).pdgId())
                            muon_pt_gen.append(genP.daughter(0).daughter(y).pt())
                            muon_eta_gen.append(genP.daughter(0).daughter(y).eta())
                            muon_phi_gen.append(genP.daughter(0).daughter(y).phi())
                                
                # Prints J/Psi or Photon parameters
                #print("Daughter 2: pdgId {0}, pt {1}, eta {2}, phi {3}, N o D {4}".format(genP.daughter(1).pdgId(),genP.daughter(1).pt(),genP.daughter(1).eta(),genP.daughter(1).phi(),genP.daughter(1).numberOfDaughters()))
                
                if genP.daughter(1).numberOfDaughters() > 0:
                    
                    for y in range(genP.daughter(1).numberOfDaughters()):
                        # Prints J/Psi or photon decay products. In this case muons or same strange photons
                        #print("Grand-Daughter 2.{0}: pdgId {1}, pt {2}, eta {3}, phi {4}, N o D {5}".format(y+1,genP.daughter(1).daughter(y).pdgId(),genP.daughter(1).daughter(y).pt(),genP.daughter(1).daughter(y).eta(),genP.daughter(1).daughter(y).phi(),genP.daughter(1).daughter(y).numberOfDaughters())) 
                        if (genP.daughter(1).daughter(y).pdgId() == 13 or genP.daughter(1).daughter(y).pdgId() == -13):
                            
                            muon_pdgId.append(genP.daughter(1).daughter(y).pdgId())
                            muon_pt_gen.append(genP.daughter(1).daughter(y).pt())
                            muon_eta_gen.append(genP.daughter(1).daughter(y).eta())
                            muon_phi_gen.append(genP.daughter(1).daughter(y).phi())
        
                
    muon_dR = []
    muon_dR_min = []
    
    for j in range(len(muon_pdgId)):
    
        muon_dR_min_calc = []
    
        for k in range(len(muon_charge)):
            if (muon_charge[k] == -1 and muon_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_pdgId[j] == -13):  
                
                muon_dR.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5)
                muon_dR_min_calc.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5)
        
            else:
                continue
        
        muon_dR_min.append(min(muon_dR_min_calc))
                
    fo=open("Muon_dR_inFile_decay_muon.txt","a")
    for z in range(len(muon_dR)):
        fo.write(str(muon_dR[z])+"\t"+str(aiziet)+"\n")
    fo.close()
    
    fom=open("Muon_dR_min_decay_muon.txt","a")
    for v in range(len(muon_dR_min)):
        fom.write(str(muon_dR_min[v])+"\t"+str(aiziet)+"\n")
    fom.close()                
                
    print("List of cross-checked muon dR is this long -> "+str(len(muon_dR)))
    print("List of cross-checked muon best min dR is this long -> "+str(len(muon_dR_min)))             
    print("GenParticles muons found -> "+str(len(muon_pdgId))+" L1T muons found -> "+str(len(muon_charge)))                           
                
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
                              
                               
