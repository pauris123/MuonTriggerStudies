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

names=pn.read_csv("Step3_FileNames.txt", sep='\t',names=["nos"]) 
print("I have the libraries and names")

# File counting stuff
notik=0
tik=1 #len(names["nos"])

for aiziet in range(notik,tik):
    print("\n {0}/{1}, {2}".format(aiziet+1,tik,wall_time(time.time()-bigBang)))

    start=time.time()       
    rootfile="/eos/cms/store/group/phys_bphys/chic_hlt/"+names["nos"][aiziet]
    events = Events("root://eoscms/"+rootfile)
    #print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))
    
    muon_charge = [] # Lists for read in L1T muons
    muon_pt = []
    muon_eta = []
    muon_phi = []
    muon_iev = []
    muon_i = []
    
    muon_pdgId = [] # Lists for decayed muons from Chi->J/Psi->MuMu in GenP  (1st task)!
    muon_pt_gen = []
    muon_eta_gen = []
    muon_phi_gen = []
    muon_iev_gen = []
    muon_i_gen = []
    
    di_muon_dR_genP = [] # dR between muons that come from Chi->J/Psi-> Mu+Mu decay in GenParticles part (2nd task)!
    di_muon_dR_genP_iev = [] # List to help identify the Chi (event) that this dR is calculated for
    di_muon_dR_genP_i = [] # List to help identify the Chi (i) that this dR is calculated for
    
    
    for iev,event in enumerate(events): # Reading in events from the file
        event.getByLabel(l1MuonLabel, l1Muons)
        event.getByLabel(genParticlesLabel, genParticles)

        # l1Muons
        for i,l1Mu in enumerate(l1Muons.product()): # Checking all L1T muons in single event
            
            muon_charge.append(l1Mu.charge())
            muon_pt.append(l1Mu.pt())
            muon_eta.append(l1Mu.eta())
            muon_phi.append(l1Mu.phi()) 
            muon_iev.append(iev)
            muon_i.append(i)   
    
         #genParticles      
        for i,genP in enumerate(genParticles.product()):
            if (genP.pdgId() == 20443 and genP.numberOfDaughters() == 2) : # Takes Chi that has 2 daughters
                
                if genP.daughter(0).numberOfDaughters() > 0: # Cheks the grand-daughters
                    
                    muon_eta_for_genP_dR_0 = [] # Lists for calculating dR between both muons that come from J/Psi decay in GenParticles
                    muon_phi_for_genP_dR_0 = [] # Lists for calculating dR between both muons that come from J/Psi decay in GenParticles
                    
                    for y in range(genP.daughter(0).numberOfDaughters()):
                        # Prints J/Psi or photon decay products. In this case muons or same strange photons

                        if (genP.daughter(0).daughter(y).pdgId() == 13 or genP.daughter(0).daughter(y).pdgId() == -13): # Cheks for decay muons
                            
                            #Takes the decay Mu parameters and adds them to the list
                            muon_pdgId.append(genP.daughter(0).daughter(y).pdgId())
                            muon_pt_gen.append(genP.daughter(0).daughter(y).pt())
                            muon_eta_gen.append(genP.daughter(0).daughter(y).eta())
                            muon_phi_gen.append(genP.daughter(0).daughter(y).phi())
                            muon_iev_gen.append(iev)
                            muon_i_gen.append(i)
                            
                            #List For calclating dR between 2 muons coming from J/Psi
                            muon_eta_for_genP_dR_0.append(genP.daughter(0).daughter(y).eta())
                            muon_phi_for_genP_dR_0.append(genP.daughter(0).daughter(y).phi())
                    
                    if len(muon_eta_for_genP_dR_0) == 2: # Check whether there is Di-muon pair from J/Psi to calculate dR for, from GenP
                        di_muon_dR_genP.append(((muon_eta_for_genP_dR_0[0]-muon_eta_for_genP_dR_0[1])**2+(muon_phi_for_genP_dR_0[0]-muon_phi_for_genP_dR_0[1])**2)**0.5)
                        di_muon_dR_genP_iev.append(iev)
                        di_muon_dR_genP_i.append(i)
            
                
                if genP.daughter(1).numberOfDaughters() > 0: # Checks all the same things for the 2nd daughter and following grand-daughters
                    
                    muon_eta_for_genP_dR_1 = [] # Lists for calculating dR between both muons that come from J/Psi decay in GenParticles
                    muon_phi_for_genP_dR_1 = [] # Lists for calculating dR between both muons that come from J/Psi decay in GenParticles
                    
                    for y in range(genP.daughter(1).numberOfDaughters()):
                        
                        if (genP.daughter(1).daughter(y).pdgId() == 13 or genP.daughter(1).daughter(y).pdgId() == -13):
                            
                            muon_pdgId.append(genP.daughter(1).daughter(y).pdgId())
                            muon_pt_gen.append(genP.daughter(1).daughter(y).pt())
                            muon_eta_gen.append(genP.daughter(1).daughter(y).eta())
                            muon_phi_gen.append(genP.daughter(1).daughter(y).phi())
                            muon_iev_gen.append(iev)
                            muon_i_gen.append(i)
                            
                            muon_eta_for_genP_dR_1.append(genP.daughter(1).daughter(y).eta())
                            muon_phi_for_genP_dR_1.append(genP.daughter(1).daughter(y).phi())
                            
                    if len(muon_eta_for_genP_dR_1) == 2:
                        di_muon_dR_genP.append(((muon_eta_for_genP_dR_1[0]-muon_eta_for_genP_dR_1[1])**2+(muon_phi_for_genP_dR_1[0]-muon_phi_for_genP_dR_1[1])**2)**0.5)
                        di_muon_dR_genP_iev.append(iev)
                        di_muon_dR_genP_i.append(i)
        
                
    
    muon_dR_min = [] # Matching genP and L1T muons by dR, the list of minimal/matching dR
    muon_matched_L1T_number = [] # Matched L1T muon parameters to genP muons. 1st L1T matched muon in this list corresponds to the first muon in genP list 
    # So the list might go - (10th L1T muon -> 1st GenP muon), (24th L1T muon -> 2nd GenP muon) and so on, so have to follow the iev and i numbers of L1T muon
    
    muon_matched_L1T_charge = []  # Lists of matched L1T muon parameters, so it is easier to calculate stuff onwards
    muon_matched_L1T_pt = []
    muon_matched_L1T_eta = []
    muon_matched_L1T_phi = []
    muon_matched_L1T_iev = []
    muon_matched_L1T_i = []
    
    
    #Calculating the dR values for each genP muon matched with every L1T muon
    for j in range(len(muon_pdgId)): # Take the genP muon
    
        muon_dR_min_calc = []
    
        for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
            if (muon_charge[k] == -1 and muon_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_pdgId[j] == -13):  
                
                muon_dR_min_calc.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5) #Calculate dR, add it to the list
        
            else:
                continue
        
        muon_dR_min.append(min(muon_dR_min_calc))  # Take the min value of dR and add it to the list for matched muon dR's
        muon_matched_L1T_number.append(muon_dR_min_calc.index(min(muon_dR_min_calc)))  # Find the index of the L1T muon, that gave the lowest dR value for the particular genP muon, and add it to the extra list.
    
    for n in range(len(muon_matched_L1T_number)): # Creates new L1T muon parameter lists, based of matching. So first element in these list is matched L1T muon paramaters, that got matched to the 1st muon in genP list. 2nd elemnt here is matched L1T muon for the 2nd genP. Every two lines are a pair. 1st and 2nd is from the 1st Chi matched, 3rd and 4th is from the 2nd Chi matched and so on.
        
        muon_matched_L1T_charge.append(muon_charge[muon_matched_L1T_number[n]])    
        muon_matched_L1T_pt.append(muon_pt[muon_matched_L1T_number[n]])
        muon_matched_L1T_eta.append(muon_eta[muon_matched_L1T_number[n]])
        muon_matched_L1T_phi.append(muon_phi[muon_matched_L1T_number[n]])
        muon_matched_L1T_iev.append(muon_iev[muon_matched_L1T_number[n]])
        muon_matched_L1T_i.append(muon_i[muon_matched_L1T_number[n]])
    
    if len(muon_pdgId) == len(muon_matched_L1T_charge): # Checking whether all genP muons got matched with something
        print("All gucci with the matching")

    
    # Here Im writing info to files, but havent started that.
            
    #fo=open("Muon_dR_inFile_decay_muon.txt","a")
    #for z in range(len(muon_dR)):
     #   fo.write(str(muon_dR[z])+"\t"+str(aiziet)+"\n")
    #fo.close()
    
    #fom=open("Muon_dR_min_decay_muon.txt","a")
    #for v in range(len(muon_dR_min)):
     #   fom.write(str(muon_dR_min[v])+"\t"+str(aiziet)+"\n")
    #fom.close()
    
    ##fon=open("Muon_dR_from_JPsi_genP.txt","a")
    ##for b in range(len(di_muon_dR_genP)):
    ##    fon.write(str(di_muon_dR_genP[b])+"\t"+str(di_muon_dR_genP_iev[b])+"\t"+str(di_muon_dR_genP_i[b])+"\n")
    ##fon.close()                
                
    #print("List of cross-checked muon dR is this long -> "+str(len(muon_dR)))
    #print("List of cross-checked muon best min dR is this long -> "+str(len(muon_dR_min)))             
    #print("GenParticles muons found -> "+str(len(muon_pdgId))+" L1T muons found -> "+str(len(muon_charge)))                           
                
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
                              
                               
