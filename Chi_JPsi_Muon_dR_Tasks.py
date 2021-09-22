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
tik=len(names["nos"])

error_list = [] # counting errors

for aiziet in range(notik,tik):
    print("\n {0}/{1}, {2}".format(aiziet+1,tik,wall_time(time.time()-bigBang)))

    start=time.time()       
    rootfile="/eos/cms/store/group/phys_bphys/chic_hlt/"+names["nos"][aiziet]
    events = Events("root://eoscms/"+rootfile)
    #print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))
    
    
    di_muon_dR_genP = [] # dR between muons that come from Chi->J/Psi-> Mu+Mu decay in GenParticles part (2nd task)!
    di_muon_dR_genP_iev = [] # List to help identify the Chi (event) that this dR is calculated for
    di_muon_dR_genP_i = [] # List to help identify the Chi (i) that this dR is calculated for
    
    muon_dR_min = [] # Matching genP and L1T muons by dR, the list of minimal/matching dR 3rd task.
    muon_dR_min_iev = []
    muon_dR_min_genP_i = []
    muon_dR_min_L1T_i = []
    
    di_muon_inv_mass = []        # Lists for dR, inv.M and in.M/dR for matched L1T di-muons reconstructed from Chi decay - The end result? 4th task
    di_muon_dR_matched_L1T = []
    inv_M_o_dR_di_muon = []
            
    for iev,event in enumerate(events): # Reading in events from the file
        event.getByLabel(l1MuonLabel, l1Muons)
        event.getByLabel(genParticlesLabel, genParticles)

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
        
        muon_dR_min_matched = []  # Same as muon_dR_min, but for each event, makes writing data more simpler in the next if and for cycles
        muon_matched_L1T_number = []
                # Matched L1T muon parameters to genP muons. 1st L1T matched muon in this list corresponds to the first muon in genP list 
    # So the list might go - (10th L1T muon -> 1st GenP muon), (24th L1T muon -> 2nd GenP muon) and so on, so have to follow the iev and i numbers of L1T muon
    
        muon_matched_L1T_charge = []  # Lists of matched L1T muon parameters, so it is easier to calculate stuff onwards
        muon_matched_L1T_pt = []
        muon_matched_L1T_eta = []
        muon_matched_L1T_phi = []
        muon_matched_L1T_iev = []
        muon_matched_L1T_i = []
        
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
        
                
    

            
    
        #Calculating the dR values for each genP muon matched with every L1T muon
        if (len(muon_pdgId) == 2 or len(muon_pdgId) == 4) and len(muon_charge) >= 2 :
            if len(muon_charge) == 2 and sum(muon_charge) != 0:
                continue
            if len(muon_charge) == 3 and abs(sum(muon_charge)) == 3:
                continue
            if len(muon_charge) == 4 and abs(sum(muon_charge)) == 4:
                continue
            if len(muon_charge) == 5 and abs(sum(muon_charge)) == 5:
                continue
            if len(muon_charge) == 6 and abs(sum(muon_charge)) == 6:
                continue
            if (muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] != muon_charge[1]):
                print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                error_list.append(i)
                continue            
            if (len(muon_charge) == 3 and ((muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] != muon_charge[1]) or (muon_eta[0] == muon_eta[2] and muon_phi[0] == muon_phi[2] and muon_charge[0] != muon_charge[2]) or (muon_eta[2] == muon_eta[1] and muon_phi[2] == muon_phi[1] and muon_charge[2] != muon_charge[1]))):
                print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                error_list.append(i)
                continue
            
            
            for j in range(len(muon_pdgId)): # Take the genP muon
    
                muon_dR_min_calc = []
                muon_dR_min_calc_fix = [] # Fixed problem with matched L1T numbers
    
                for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
                    
                    muon_dR_min_calc_fix.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5) #Temporary fix
                    
                    if (muon_charge[k] == -1 and muon_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_pdgId[j] == -13):  
                
                        muon_dR_min_calc.append(((muon_eta[k]-muon_eta_gen[j])**2+(muon_phi[k]-muon_phi_gen[j])**2)**0.5) #Calculate dR, add it to the list
        
                    else:
                        continue
        
                muon_dR_min.append(min(muon_dR_min_calc))  # Take the min value of dR and add it to the list for matched muon dR's
                muon_dR_min_matched.append(min(muon_dR_min_calc))
                muon_matched_L1T_number.append(muon_dR_min_calc_fix.index(min(muon_dR_min_calc)))      
                # Find the index of the L1T muon, that gave the lowest dR value for the particular genP muon, and add it to the extra list.
                muon_dR_min_iev.append(iev)
                muon_dR_min_genP_i.append([muon_i_gen[j]])
            
            if muon_matched_L1T_number[0] == muon_matched_L1T_number[1]:
                print("Checking for repeated matching of the same muon, double matching... Need a fix")
                
            for n in range(len(muon_matched_L1T_number)):
                        
                muon_matched_L1T_charge.append(muon_charge[muon_matched_L1T_number[n]])    
                muon_matched_L1T_pt.append(muon_pt[muon_matched_L1T_number[n]])
                muon_matched_L1T_eta.append(muon_eta[muon_matched_L1T_number[n]])
                muon_matched_L1T_phi.append(muon_phi[muon_matched_L1T_number[n]])
                muon_matched_L1T_iev.append(muon_iev[muon_matched_L1T_number[n]])
                muon_matched_L1T_i.append(muon_i[muon_matched_L1T_number[n]])
                
                muon_dR_min_L1T_i.append(muon_i[muon_matched_L1T_number[n]])
                
            if len(muon_pdgId) != len(muon_matched_L1T_charge):
                print("Something is not gucci with the matching count")
                
            if len(muon_pdgId) == len(muon_matched_L1T_charge): # Checking whether all genP muons got matched with something. 
                
                if (muon_matched_L1T_charge[0]+muon_matched_L1T_charge[1] != 0):
                    print("Something is wrong with the matched muon charges or the numbering in the list")
                
                if len(muon_matched_L1T_charge) == 2: # Calculating the dR and invariant masses for L1T muon pairs
                
                    muon1 = ROOT.TLorentzVector()
                    muon2 = ROOT.TLorentzVector()

                    muon1.SetPtEtaPhiM(muon_matched_L1T_pt[0],
                                     muon_matched_L1T_eta[0],
                                     muon_matched_L1T_phi[0],0.10566)
                    muon2.SetPtEtaPhiM(muon_matched_L1T_pt[1],
                                     muon_matched_L1T_eta[1],
                                     muon_matched_L1T_phi[1],0.10566)
            
                    di_muon_inv_mass.append((muon1+muon2).M())
                    di_muon_dR_matched_L1T.append(((muon_matched_L1T_eta[0]-muon_matched_L1T_eta[1])**2+(muon_matched_L1T_phi[0]-muon_matched_L1T_phi[1])**2)**0.5)
                    inv_M_o_dR_di_muon.append(((muon1+muon2).M()/((muon_matched_L1T_eta[0]-muon_matched_L1T_eta[1])**2+(muon_matched_L1T_phi[0]-muon_matched_L1T_phi[1])**2)**0.5))
                else:
                    print("We have more than 2 matched L1T muons, probably 4, need to upgrade the code for dR and inv.M/dR calculation")
    
    # Here Im writing info to files.
            
    
    fon=open("Muon_dR_from_JPsi_genP.txt","a")
    for b in range(len(di_muon_dR_genP)):
        fon.write(str(di_muon_dR_genP[b])+"\t"+str(di_muon_dR_genP_iev[b])+"\t"+str(di_muon_dR_genP_i[b])+"\t"+str(aiziet)+"\n")
    fon.close()
    
    fob=open("Muon_dR_from_Matching_genP_L1T.txt","a")
    for h in range(len(muon_dR_min)):
        fob.write(str(muon_dR_min[h])+"\t"+str(muon_dR_min_iev[h])+"\t"+str(muon_dR_min_genP_i[h])+"\t"+str(muon_dR_min_L1T_i[h])+"\t"+str(aiziet)+"\n")
    fob.close()
    
    fov=open("Di_Muon_dR_invM_invM_o_dR.txt","a")
    for d in range(len(di_muon_dR_matched_L1T)):
        fov.write(str(di_muon_dR_matched_L1T[d])+"\t"+str(di_muon_inv_mass[d])+"\t"+str(inv_M_o_dR_di_muon[d])+"\t"+str(aiziet)+"\n")
    fov.close()                
                
    print("List of di-muon dR from J/Psi in genP is this long -> "+str(len(di_muon_dR_genP)))
    print("List of matched genP->L1T muons from pairs is this long -> "+str(len(muon_dR_min)))             
    print("List of di-muon dR from J/Psi in matched L1T muons is this long -> "+str(len(di_muon_dR_matched_L1T)))                           
                
print("We had {0} strange errors, with same parameters for muon and anti-muon".format(len(error_list)))
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
                              
                               
