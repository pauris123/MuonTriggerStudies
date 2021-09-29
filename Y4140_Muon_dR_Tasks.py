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

names=pn.read_csv("filenames_Y4140.txt", sep='\t',names=["nos"]) 
print("I have the libraries and names")

# File counting stuff
notik=0
tik=len(names["nos"])

error_list = [] # counting errors
shady_match_list = [] # list for counting how many shady matchings with 1 muon from J/Psi and other from Phi we get.

for aiziet in range(notik,tik):
    print("\n {0}/{1}, {2}".format(aiziet+1,tik,wall_time(time.time()-bigBang)))

    start=time.time()       
    rootfile="//cmsxrootd.fnal.gov///"+names["nos"][aiziet]
    events = Events("root:"+rootfile)
    #print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))
    
    
    
    di_muon_dR_JPsi_genP = [] # dR between muons that come from J/Psi to 2Mu
    di_muon_dR_JPsi_genP_iev = [] # List to help identify the Chi? (event) that this dR is calculated for
    di_muon_dR_JPsi_genP_i = [] # List to help identify the Chi? (i) that this dR is calculated for
    di_muon_invM_JPsi_genP = []
    di_muon_invM_o_dR_JPsi_genP = []
    
    di_muon_dR_Phi_genP = [] # dR between muons that come from Phi to 2Mu
    di_muon_dR_Phi_genP_iev = [] # List to help identify the Chi? (event) that this dR is calculated for
    di_muon_dR_Phi_genP_i = [] # List to help identify the Chi? (i) that this dR is calculated for
    di_muon_invM_Phi_genP = []
    di_muon_invM_o_dR_Phi_genP = []
    
    
    muon_dR_min_JPsi = [] # Matching genP and L1T muons by dR, the list of minimal/matching dR for J/Psi.
    muon_dR_min_JPsi_iev = []
    muon_dR_min_JPsi_genP_i = []
    
    muon_dR_min_Phi = [] # Matching genP and L1T muons by dR, the list of minimal/matching dR for Phi.
    muon_dR_min_Phi_iev = []
    muon_dR_min_Phi_genP_i = []
    
    
    muon_dR_min_JPsi_L1T_i = []
    muon_dR_min_Phi_L1T_i = []
    
    
    di_muon_inv_mass_JPsi = []        # Lists for dR, inv.M and in.M/dR for matched L1T di-muons reconstructed from Chi decay 
    di_muon_dR_matched_L1T_JPsi = []
    inv_M_o_dR_di_muon_JPsi = []
    
    di_muon_inv_mass_Phi = []        # Lists for dR, inv.M and in.M/dR for matched L1T di-muons reconstructed from Chi decay
    di_muon_dR_matched_L1T_Phi = []
    inv_M_o_dR_di_muon_Phi = []
            
    for iev,event in enumerate(events): # Reading in events from the file
        event.getByLabel(l1MuonLabel, l1Muons)
        event.getByLabel(genParticlesLabel, genParticles)

        muon_charge = [] # Lists for read in L1T muons
        muon_pt = []
        muon_eta = []
        muon_phi = []
        muon_iev = []
        muon_i = []
    
        #muon_pdgId = [] # Lists for decayed muons from Chi->J/Psi->MuMu in GenP  (1st task)!
        #muon_pt_gen = []
        #muon_eta_gen = []
        #muon_phi_gen = []
        #muon_iev_gen = []
        #muon_i_gen = []
        
        muon_JPsi_pdgId = [] # Lists for decayed muons from J/Psi to 2 Mu, for redundancy
        muon_pt_JPsi_gen = []
        muon_eta_JPsi_gen = []
        muon_phi_JPsi_gen = []
        muon_iev_JPsi_gen = []
        muon_i_JPsi_gen = []
        
        muon_Phi_pdgId = [] # Lists for decayed muons from Phi to 2 Mu, for redundancy
        muon_pt_Phi_gen = []
        muon_eta_Phi_gen = []
        muon_phi_Phi_gen = []
        muon_iev_Phi_gen = []
        muon_i_Phi_gen = []
        
        muon_dR_min_matched_JPsi = []
        muon_matched_JPsi_L1T_number = []
        
        muon_dR_min_matched_Phi = []  # Same as muon_dR_min, but for each event, makes writing data more simpler in the next if and for cycles
        muon_matched_Phi_L1T_number = []
                # Matched L1T muon parameters to genP muons. 1st L1T matched muon in this list corresponds to the first muon in genP list 
    # So the list might go - (10th L1T muon -> 1st GenP muon), (24th L1T muon -> 2nd GenP muon) and so on, so have to follow the iev and i numbers of L1T muon
    
        muon_matched_JPsi_L1T_charge = []  # Lists of matched J/Psi L1T muon parameters, so it is easier to calculate stuff onwards
        muon_matched_JPsi_L1T_pt = []
        muon_matched_JPsi_L1T_eta = []
        muon_matched_JPsi_L1T_phi = []
        muon_matched_JPsi_L1T_iev = []
        muon_matched_JPsi_L1T_i = []
        
        muon_matched_Phi_L1T_charge = []  # Lists of matched Phi L1T muon parameters, so it is easier to calculate stuff onwards
        muon_matched_Phi_L1T_pt = []
        muon_matched_Phi_L1T_eta = []
        muon_matched_Phi_L1T_phi = []
        muon_matched_Phi_L1T_iev = []
        muon_matched_Phi_L1T_i = []
        
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
            if (genP.pdgId() == 100443 and genP.numberOfDaughters() == 2) : # Takes Mother particle that has 2 daughters
                
                #if genP.daughter(0).pdgId() == 443: # Checing whether the first daughter is J/Psi or Phi
                    #print("J/Psi is 1st daughter")
                #elif genP.daughter(0).pdgId() == 333:
                    #print("Phi is 1st daughter")
                
                #if genP.daughter(1).pdgId() == 443: # Checing whether the second daughter is J/Psi or Phi
                    #print("J/Psi is 2nd daughter")
                #elif genP.daughter(1).pdgId() == 333:
                    #print("Phi is 2nd daughter")
                
                if genP.daughter(0).numberOfDaughters() > 0 and genP.daughter(0).pdgId() == 443: # Cheks for grand-daughters
                    
                    muon_eta_for_genP_dR_0 = [] # Lists for calculating di-muon dR that come from 1st daughter decay in GenParticles
                    muon_phi_for_genP_dR_0 = []
                    muon_pt_for_genP_dR_0 = []
                    
                    for y in range(genP.daughter(0).numberOfDaughters()):
                        # Prints J/Psi or photon decay products. In this case muons or same strange photons

                        if (genP.daughter(0).daughter(y).pdgId() == 13 or genP.daughter(0).daughter(y).pdgId() == -13): # Cheks for decay muons
                            
                            #Takes the decay Mu parameters and adds them to the list
                            #muon_pdgId.append(genP.daughter(0).daughter(y).pdgId())
                            #muon_pt_gen.append(genP.daughter(0).daughter(y).pt())
                            #muon_eta_gen.append(genP.daughter(0).daughter(y).eta())
                            #muon_phi_gen.append(genP.daughter(0).daughter(y).phi())
                            #muon_iev_gen.append(iev)
                            #muon_i_gen.append(i)
                            
                            muon_JPsi_pdgId.append(genP.daughter(0).daughter(y).pdgId()) # Lists for decayed muons from Phi to 2 Mu, for redundancy
                            muon_pt_JPsi_gen.append(genP.daughter(0).daughter(y).pt())
                            muon_eta_JPsi_gen.append(genP.daughter(0).daughter(y).eta())
                            muon_phi_JPsi_gen.append(genP.daughter(0).daughter(y).phi())
                            muon_iev_JPsi_gen.append(iev)
                            muon_i_JPsi_gen.append(i)
                            
                            #List For calclating dR between 2 muons coming from J/Psi
                            muon_eta_for_genP_dR_0.append(genP.daughter(0).daughter(y).eta())
                            muon_phi_for_genP_dR_0.append(genP.daughter(0).daughter(y).phi())
                            muon_pt_for_genP_dR_0.append(genP.daughter(0).daughter(y).pt())
                            
                    if len(muon_eta_for_genP_dR_0) == 2: # Check whether there is Di-muon pair from J/Psi to calculate dR for, from GenP
                        
                        phi_check_JPsi_genP = []
                        
                        if abs(muon_phi_for_genP_dR_0[0]-muon_phi_for_genP_dR_0[1]) > 3.1415926:
                            phi_check_JPsi_genP.append((6.2831852-abs(muon_phi_for_genP_dR_0[0]-muon_phi_for_genP_dR_0[1])))
                        else:
                            phi_check_JPsi_genP.append((muon_phi_for_genP_dR_0[0]-muon_phi_for_genP_dR_0[1]))                        
                        
                        
                        di_muon_dR_JPsi_genP.append(((muon_eta_for_genP_dR_0[0]-muon_eta_for_genP_dR_0[1])**2+(phi_check_JPsi_genP[0])**2)**0.5)
                        di_muon_dR_JPsi_genP_iev.append(iev)
                        di_muon_dR_JPsi_genP_i.append(i)
            
                        muon1_JPsi_genP = ROOT.TLorentzVector()
                        muon2_JPsi_genP = ROOT.TLorentzVector()

                        muon1_JPsi_genP.SetPtEtaPhiM(muon_pt_for_genP_dR_0[0],
                                         muon_eta_for_genP_dR_0[0],
                                         muon_phi_for_genP_dR_0[0],0.10566)
                        muon2_JPsi_genP.SetPtEtaPhiM(muon_pt_for_genP_dR_0[1],
                                         muon_eta_for_genP_dR_0[1],
                                         muon_phi_for_genP_dR_0[1],0.10566)
            
                        di_muon_invM_JPsi_genP.append((muon1_JPsi_genP+muon2_JPsi_genP).M())
                        di_muon_invM_o_dR_JPsi_genP.append(((muon1_JPsi_genP+muon2_JPsi_genP).M()/((muon_eta_for_genP_dR_0[0]-muon_eta_for_genP_dR_0[1])**2+(phi_check_JPsi_genP[0])**2)**0.5))
                
                if genP.daughter(1).numberOfDaughters() > 0 and genP.daughter(1).pdgId() == 333: # Checks all the same things for the 2nd daughter and following grand-daughters
                    
                    muon_eta_for_genP_dR_1 = [] # Lists for calculating dR between both muons that come from 1st daughter decay in GenParticles (should be Phi)
                    muon_phi_for_genP_dR_1 = []
                    muon_pt_for_genP_dR_1 = []
                    
                    for y in range(genP.daughter(1).numberOfDaughters()):
                        
                        if (genP.daughter(1).daughter(y).pdgId() == 13 or genP.daughter(1).daughter(y).pdgId() == -13):
                            
                            #muon_pdgId.append(genP.daughter(1).daughter(y).pdgId())
                            #muon_pt_gen.append(genP.daughter(1).daughter(y).pt())
                            #muon_eta_gen.append(genP.daughter(1).daughter(y).eta())
                            #muon_phi_gen.append(genP.daughter(1).daughter(y).phi())
                            #muon_iev_gen.append(iev)
                            #muon_i_gen.append(i)
                            
                            muon_Phi_pdgId.append(genP.daughter(1).daughter(y).pdgId()) # Lists for decayed muons from J/Psi to 2 Mu, for redundancy
                            muon_pt_Phi_gen.append(genP.daughter(1).daughter(y).pt())
                            muon_eta_Phi_gen.append(genP.daughter(1).daughter(y).eta())
                            muon_phi_Phi_gen.append(genP.daughter(1).daughter(y).phi())
                            muon_iev_Phi_gen.append(iev)
                            muon_i_Phi_gen.append(i)
                            
                            muon_eta_for_genP_dR_1.append(genP.daughter(1).daughter(y).eta())
                            muon_phi_for_genP_dR_1.append(genP.daughter(1).daughter(y).phi())
                            muon_pt_for_genP_dR_1.append(genP.daughter(1).daughter(y).pt())
                            
                    if len(muon_eta_for_genP_dR_1) == 2:
                        
                        phi_check_Phi_genP = []
                        
                        if abs(muon_phi_for_genP_dR_1[0]-muon_phi_for_genP_dR_1[1]) > 3.1415926:
                            phi_check_Phi_genP.append((6.2831852-abs(muon_phi_for_genP_dR_1[0]-muon_phi_for_genP_dR_1[1])))
                        else:
                            phi_check_Phi_genP.append((muon_phi_for_genP_dR_1[0]-muon_phi_for_genP_dR_1[1]))                        
                        
                        
                        di_muon_dR_Phi_genP.append(((muon_eta_for_genP_dR_1[0]-muon_eta_for_genP_dR_1[1])**2+(phi_check_Phi_genP[0])**2)**0.5)
                        di_muon_dR_Phi_genP_iev.append(iev)
                        di_muon_dR_Phi_genP_i.append(i)
        
                        muon1_Phi_genP = ROOT.TLorentzVector()
                        muon2_Phi_genP = ROOT.TLorentzVector()

                        muon1_Phi_genP.SetPtEtaPhiM(muon_pt_for_genP_dR_0[0],
                                         muon_eta_for_genP_dR_0[0],
                                         muon_phi_for_genP_dR_0[0],0.10566)
                        muon2_Phi_genP.SetPtEtaPhiM(muon_pt_for_genP_dR_0[1],
                                         muon_eta_for_genP_dR_0[1],
                                         muon_phi_for_genP_dR_0[1],0.10566)
            
                        di_muon_invM_Phi_genP.append((muon1_Phi_genP+muon2_Phi_genP).M())
                        di_muon_invM_o_dR_Phi_genP.append(((muon1_Phi_genP+muon2_Phi_genP).M()/((muon_eta_for_genP_dR_1[0]-muon_eta_for_genP_dR_1[1])**2+(phi_check_Phi_genP[0])**2)**0.5))                
    

            
    
        #Calculating the dR values for each genP muon matched with every L1T muon
        if (len(muon_JPsi_pdgId) == 2 and len(muon_Phi_pdgId) == 2 and len(muon_charge) >= 2) :  # (len(muon_pdgId) == 2 or len(muon_pdgId) == 4)
            if len(muon_charge) == 2 and sum(muon_charge) != 0:
                continue
            if len(muon_charge) == 3 and abs(sum(muon_charge)) == 3:
                continue
            if len(muon_charge) == 4 and abs(sum(muon_charge)) == 4 :
                continue
            if len(muon_charge) == 5 and abs(sum(muon_charge)) == 5 :
                continue
            if len(muon_charge) == 6 and abs(sum(muon_charge)) == 6 :
                continue
            if (muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] != muon_charge[1]):
                print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                error_list.append(i)
                continue            
            if (len(muon_charge) == 3 and ((muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] != muon_charge[1]) or (muon_eta[0] == muon_eta[2] and muon_phi[0] == muon_phi[2] and muon_charge[0] != muon_charge[2]) or (muon_eta[2] == muon_eta[1] and muon_phi[2] == muon_phi[1] and muon_charge[2] != muon_charge[1]))):
                print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                error_list.append(i)
                continue
            if (len(muon_charge) == 4 and ((muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] != muon_charge[1]) or \
            (muon_eta[0] == muon_eta[2] and muon_phi[0] == muon_phi[2] and muon_charge[0] != muon_charge[2]) or \
            (muon_eta[0] == muon_eta[3] and muon_phi[0] == muon_phi[3] and muon_charge[0] != muon_charge[3]) or \
            (muon_eta[1] == muon_eta[2] and muon_phi[1] == muon_phi[2] and muon_charge[1] != muon_charge[2]) or \
            (muon_eta[1] == muon_eta[3] and muon_phi[1] == muon_phi[3] and muon_charge[1] != muon_charge[3]) or \
            (muon_eta[2] == muon_eta[3] and muon_phi[2] == muon_phi[3] and muon_charge[2] != muon_charge[3]))):
                print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                error_list.append(i)
                continue
            
            
            muon_matched_L1T_JPsi_number_matched = []
            muon_dR_min_JPsi_genP_i_matched = []            
            muon_matched_L1T_Phi_number_matched = []
            muon_dR_min_Phi_genP_i_matched = []
                        
            for j in range(len(muon_JPsi_pdgId)): # Check matching dR for J/Psi muons
    
                muon_dR_min_calc = []
                muon_dR_min_calc_fix = [] # Fixed problem with matched L1T numbers
    
                for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
                    
                    phi_check_JPsi_matching = []
                        
                    if abs(muon_phi[k]-muon_phi_JPsi_gen[j]) > 3.1415926:
                        phi_check_JPsi_matching.append((6.2831852-abs(muon_phi[k]-muon_phi_JPsi_gen[j])))
                    else:
                        phi_check_JPsi_matching.append((muon_phi[k]-muon_phi_JPsi_gen[j]))
                                            
                    
                    muon_dR_min_calc_fix.append(((muon_eta[k]-muon_eta_JPsi_gen[j])**2+(phi_check_JPsi_matching[0])**2)**0.5) #Temporary fix
                    
                    if (muon_charge[k] == -1 and muon_JPsi_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_JPsi_pdgId[j] == -13):  
                
                        muon_dR_min_calc.append(((muon_eta[k]-muon_eta_JPsi_gen[j])**2+(phi_check_JPsi_matching[0])**2)**0.5) #Calculate dR, add it to the list
        
                    else:
                        continue
                
                    # J/Psi Matching
                muon_dR_min_matched_JPsi.append(min(muon_dR_min_calc))
                muon_matched_L1T_JPsi_number_matched.append(muon_dR_min_calc_fix.index(min(muon_dR_min_calc)))
                muon_dR_min_JPsi_genP_i_matched.append([muon_i_JPsi_gen[j]])                                  
            
            
            for j in range(len(muon_Phi_pdgId)): # Check matching dR for Phi muons
    
                muon_dR_min_calc = []
                muon_dR_min_calc_fix = [] # Fixed problem with matched L1T numbers
    
                for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
                    
                    phi_check_Phi_matching = []
                        
                    if abs(muon_phi[k]-muon_phi_Phi_gen[j]) > 3.1415926:
                        phi_check_Phi_matching.append((6.2831852-abs(muon_phi[k]-muon_phi_Phi_gen[j])))
                    else:
                        phi_check_Phi_matching.append((muon_phi[k]-muon_phi_Phi_gen[j]))                    
                    
                    muon_dR_min_calc_fix.append(((muon_eta[k]-muon_eta_Phi_gen[j])**2+(phi_check_Phi_matching[0])**2)**0.5) #Temporary fix
                    
                    if (muon_charge[k] == -1 and muon_Phi_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_Phi_pdgId[j] == -13):  
                
                        muon_dR_min_calc.append(((muon_eta[k]-muon_eta_Phi_gen[j])**2+(phi_check_Phi_matching[0])**2)**0.5) #Calculate dR, add it to the list
        
                    else:
                        continue                

                    # Phi Matching
                muon_dR_min_matched_Phi.append(min(muon_dR_min_calc))
                muon_matched_L1T_Phi_number_matched.append(muon_dR_min_calc_fix.index(min(muon_dR_min_calc)))
                muon_dR_min_Phi_genP_i_matched.append([muon_i_Phi_gen[j]])                
                
            
            #print("JPsi matching dR -> "+str(muon_dR_min_matched_JPsi))
            #print("Phi matching dR -> "+str(muon_dR_min_matched_Phi))
            
            if len(muon_charge) == 2 or len(muon_charge) == 3: # When we can use only the best fitting pair out of J/Psi or Phi
            # We chose which pair fits better
                if (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) < (muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]):
                # We can add matching dR cut here
                    for l in range(len(muon_dR_min_matched_JPsi)):
                        
                        muon_dR_min_JPsi.append(muon_dR_min_matched_JPsi[l])
                        muon_matched_JPsi_L1T_number.append(muon_matched_L1T_JPsi_number_matched[l])   
                        muon_dR_min_JPsi_iev.append(iev)
                        muon_dR_min_JPsi_genP_i.append(muon_dR_min_JPsi_genP_i_matched[l])
                    
                    if muon_matched_JPsi_L1T_number[0] == muon_matched_JPsi_L1T_number[1]:
                        print("Checking for repeated matching of the same muon, double matching... Need a fix (2-3 muons)")
                    
                    #if ((muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) > (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_Phi[0])) or \
                    #((muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) > (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_Phi[1])) or \
                    #((muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) > (muon_dR_min_matched_JPsi[1]+muon_dR_min_matched_Phi[0])) or \
                    #((muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) > (muon_dR_min_matched_JPsi[1]+muon_dR_min_matched_Phi[1])):
                        # Need to check the "if" statement above, seems a bit shady
                        #print("One muon from J/Psi and one from Phi fits better in terms of dR match")
                    if (muon_dR_min_matched_JPsi[0] > 2 or muon_dR_min_matched_JPsi[1] > 2):
                        shady_match_list.append(l)
                        del muon_dR_min_JPsi[-1]
                        del muon_dR_min_JPsi[-1]
                        del muon_matched_JPsi_L1T_number[-1]
                        del muon_matched_JPsi_L1T_number[-1]
                        del muon_dR_min_JPsi_iev[-1]
                        del muon_dR_min_JPsi_iev[-1]
                        del muon_dR_min_JPsi_genP_i[-1]
                        del muon_dR_min_JPsi_genP_i[-1]
            
                else:
                # We can add matching dR cut here
                    for w in range(len(muon_dR_min_matched_Phi)):
                    
                        muon_dR_min_Phi.append(muon_dR_min_matched_Phi[w])
                        muon_matched_Phi_L1T_number.append(muon_matched_L1T_Phi_number_matched[w])   
                        muon_dR_min_Phi_iev.append(iev)
                        muon_dR_min_Phi_genP_i.append(muon_dR_min_Phi_genP_i_matched[w])   
                        
                    if muon_matched_Phi_L1T_number[0] == muon_matched_Phi_L1T_number[1]:
                        print("Checking for repeated matching of the same muon, double matching... Need a fix (2-3 muons)")                                     
            
                    #if ((muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]) > (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_Phi[0])) or \
                    #((muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]) > (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_Phi[1])) or \
                    #((muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]) > (muon_dR_min_matched_JPsi[1]+muon_dR_min_matched_Phi[0])) or \
                    #((muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]) > (muon_dR_min_matched_JPsi[1]+muon_dR_min_matched_Phi[1])):
                        #print("One muon from J/Psi and one from Phi fits better in terms of dR match")
                        # Need to check the "if" statement above, seems a bit shady
                    if (muon_dR_min_matched_Phi[0] > 2 or muon_dR_min_matched_Phi[1] > 2):    
                        shady_match_list.append(w)
                        del muon_dR_min_Phi[-1]
                        del muon_dR_min_Phi[-1]
                        del muon_matched_Phi_L1T_number[-1]
                        del muon_matched_Phi_L1T_number[-1]
                        del muon_dR_min_Phi_iev[-1]
                        del muon_dR_min_Phi_iev[-1]
                        del muon_dR_min_Phi_genP_i[-1]
                        del muon_dR_min_Phi_genP_i[-1]            
            
            else: # If I can find all 4 muons coming from J/Psi or Phi
                #taisam matching dR < 2 cut
                if muon_dR_min_matched_JPsi[0] > 2 or muon_dR_min_matched_JPsi[1] > 2:
                    if muon_dR_min_matched_Phi[0] > 2 or muon_dR_min_matched_Phi[1] > 2:
                        continue
                    else:
                        for l in range(len(muon_dR_min_matched_Phi)):             
            
                            muon_dR_min_Phi.append(muon_dR_min_matched_Phi[l])
                            muon_matched_Phi_L1T_number.append(muon_matched_L1T_Phi_number_matched[l])   
                            muon_dR_min_Phi_iev.append(iev)
                            muon_dR_min_Phi_genP_i.append(muon_dR_min_Phi_genP_i_matched[l])
                        print("dR cut condition applied for 4muon case")
                        continue
                
                if muon_dR_min_matched_Phi[0] > 2 or muon_dR_min_matched_Phi[1] > 2:
                    if muon_dR_min_matched_JPsi[0] > 2 or muon_dR_min_matched_JPsi[1] > 2:
                        continue
                    else:
                        for l in range(len(muon_dR_min_matched_JPsi)):             
            
                            muon_dR_min_JPsi.append(muon_dR_min_matched_JPsi[l])
                            muon_matched_JPsi_L1T_number.append(muon_matched_L1T_JPsi_number_matched[l])   
                            muon_dR_min_JPsi_iev.append(iev)
                            muon_dR_min_JPsi_genP_i.append(muon_dR_min_JPsi_genP_i_matched[l]) 
                        print("dR cut condition applied for 4muon case")
                        continue                
                
                for l in range(len(muon_dR_min_matched_JPsi)):
                        
                    muon_dR_min_JPsi.append(muon_dR_min_matched_JPsi[l])
                    muon_matched_JPsi_L1T_number.append(muon_matched_L1T_JPsi_number_matched[l])   
                    muon_dR_min_JPsi_iev.append(iev)
                    muon_dR_min_JPsi_genP_i.append(muon_dR_min_JPsi_genP_i_matched[l])                
            
                    muon_dR_min_Phi.append(muon_dR_min_matched_Phi[l])
                    muon_matched_Phi_L1T_number.append(muon_matched_L1T_Phi_number_matched[l])   
                    muon_dR_min_Phi_iev.append(iev)
                    muon_dR_min_Phi_genP_i.append(muon_dR_min_Phi_genP_i_matched[l])
                    
                if (muon_matched_JPsi_L1T_number[0] == muon_matched_JPsi_L1T_number[1]) or (muon_matched_Phi_L1T_number[0] == muon_matched_Phi_L1T_number[1]) or \
                (muon_matched_JPsi_L1T_number[0] == muon_matched_Phi_L1T_number[0]) or (muon_matched_JPsi_L1T_number[0] == muon_matched_Phi_L1T_number[1]) or \
                (muon_matched_JPsi_L1T_number[1] == muon_matched_Phi_L1T_number[0]) or (muon_matched_JPsi_L1T_number[1] == muon_matched_Phi_L1T_number[1]):
                    print("Checking for repeated matching of the same muon, double matching... Need a fix (4 muons)")
                    print(str(muon_matched_JPsi_L1T_number)+" "+str(muon_matched_Phi_L1T_number))
                    print(str(muon_dR_min_matched_JPsi)+" "+str(muon_dR_min_matched_Phi))
                    
                    if (muon_dR_min_matched_JPsi[0]+muon_dR_min_matched_JPsi[1]) < (muon_dR_min_matched_Phi[0]+muon_dR_min_matched_Phi[1]):
                        # J/Psi fits better, need to refit Phi
                        del muon_dR_min_Phi[-1]
                        del muon_dR_min_Phi[-1]
                        del muon_matched_Phi_L1T_number[-1]
                        del muon_matched_Phi_L1T_number[-1]
                        del muon_dR_min_Phi_genP_i[-1]
                        del muon_dR_min_Phi_genP_i[-1]                    
                    
                        if ((len(muon_charge) == 4 and abs(sum(muon_charge)) == 2) or (len(muon_charge)==5 and abs(sum(muon_charge)) == 3) or (len(muon_charge)==6 and abs(sum(muon_charge)) == 4)):
                            print("Could fit only J/Psi, Phi discarded")
                            continue
                        
                        if (len(muon_charge) == 4 and ((muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] == muon_charge[1]) or \
                        (muon_eta[0] == muon_eta[2] and muon_phi[0] == muon_phi[2] and muon_charge[0] == muon_charge[2]) or \
                        (muon_eta[0] == muon_eta[3] and muon_phi[0] == muon_phi[3] and muon_charge[0] == muon_charge[3]) or \
                        (muon_eta[1] == muon_eta[2] and muon_phi[1] == muon_phi[2] and muon_charge[1] == muon_charge[2]) or \
                        (muon_eta[1] == muon_eta[3] and muon_phi[1] == muon_phi[3] and muon_charge[1] == muon_charge[3]) or \
                        (muon_eta[2] == muon_eta[3] and muon_phi[2] == muon_phi[3] and muon_charge[2] == muon_charge[3]))):
                            print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                            error_list.append(i)
                            continue                        
                        
                        
                        
                        muon_matched_L1T_Phi_number_matched_4mu = []
                        
                        muon_numbers = range(len(muon_charge))
                        muon_numbers.remove(muon_matched_L1T_JPsi_number_matched[0])
                        muon_numbers.remove(muon_matched_L1T_JPsi_number_matched[1])
                                                
                        for j in range(len(muon_Phi_pdgId)): # Check matching dR for Phi muons
                            
                            muon_dR_min_calc_4mu = []
                            muon_dR_min_calc_fix_4mu = [] # Fixed problem with matched L1T numbers
                                
                            for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
                    
                                phi_check_Phi_matching_4mu = []
                        
                                if abs(muon_phi[k]-muon_phi_Phi_gen[j]) > 3.1415926:
                                    phi_check_Phi_matching_4mu.append((6.2831852-abs(muon_phi[k]-muon_phi_Phi_gen[j])))
                                else:
                                    phi_check_Phi_matching_4mu.append((muon_phi[k]-muon_phi_Phi_gen[j]))                    
                    
                                muon_dR_min_calc_fix_4mu.append(((muon_eta[k]-muon_eta_Phi_gen[j])**2+(phi_check_Phi_matching_4mu[0])**2)**0.5) #Temporary fix
                                
                                if k in muon_numbers:
                    
                                    if (muon_charge[k] == -1 and muon_Phi_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_Phi_pdgId[j] == -13):  
                
                                        muon_dR_min_calc_4mu.append(((muon_eta[k]-muon_eta_Phi_gen[j])**2+(phi_check_Phi_matching_4mu[0])**2)**0.5) 
                                        #Calculate dR, add it to the list
                                else:
                                    continue                      
                            muon_matched_L1T_Phi_number_matched_4mu.append(muon_dR_min_calc_fix_4mu.index(min(muon_dR_min_calc_4mu)))
                            muon_dR_min_Phi.append(min(muon_dR_min_calc_4mu))
                            muon_matched_Phi_L1T_number.append(muon_dR_min_calc_fix_4mu.index(min(muon_dR_min_calc_4mu)))
                            muon_dR_min_Phi_genP_i.append([muon_i_Phi_gen[j]])
                        print("J/Psi is fine. Now it should be fixed for Phi")
                        
                        print("Fixed nubemrs"+str(muon_matched_JPsi_L1T_number)+" "+str(muon_matched_Phi_L1T_number))
                        if len(muon_matched_L1T_Phi_number_matched_4mu) != 2:
                            print("We have rematched less or more than 2 muons in 4Mu rematching part of the code...")


                    
                    else: # Phi fits better, need to refit J/Psi
                        del muon_dR_min_JPsi[-1]
                        del muon_dR_min_JPsi[-1]
                        del muon_matched_JPsi_L1T_number[-1]
                        del muon_matched_JPsi_L1T_number[-1]
                        del muon_dR_min_JPsi_genP_i[-1]
                        del muon_dR_min_JPsi_genP_i[-1]                

                        if ((len(muon_charge) == 4 and abs(sum(muon_charge)) == 2) or (len(muon_charge)==5 and abs(sum(muon_charge)) == 3) or (len(muon_charge)==6 and abs(sum(muon_charge)) == 4)):
                            print("Could fit only Phi, JPsi discarded")
                            continue
                        
                        if (len(muon_charge) == 4 and ((muon_eta[0] == muon_eta[1] and muon_phi[0] == muon_phi[1] and muon_charge[0] == muon_charge[1]) or \
                        (muon_eta[0] == muon_eta[2] and muon_phi[0] == muon_phi[2] and muon_charge[0] == muon_charge[2]) or \
                        (muon_eta[0] == muon_eta[3] and muon_phi[0] == muon_phi[3] and muon_charge[0] == muon_charge[3]) or \
                        (muon_eta[1] == muon_eta[2] and muon_phi[1] == muon_phi[2] and muon_charge[1] == muon_charge[2]) or \
                        (muon_eta[1] == muon_eta[3] and muon_phi[1] == muon_phi[3] and muon_charge[1] == muon_charge[3]) or \
                        (muon_eta[2] == muon_eta[3] and muon_phi[2] == muon_phi[3] and muon_charge[2] == muon_charge[3]))):
                            print("I did an upsi with same parameters for muon and anti-muon #Rebel")
                            error_list.append(i)
                            continue                        
                        
                        muon_matched_L1T_JPsi_number_matched_4mu = []
                        
                        muon_numbers = range(len(muon_charge))
                        muon_numbers.remove(muon_matched_L1T_Phi_number_matched[0])
                        muon_numbers.remove(muon_matched_L1T_Phi_number_matched[1])
                                                
                        for j in range(len(muon_JPsi_pdgId)): # Check matching dR for Phi muons
    
                            muon_dR_min_calc_4mu = []
                            muon_dR_min_calc_fix_4mu = [] # Fixed problem with matched L1T numbers
                                
                            for k in range(len(muon_charge)): # Check the the chosen genP muon with all the usable L1T muons
                    
                                phi_check_JPsi_matching_4mu = []
                        
                                if abs(muon_phi[k]-muon_phi_JPsi_gen[j]) > 3.1415926:
                                    phi_check_JPsi_matching_4mu.append((6.2831852-abs(muon_phi[k]-muon_phi_JPsi_gen[j])))
                                else:
                                    phi_check_JPsi_matching_4mu.append((muon_phi[k]-muon_phi_JPsi_gen[j]))                    
                    
                                muon_dR_min_calc_fix_4mu.append(((muon_eta[k]-muon_eta_JPsi_gen[j])**2+(phi_check_JPsi_matching_4mu[0])**2)**0.5) #Temporary fix
                    
                                if k in muon_numbers:    
                                    
                                    if (muon_charge[k] == -1 and muon_JPsi_pdgId[j] == 13) or (muon_charge[k] == 1 and muon_JPsi_pdgId[j] == -13):  
                
                                        muon_dR_min_calc_4mu.append(((muon_eta[k]-muon_eta_JPsi_gen[j])**2+(phi_check_JPsi_matching_4mu[0])**2)**0.5) 
                                        #Calculate dR, add it to the list        
                                else:
                                    continue                      
                            muon_matched_L1T_JPsi_number_matched_4mu.append(muon_dR_min_calc_fix_4mu.index(min(muon_dR_min_calc_4mu)))
                            muon_dR_min_JPsi.append(min(muon_dR_min_calc_4mu))
                            muon_matched_JPsi_L1T_number.append(muon_dR_min_calc_fix_4mu.index(min(muon_dR_min_calc_4mu)))
                            muon_dR_min_JPsi_genP_i.append([muon_i_Phi_gen[j]])
                        print("Phi is fine. Now it should be fixed for J/Psi")

                        print("Fixed nubemrs"+str(muon_matched_JPsi_L1T_number)+" "+str(muon_matched_Phi_L1T_number))
                        if len(muon_matched_L1T_JPsi_number_matched_4mu) != 2:
                            print("We have rematched less or more than 2 muons in 4Mu rematching part of the code...")            
            
            
            if len(muon_matched_JPsi_L1T_number) == 2:
                for m in range(len(muon_matched_JPsi_L1T_number)): # Matching JPsi from L1T muons
                        
                    muon_matched_JPsi_L1T_charge.append(muon_charge[muon_matched_JPsi_L1T_number[m]])    
                    muon_matched_JPsi_L1T_pt.append(muon_pt[muon_matched_JPsi_L1T_number[m]])
                    muon_matched_JPsi_L1T_eta.append(muon_eta[muon_matched_JPsi_L1T_number[m]])
                    muon_matched_JPsi_L1T_phi.append(muon_phi[muon_matched_JPsi_L1T_number[m]])
                    muon_matched_JPsi_L1T_iev.append(muon_iev[muon_matched_JPsi_L1T_number[m]])
                    muon_matched_JPsi_L1T_i.append(muon_i[muon_matched_JPsi_L1T_number[m]])
                
                    #muon_dR_min_JPsi_L1T_i.append(muon_i[muon_matched_L1T_number[m]])
                
            if len(muon_matched_Phi_L1T_number) == 2:
                for n in range(len(muon_matched_Phi_L1T_number)): # Matching Phi from L1T muons
                        
                    muon_matched_Phi_L1T_charge.append(muon_charge[muon_matched_Phi_L1T_number[n]])    
                    muon_matched_Phi_L1T_pt.append(muon_pt[muon_matched_Phi_L1T_number[n]])
                    muon_matched_Phi_L1T_eta.append(muon_eta[muon_matched_Phi_L1T_number[n]])
                    muon_matched_Phi_L1T_phi.append(muon_phi[muon_matched_Phi_L1T_number[n]])
                    muon_matched_Phi_L1T_iev.append(muon_iev[muon_matched_Phi_L1T_number[n]])
                    muon_matched_Phi_L1T_i.append(muon_i[muon_matched_Phi_L1T_number[n]])
                
                    #muon_dR_min_Phi_L1T_i.append(muon_i[muon_matched_L1T_number[n]])    
                
            #if len(muon_pdgId) != (len(muon_matched_JPsi_L1T_charge)+len(muon_matched_Phi_L1T_charge)):
                #print("Something is not gucci with the matching count")
                
            if ((len(muon_matched_JPsi_L1T_charge)+len(muon_matched_Phi_L1T_charge)) == 2) or ((len(muon_matched_JPsi_L1T_charge)+len(muon_matched_Phi_L1T_charge)) == 4): 
            # Checking whether all genP muons got matched with something. 
                
                
                if len(muon_matched_JPsi_L1T_charge) == 2: # Calculating the dR and invariant masses for L1T muon pairs
                
                    if (muon_matched_JPsi_L1T_charge[0]+muon_matched_JPsi_L1T_charge[1] != 0):
                        print("Something is wrong with the matched J/Psi muon charges or the numbering in the list")                    
                    
                    phi_check_JPsi_L1T = []
                        
                    if abs(muon_matched_JPsi_L1T_phi[0]-muon_matched_JPsi_L1T_phi[1]) > 3.1415926:
                        phi_check_JPsi_L1T.append((6.2831852-abs(muon_matched_JPsi_L1T_phi[0]-muon_matched_JPsi_L1T_phi[1])))
                    else:
                        phi_check_JPsi_L1T.append((muon_matched_JPsi_L1T_phi[0]-muon_matched_JPsi_L1T_phi[1]))                         
                    
                    
                    muon1_JPsi = ROOT.TLorentzVector()
                    muon2_JPsi = ROOT.TLorentzVector()

                    muon1_JPsi.SetPtEtaPhiM(muon_matched_JPsi_L1T_pt[0],
                                     muon_matched_JPsi_L1T_eta[0],
                                     muon_matched_JPsi_L1T_phi[0],0.10566)
                    muon2_JPsi.SetPtEtaPhiM(muon_matched_JPsi_L1T_pt[1],
                                     muon_matched_JPsi_L1T_eta[1],
                                     muon_matched_JPsi_L1T_phi[1],0.10566)
            
                    di_muon_inv_mass_JPsi.append((muon1_JPsi+muon2_JPsi).M())
                    di_muon_dR_matched_L1T_JPsi.append(((muon_matched_JPsi_L1T_eta[0]-muon_matched_JPsi_L1T_eta[1])**2+(phi_check_JPsi_L1T[0])**2)**0.5)
                    inv_M_o_dR_di_muon_JPsi.append(((muon1_JPsi+muon2_JPsi).M()/((muon_matched_JPsi_L1T_eta[0]-muon_matched_JPsi_L1T_eta[1])**2+(phi_check_JPsi_L1T[0])**2)**0.5))

                
                if len(muon_matched_Phi_L1T_charge) == 2: # Calculating the dR and invariant masses for L1T muon pairs
                
                    if (muon_matched_Phi_L1T_charge[0]+muon_matched_Phi_L1T_charge[1] != 0):
                        print("Something is wrong with the matched Phi muon charges or the numbering in the list")                
                
                    phi_check_Phi_L1T = []
                        
                    if abs(muon_matched_Phi_L1T_phi[0]-muon_matched_Phi_L1T_phi[1]) > 3.1415926:
                        phi_check_Phi_L1T.append((6.2831852-abs(muon_matched_Phi_L1T_phi[0]-muon_matched_Phi_L1T_phi[1])))
                    else:
                        phi_check_Phi_L1T.append((muon_matched_Phi_L1T_phi[0]-muon_matched_Phi_L1T_phi[1]))                      
                    
                    
                    muon1_Phi = ROOT.TLorentzVector()
                    muon2_Phi = ROOT.TLorentzVector()

                    muon1_Phi.SetPtEtaPhiM(muon_matched_Phi_L1T_pt[0],
                                     muon_matched_Phi_L1T_eta[0],
                                     muon_matched_Phi_L1T_phi[0],0.10566)
                    muon2_Phi.SetPtEtaPhiM(muon_matched_Phi_L1T_pt[1],
                                     muon_matched_Phi_L1T_eta[1],
                                     muon_matched_Phi_L1T_phi[1],0.10566)
            
                    di_muon_inv_mass_Phi.append((muon1_Phi+muon2_Phi).M())
                    di_muon_dR_matched_L1T_Phi.append(((muon_matched_Phi_L1T_eta[0]-muon_matched_Phi_L1T_eta[1])**2+(phi_check_Phi_L1T[0])**2)**0.5)
                    inv_M_o_dR_di_muon_Phi.append(((muon1_Phi+muon2_Phi).M()/((muon_matched_Phi_L1T_eta[0]-muon_matched_Phi_L1T_eta[1])**2+(phi_check_Phi_L1T[0])**2)**0.5))

    
    
    # Here Im writing info into files...
            
    
    fon=open("Muon_dR_from_JPsi_Phi_genP_Y.txt","a")
    for b in range(len(di_muon_dR_JPsi_genP)):
        fon.write(str(di_muon_dR_JPsi_genP[b])+"\t"+str(di_muon_invM_JPsi_genP[b])+"\t"+str(di_muon_invM_o_dR_JPsi_genP[b])+"\t"+str(di_muon_dR_Phi_genP[b])+"\t"+str(di_muon_invM_Phi_genP[b])+"\t"+str(di_muon_invM_o_dR_Phi_genP[b])+"\t"+str(di_muon_dR_JPsi_genP_iev[b])+"\t"+str(aiziet)+"\n")
    fon.close()
    
    fob=open("Muon_dR_from_Matching_genP_L1T_JPsi_Y.txt","a")
    for h in range(len(muon_dR_min_JPsi)):
        fob.write(str(muon_dR_min_JPsi[h])+"\t"+str(aiziet)+"\n")
    fob.close()
    
    foc=open("Muon_dR_from_Matching_genP_L1T_Phi_Y.txt","a")
    for h in range(len(muon_dR_min_Phi)):
        foc.write(str(muon_dR_min_Phi[h])+"\t"+str(aiziet)+"\n")
    foc.close()
    
    fov=open("Di_Muon_dR_invM_invM_o_dR_JPsi_Y.txt","a")
    for d in range(len(di_muon_dR_matched_L1T_JPsi)):
        fov.write(str(di_muon_dR_matched_L1T_JPsi[d])+"\t"+str(di_muon_inv_mass_JPsi[d])+"\t"+str(inv_M_o_dR_di_muon_JPsi[d])+"\t"+str(aiziet)+"\n")
    fov.close()  
    
    fol=open("Di_Muon_dR_invM_invM_o_dR_Phi_Y.txt","a")
    for s in range(len(di_muon_dR_matched_L1T_Phi)):
        fol.write(str(di_muon_dR_matched_L1T_Phi[s])+"\t"+str(di_muon_inv_mass_Phi[s])+"\t"+str(inv_M_o_dR_di_muon_Phi[s])+"\t"+str(aiziet)+"\n")
    fol.close()              
                
    print("List of di-muon dR from J/Psi and Phi in genP is this long -> "+str(len(di_muon_dR_JPsi_genP))+" "+str(len(di_muon_dR_Phi_genP)))
    print("List of matched J/Psi and Phi genP->L1T muons from pairs is this long -> "+str(len(muon_dR_min_JPsi))+" "+str(len(muon_dR_min_Phi)))             
    print("List of di-muon dR from J/Psi and Phi in matched L1T muons is this long -> "+str(len(di_muon_dR_matched_L1T_JPsi))+" "+str(len(di_muon_dR_matched_L1T_Phi)))                           
    
print("We had {0} strange matches, 1 from J/Psi, 1 from Phi would fit better. We removed them".format(len(shady_match_list)))              
print("We had {0} strange errors, with same parameters for muon and anti-muon".format(len(error_list)))
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
                              
                                                                                                                   
