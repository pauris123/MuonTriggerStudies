import time
bigBang=time.time()
import uproot
import ROOT
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

names=pn.read_csv("Step3_FileNames.txt", sep='\t',names=["nos"]) 
print("I have the libraries and names")

notik=0
tik=len(names["nos"])

for aiziet in range(notik,tik):
    print("\n {0}/{1}, {2}".format(aiziet+1,tik,wall_time(time.time()-bigBang)))

    start=time.time()
        
    rootfile="/eos/cms/store/group/phys_bphys/chic_hlt/"+names["nos"][aiziet]
    tree = uproot.open(rootfile)["Events"]
    print(" {0}\tI have the tree".format(l_wall_time(time.time()-start)))


    #subsection = ["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.ptUnconstrained_","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.etaAtVtx_",
    #             "l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.phiAtVtx_"]

    subsection = ["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta",
                 "l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"]


    ds = tree.pandas.df(subsection, entrystop=None)
    print(" {0}\tI have the branches".format(l_wall_time(time.time()-start)))

    muon_charge = []
    muon_pt = []
    muon_eta = []
    muon_phi = []
    
    for i in range(len(ds)):
        try:
              
            if len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 1 :
                
                muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][0])
                muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][0])
                muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][0])
                muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][0])                
            
            elif len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 2 :
                for y in range(2):
                    
                    muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][y])
                    muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][y])
                    muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][y])
                    muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][y]) 
    
            elif len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 3 :
                for y in range(3):
                    
                    muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][y])
                    muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][y])
                    muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][y])
                    muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][y])
    
            elif len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 4 :
                for y in range(4):
                    
                    muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][y])
                    muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][y])
                    muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][y])
                    muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][y])
    
            elif len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 5 :
                for y in range(5):
                    
                    muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][y])
                    muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][y])
                    muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][y])
                    muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][y])
                    
            elif len(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i]) == 6 :
                for y in range(6):
                    
                    muon_charge.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_"][i][y])
                    muon_pt.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"][i][y])
                    muon_eta.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta"][i][y])
                    muon_phi.append(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"][i][y])
            
            #print(str(len(muon_charge)))  
            
        except:
            continue
    
    print("List of muon charge is this long -> "+str(len(muon_charge)))     
    
    muon_mass = []
    muon_InvMoDR = []
    
    for k in range(len(muon_charge)-1):
    
        for j in range(k+1,len(muon_charge)):
            if muon_charge[k] + muon_charge[j] == 1:  
                if (((muon_eta[k]-muon_eta[j])**2+(muon_phi[k]-muon_phi[j])**2)**0.5<0.3 and ((muon_eta[k]-muon_eta[j])**2+(muon_phi[k]-muon_phi[j])**2)**0.5>0.02):
                    muon1 = ROOT.TLorentzVector()
                    muon2 = ROOT.TLorentzVector()

                    muon1.SetPtEtaPhiM(muon_pt[k],
                                     muon_eta[k],
                                     muon_phi[k],0.10566)
                    muon2.SetPtEtaPhiM(muon_pt[j],
                                     muon_eta[j],
                                     muon_phi[j],0.10566)
            
                    muon_mass.append((muon1+muon2).M())
                    muon_InvMoDR.append(((muon1+muon2).M()/((muon_eta[k]-muon_eta[j])**2+(muon_phi[k]-muon_phi[j])**2)**0.5))
        
                else:
                    continue
            else:
                continue
                
    fo=open("Muon_invM_DR_inFile.txt","a")
    for z in range(len(muon_mass)):
        fo.write(str(muon_mass[z])+"\t"+str(muon_InvMoDR[z])+"\n")
    fo.close()                
                
    print("List of muon dimass is this long -> "+str(len(muon_mass)))            
                
                
                
                
                
print("{0}\tAll done".format(wall_time(time.time()-bigBang)))
