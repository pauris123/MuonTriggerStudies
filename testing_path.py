import time
import uproot
import ROOT
import matplotlib.pyplot as plt
import numpy as np
import pandas as pn


print("I have the libraries and names")

rootfile="/eos/cms/store/group/phys_bphys/chic_hlt/crab_chic_PUFlat55To75_2021_20210827_172116/210827_152158/0000/step3_inAODSIM_10.root"
tree = uproot.open(rootfile)["Events"]

print("Got the root file")


subsection_L1T = ["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fEta",
                 "l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.hwCharge_","l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPhi"]

subsection_gen = ["recoGenParticles_genParticles__HLT.obj.m_state.p4Polar_.fCoordinates.fPt","recoGenParticles_genParticles__HLT.obj.m_state.p4Polar_.fCoordinates.fEta",
                  "recoGenParticles_genParticles__HLT.obj.m_state.pdgId_","recoGenParticles_genParticles__HLT.obj.m_state.p4Polar_.fCoordinates.fPhi"]
               
#subsection_gen = ["fPt","fEta","pdgId_","fPhi"]
                 
ds = tree.pandas.df(subsection_L1T, entrystop=None)
dg = tree.pandas.df(subsection_gen, entrystop=None)

print(ds["l1tMuonBXVector_gmtStage2Digis_Muon_RECO.obj.data_.m_state.p4Polar_.fCoordinates.fPt"])
print(dg["recoGenParticles_genParticles__HLT.obj.m_state.p4Polar_.fCoordinates.fPt"][:100])

print("Test done")
