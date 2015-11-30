from HighMassVHbbAnalysis.Heppy.vhbbobj import ptRel
from PhysicsTools.HeppyCore.utils.deltar import deltaR,deltaPhi
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Jet


#from math import *
import ROOT
import array
import math


def deltaR2(a, b):
    phi1 = float(a.phi())
    eta1 = float(a.eta())

    phi2 = float(b.phi())
    eta2 = float(b.eta())

    return pow(abs(abs(phi1 - phi2) - math.pi) - math.pi, 2) + pow(eta1 - eta2, 2)




class JetRegressionAK08 :
    def __init__(self,weightfile,name) :
#      weights = { "../weights/Zll_weights_phys14.xml" , "../weights/Wln_weights_phys14.xml", "../weights/Znn_weights_phys14.xml"}		
#      reader_name = {"jet0Regression_zll", "jet0Regression_wln","jet0Regression_znn"}	
#      for i in range(0,3):
        reader = ROOT.TMVA.Reader()

       
        self.FatjetAK08ungroomed_pt = array.array('f',[0])
        self.FatjetAK08pruned_pt = array.array('f',[0])
        self.FatjetAK08prunedCal_pt = array.array('f',[0])
        self.FatjetAK08prunedCal_eta = array.array('f',[0])
        self.FatjetAK08ungroomed_vertexNTracks = array.array('f',[0])
        self.FatjetAK08ungroomed_SV_mass_0 =   array.array('f',[0])
        self.FatjetAK08ungroomed_SV_EnergyRatio_0 =   array.array('f',[0])
        self.FatjetAK08ungroomed_SV_EnergyRatio_1 =   array.array('f',[0])
        self.FatjetAK08ungroomed_PFLepton_ptrel =   array.array('f',[0])
        self.FatjetAK08ungroomed_nSL =   array.array('f',[0])        


        reader.AddVariable("FatjetAK08ungroomed_pt[0]",self.FatjetAK08ungroomed_pt)
        reader.AddVariable("FatjetAK08pruned_pt[0]",self.FatjetAK08pruned_pt)
        reader.AddVariable("FatjetAK08prunedCal_pt[0]",self.FatjetAK08prunedCal_pt)
        reader.AddVariable("FatjetAK08prunedCal_eta[0]",self.FatjetAK08prunedCal_eta)
        reader.AddVariable("FatjetAK08ungroomed_vertexNTracks[0]",self.FatjetAK08ungroomed_vertexNTracks)
        reader.AddVariable("FatjetAK08ungroomed_SV_mass_0[0]",self.FatjetAK08ungroomed_SV_mass_0)
        reader.AddVariable("FatjetAK08ungroomed_SV_EnergyRatio_0[0]",self.FatjetAK08ungroomed_SV_EnergyRatio_0)
        reader.AddVariable("FatjetAK08ungroomed_SV_EnergyRatio_1[0]",self.FatjetAK08ungroomed_SV_EnergyRatio_1)
        reader.AddVariable("FatjetAK08ungroomed_PFLepton_ptrel[0]",self.FatjetAK08ungroomed_PFLepton_ptrel)
        reader.AddVariable("FatjetAK08ungroomed_nSL[0]",self.FatjetAK08ungroomed_nSL)
        reader.BookMVA(name,weightfile)
        self.reader=reader
        self.name=name

    def evaluateRegressionAK08(self, event):
#self.readCollections( event.input )
        reg_fj = []
        for ung_fj in getattr(event, "ak08"):         

            # We need the closest ungroomed fatjet to get the JEC:                                                                                                                                              # - Make a list of pairs: deltaR(ungroomed fj, groomed fj) for all ungroomed fatjets                                                                                                                # - Sort by deltaR                                                                                                                                                                                  # - And take the minimum            
            if len(getattr(event, "ak08pruned")):
                closest_pr_fj_and_dr = sorted( [(pr_fj, deltaR2(ung_fj, pr_fj)) for pr_fj in getattr(event, "ak08pruned")],  key=lambda x:x[1])[0] 
            else:
                print "WARNING: No pruned fatjets found in event with ungroomed fatjet. Skipping"
                continue

            # Use the jet cone size for matching                                                                                                                                                    
            minimal_dr_groomed_ungroomed = 0.8
            if closest_pr_fj_and_dr[1] > minimal_dr_groomed_ungroomed:
                print "WARNING: No pruned fatjet found close to ungroomed fatjet. Skipping"
                continue
            
            pr_jet = Jet(closest_pr_fj_and_dr[0])
            
            if len(getattr(event, "ak08prunedcal")):
                closest_cal_fj_and_dr = sorted( [(cal_fj, deltaR2(ung_fj, cal_fj)) for cal_fj in getattr(event, "ak08prunedcal")],  key=lambda x:x[1])[0] 
            else:
                print "WARNING: No calib groomed fatjets found in event with ungroomed fatjet. Skipping"
                continue

            if closest_cal_fj_and_dr[1] > minimal_dr_groomed_ungroomed:
                print "WARNING: No calib fatjet found close to ungroomed fatjet. Skipping"
                continue


            cal_jet = Jet(closest_cal_fj_and_dr[0])

            
#        for j in event.FatjetAK08pruned:

            self.FatjetAK08ungroomed_pt[0] = ung_fj.pt()
            print 'ung_fj.pt() ', ung_fj.pt()
            self.FatjetAK08pruned_pt[0] = pr_jet.pt()
            print 'pr_jet.pt() ',  pr_jet.pt()
            self.FatjetAK08prunedCal_pt[0] = cal_jet.pt()
            print 'cal_jet.pt() ', cal_jet.pt()
            self.FatjetAK08prunedCal_eta[0] = cal_jet.eta()
            self.FatjetAK08ungroomed_vertexNTracks[0] = ung_fj.vertexNTracks
            self.FatjetAK08ungroomed_SV_mass_0[0] = ung_fj.SV_mass_0
            self.FatjetAK08ungroomed_SV_EnergyRatio_0[0] = ung_fj.SV_EnergyRatio_0
            self.FatjetAK08ungroomed_SV_EnergyRatio_1[0] = ung_fj.SV_EnergyRatio_1
            self.FatjetAK08ungroomed_PFLepton_ptrel[0] = ung_fj.PFLepton_ptrel
            self.FatjetAK08ungroomed_nSL[0] = ung_fj.nSL
           
             # Need to do a deep-copy. Otherwise the original jet will be modified                                                                                                                  
            reg_groomed_fj = PhysicsObject(cal_jet).__copy__()
            reg_groomed_fj.scaleEnergy(self.reader.EvaluateRegression(self.name)[0])


            reg_fj.append(reg_groomed_fj)
            print 'reg_groomed_fj.pt() ', reg_groomed_fj.pt()
  
        
        setattr(event, 'ak08prunedreg',  reg_fj )
           
                #j.pt_reg = self.reader.EvaluateRegression(self.name)[0]
		
            


