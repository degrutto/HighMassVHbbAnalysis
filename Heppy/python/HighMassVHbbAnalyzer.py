from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.objects.JetAnalyzer import cleanJetsAndLeptons
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR,deltaPhi
from copy import deepcopy
from math import *
import itertools
import ROOT

def Boost(self,boost):
   bx=boost.X()
   by=boost.Y()
   bz=boost.Z()
   b2 = bx*bx + by*by + bz*bz; 
   gamma = 1.0 / sqrt(1.0 - b2);
   bp = bx*self.X() + by*self.Y() + bz*self.Z();
   gamma2 =  (gamma - 1.0)/b2 if b2 >0 else  0.0;
   #print gamma2,gamma,bp,self.X(),self.Y(),self.Z(),self.T()
   self.SetXYZT(self.X() + gamma2*bp*bx + gamma*bx*self.T(),
   self.Y() + gamma2*bp*by + gamma*by*self.T(),
   self.Z() + gamma2*bp*bz + gamma*bz*self.T(),
   (gamma*(self.T() + bp)))
   return self

class HighMassVHbbAnalyzer( Analyzer ):
    '''Analyze VH events
    '''
    # === Adding ID function for muons ===
    def isGoodMuon(self,lepton):#self,lepton):
      #if((lepton.isGlobalMuon()) and (lepton.globalTrack().hitPattern().numberOfValidMuonHits()>0) and (lepton.numberOfMatchedStations()>1) and ((lepton.muonBestTrack().ptError()/lepton.muonBestTrack().pt())<0.3) and (fabs(lepton.muonBestTrack().dxy()<0.2 and lepton.muonBestTrack().isNonnull())) and (fabs(lepton.muonBestTrack().dz())<0.5 and lepton.muonBestTrack().isNonnull()) and (lepton.innerTrack().hitPattern().numberOfValidPixelHits()>0) and (lepton.innerTrack().hitPattern().trackerLayersWithMeasurement()>5)):
      if(lepton.isGlobalMuon() and lepton.globalTrack().hitPattern().numberOfValidMuonHits() > 0 and lepton.numberOfMatchedStations() > 1 and lepton.muonBestTrack().ptError()/lepton.muonBestTrack().pt()<0.3 and lepton.dB()< 0.2 and lepton.innerTrack().hitPattern().numberOfValidPixelHits()>0 and lepton.innerTrack().hitPattern().trackerLayersWithMeasurement()>5):
        #print "We have a GOOD Muon!"
        return True
      else:
        #print "We have a BAD Muon"
        return False

    # === Adding ID function for electrons ===
    def isGoodElectron(self,lepton):#self,lepton):
      if ((abs(lepton.superCluster().eta())<1.4442 and lepton.ecalDriven() and abs(lepton.deltaEtaSeedClusterTrackAtVtx())<0.004 and abs(lepton.deltaPhiSuperClusterTrackAtVtx())< 0.06 and (lepton.hadronicOverEm()<1.0/lepton.superCluster().energy()+0.05) and (lepton.e2x5Max()/lepton.e5x5()>0.94 or lepton.e1x5()/lepton.e5x5()>0.83) and abs(lepton.dxy())<0.02) or (abs(lepton.superCluster().eta())>1.566 and abs(lepton.superCluster().eta())<2.5 and lepton.ecalDriven() and abs(lepton.deltaEtaSeedClusterTrackAtVtx())<0.006 and abs(lepton.deltaPhiSuperClusterTrackAtVtx())< 0.06 and (lepton.hadronicOverEm()<5.0/lepton.superCluster().energy()+0.05) and abs(lepton.dxy())<0.05 and lepton.full5x5_sigmaIetaIeta()<0.03)):
        return True
      #if((lepton.pt()>35) and (lepton.ecalDrivenSeed()==1)):#isEcalDriven()==1)):
      #  if(abs(lepton.eta())<1.4442):                             #Barrel
      #     if((abs(lepton.deltaEtaSuperClusterTrackAtVtx())<0.004) and (abs(lepton.deltaPhiSuperClusterTrackAtVtx())<0.06) and (lepton.hadronicOverEm()<(1/lepton.pt())+0.05) and (((lepton.e2x5Max()/lepton.e5x5())>0.94) or ((lepton.e1x5()/lepton.e5x5())>0.83)) and (lepton.dr03TkSumPt()<5) and (lepton.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS)<=1) and  (abs(lepton.dxy())<0.02)): #and ((lepton.dr03EcalRecHitSumEt()+lepton.dr03HcalDepth1TowerSumEt())<(2+0.03*(lepton.pt())+0.28*(lepton.rho())))):
      #       print "Barrel good ele"
      #       return True
      #  elif(abs(lepton.eta())>1.566 and abs(lepton.eta())<2.5): #EndCap
      #     if((abs(lepton.deltaEtaSuperClusterTrackAtVtx())<0.006) and (abs(lepton.deltaPhiSuperClusterTrackAtVtx())<0.06) and (lepton.hadronicOverEm()<(5/lepton.pt())+0.05) and (lepton.full5x5_sigmaIetaIeta()>0.03) and (lepton.dr03TkSumPt()<5) and (lepton.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS)<=1) and  (abs(lepton.dxy())<0.05)): #and ((lepton.dr03EcalRecHitSumEt()+lepton.dr03HcalDepth1TowerSumEt())<(2+0.03*(lepton.pt())+0.28*(lepton.rho())))):
      #       print "ENDCAP GOOD Ele!"
      #       return True
      else:
        print "We have a BAD Ele"
        return False
      return False

    # === Adding ID function for jets ===
    def isGoodLooseJet(self,jet):#self,lepton):
      if(((jet.neutralHadronEnergyFraction()<0.99) and (jet.neutralEmEnergyFraction()<0.99) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1)) and ((abs(jet.eta())<=2.4) and (jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.99) or (abs(jet.eta())>2.4)) and (abs(jet.eta())<=3.0)):
        return True #for |eta|<=3
      elif((jet.neutralEmEnergyFraction()<0.90) and (jet.neutralMultiplicity()>10) and (abs(jet.eta())>3.0)):
        return True #for |eta|>3
      elif((jet.neutralHadronEnergyFraction()<0.99) and (jet.neutralEmEnergyFraction()<0.99) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1)):
        if(abs(jet.eta())>=-3.0 and abs(jet.eta())<=3.0):
           return True
        elif(abs(jet.eta())>=-2.4 and abs(jet.eta())<=2.4):
          if((jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.99)):
            return True
          else:
            return False
        else:
          return False
      else:
        print "We have a BAD Jet"
        return False

    def isGoodTightJet(self,jet):#self,lepton):
      if(((jet.neutralHadronEnergyFraction()<0.90) and (jet.neutralEmEnergyFraction()<0.90) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1)) and ((abs(jet.eta())<=2.4) and (jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.99) or (abs(jet.eta())>2.4)) and (abs(jet.eta())<=3.0)):
        return True #for |eta|<=3
      elif((jet.neutralEmEnergyFraction()<0.90) and (jet.neutralMultiplicity()>10) and (abs(jet.eta())>3.0)):
          return True #for |eta|>3
      elif((jet.neutralHadronEnergyFraction()<0.90) and (jet.neutralEmEnergyFraction()<0.90) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1)):
        if(abs(jet.eta())>=-3.0 and abs(jet.eta())<=3.0):
           return True
        elif(abs(jet.eta())>=-2.4 and abs(jet.eta())<=2.4):
          if((jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.99)):
            return True
          else:
            return False
        else:
          return False
      else:
        print "We have a BAD Jet"
        return False

    def isGoodTightLepVetoJet(self,jet):#self,lepton):
       if(((jet.neutralHadronEnergyFraction()<0.90) and (jet.neutralEmEnergyFraction()<0.90) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1) and (jet.muonEnergyFraction    ()<0.8)) and ((abs(jet.eta())<=2.4) and (jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.90) or (abs(jet.eta())>2.4)) and (abs(jet.eta())<=3.0)):
         return True #for |eta|<=3
       elif((jet.neutralHadronEnergyFraction()<0.90) and (jet.neutralEmEnergyFraction()<0.90) and ((jet.chargedMultiplicity()+jet.neutralMultiplicity())>1) and (jet.muonEnergyFraction()<0.8)):
         if(abs(jet.eta())>=-3.0 and abs(jet.eta())<=3.0):
            return True
         elif(abs(jet.eta())>=-2.4 and abs(jet.eta())<=2.4):
           if((jet.chargedHadronEnergyFraction()>0) and (jet.chargedMultiplicity()>0) and (jet.chargedEmEnergyFraction()<0.90)):
             return True
           else:
             return False
         else:
           return False
       else:
         print "We have a BAD Jet"
         return False


    def declareHandles(self):
        super(HighMassVHbbAnalyzer, self).declareHandles()
        if self.cfg_comp.isMC:
            self.handles['GenInfo'] = AutoHandle( ('generator','',''), 'GenEventInfoProduct' )
   


    def beginLoop(self,setup):
        super(HighMassVHbbAnalyzer,self).beginLoop(setup)
        if "outputfile" in setup.services :
            setup.services["outputfile"].file.cd()
            self.inputCounter = ROOT.TH1F("Count","Count",1,0,2)
            self.inputCounterWeighted = ROOT.TH1F("CountWeighted","Count with gen weight and pu weight",1,0,2)
            self.inputCounterPosWeight = ROOT.TH1F("CountPosWeight","Count genWeight>0",1,0,2)
            self.inputCounterNegWeight = ROOT.TH1F("CountNegWeight","Count genWeight<0",1,0,2)
   


#    def doFakeMET(self,event):
#	#fake MET from Zmumu
#	event.fakeMET = ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.)
#	event.fakeMET.sumet = 0
#	if event.Vtype == 0 :
#		event.fakeMET=event.met.p4() + event.V
#                event.fakeMET.sumet = event.met.sumEt() - event.V.pt()

    def doHtMhtJets30(self,event):
        ## with Central Jets
        objects30 = [ j for j in event.cleanJets if j.pt() > 30 ] + event.selectedLeptons
        event.htJet30 = sum([x.pt() for x in objects30])
        event.mhtJet30vec = ROOT.reco.Particle.LorentzVector(-1.*(sum([x.px() for x in objects30])) , -1.*(sum([x.py() for x in objects30])), 0, 0 )             
        event.mhtJet30 = event.mhtJet30vec.pt()
        event.mhtPhiJet30 = event.mhtJet30vec.phi()







    def process(self, event):
	#print "Event number",event.iEv
        self.readCollections( event.input )
        self.inputCounter.Fill(1)
        if self.cfg_comp.isMC:
            genWeight = self.handles['GenInfo'].product().weight()
            self.inputCounterWeighted.Fill(1,copysign(1.0,genWeight)*event.puWeight)
            if genWeight > 0:
                self.inputCounterPosWeight.Fill(1)
            elif genWeight < 0:
                self.inputCounterNegWeight.Fill(1)

#        self.initOutputs(event)

#	event.pfCands = self.handles['pfCands'].product()
# 	met = event.met
        #self.addNewBTag(event)
        #Add CSV ranking                                                                                                                                                                             
        csvSortedJets=sorted(event.cleanJetsAll, key =  lambda jet : jet.btag(getattr(self.cfg_ana,"btagDiscriminator",'pfCombinedInclusiveSecondaryVertexV2BJetTags')),reverse=True)
        for j in event.cleanJetsAll:
              j.btagIdx=csvSortedJets.index(j)
        for j in event.discardedJets:
              j.btagIdx=-1

#	self.doFakeMET(event)
	self.doHtMhtJets30(event)
      
	#substructure threshold, make configurable
	ssTrheshold = 200.
	# filter events with less than 2 jets with pt 20

#        self.fillTauIndices(event)

#        event.jetsForHiggs = [x for x in event.cleanJets if self.cfg_ana.JetsPreSelection(x) ]

     

	#if not  len(event.jetsForHiggs) >= 2 : # and event.jetsForHiggs[1] > 20.) : # or(len(event.cleanJets) == 1 and event.cleanJets[0] > ssThreshold ) ) :
        #   return self.cfg_ana.passall
        #if event.Vtype < 0 and not ( sum(x.pt() > 30 for x in event.jetsForHiggs) >= 4 or sum(x.pt() for x in event.jetsForHiggs[:4]) > 160 ):
        #   return self.cfg_ana.passall
     

        if not len(getattr(event, "ak08pruned")) or (len(getattr(event, "ak08pruned"))>0 and getattr(event, "ak08pruned")[0].pt()<ssTrheshold  ) : 
            return self.cfg_ana.passall          


        ## Clean Jets from leptons
        leptons = []
        if hasattr(event, 'inclusiveLeptons'):
          for l in event.inclusiveLeptons: #added
            if(abs(l.pdgId())==13): #added
              print "MUONE" #added
              #leptons.isMyGoodMuon = self.isGoodMuon(l) #added
              #event.inclusiveLeptons.isMyGoodMuon = self.isGoodMuon(l) #added
              l.isMyGoodMuon = self.isGoodMuon(l) #added
              #self.isMyGoodMuon = self.isGoodMuon(l) #added
              print "Result", l.isMyGoodMuon #added
            elif(abs(l.pdgId())==11): #added
              if(abs(l.eta())<1.442):
                l.isBarrelEle = True
                l.isEndCapEle = False
              elif(abs(l.eta())>1.566 and abs(l.eta())<2.5):
                l.isBarrelEle = False
                l.isEndCapEle = True
              print "Elettrone!"      #added
              l.isMyGoodElectron = self.isGoodElectron(l)  #added
              print "Result", l.isMyGoodElectron #added
        #if hasattr(event, 'inclusiveLeptons'): #removed
          leptons = [ l for l in event.inclusiveLeptons if ( l.pt()>80 and  ( (abs(l.pdgId()) ==13 and l.relIso03<0.4 ) or (abs(l.pdgId())==11 and l.relIso04 <0.4) ) ) ]
              
        if hasattr(event, 'ak08'):
            AK08Jets = [ j for j in event.ak08  ]
        if hasattr(event, 'ak08pruned'):
            AK08prunedJets = [ j for j in event.ak08pruned  ]

        if hasattr(event, 'ak08prunedsubjets'):
            AK08prunedsubJets = [ j for j in event.ak08prunedsubjets  ]

        if hasattr(event, 'ak08prunedcal'):
            AK08prunedcalJets = [ j for j in event.ak08prunedcal  ]
        if hasattr(event, 'ak08prunedreg'):
            AK08prunedregJets = [ j for j in event.ak08prunedreg  ]

 
        event.cleanAK08JetsAll = []
        event.cleanAK08JetsAll, cleanLeptons = cleanJetsAndLeptons(AK08Jets, leptons, 0.8 , lambda jet,lepton : lepton)

        event.cleanAK08prunedJetsAll = []
        event.cleanAK08prunedJetsAll, cleanLeptons = cleanJetsAndLeptons(AK08prunedJets, leptons, 0.8 , lambda jet,lepton : lepton)

        event.cleanAK08prunedcalJetsAll = []
        event.cleanAK08prunedcalJetsAll, cleanLeptons = cleanJetsAndLeptons(AK08prunedcalJets, leptons, 0.8 , lambda jet,lepton : lepton)

        event.cleanAK08prunedregJetsAll = []
        event.cleanAK08prunedregJetsAll, cleanLeptons = cleanJetsAndLeptons(AK08prunedregJets, leptons, 0.8 , lambda jet,lepton : lepton)

        event.cleanAK08prunedregJetsAll = []
        event.cleanAK08prunedregJetsAll, cleanLeptons = cleanJetsAndLeptons(AK08prunedregJets, leptons, 0.8 , lambda jet,lepton : lepton)

        event.cleanAK08prunedsubJetsAll = []
        event.cleanAK08prunedsubJetsAll, cleanLeptons = cleanJetsAndLeptons(AK08prunedsubJets, leptons, 0.4 , lambda jet,lepton : lepton)



        return True

        
