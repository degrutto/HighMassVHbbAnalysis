#! /usr/bin/env python
fastObjects=True

#Switch to True to produce x1, x2, id1, id2, pdf scale
doPDFVars = False

import ROOT
from DataFormats.FWLite import *
import PhysicsTools.HeppyCore.framework.config as cfg
from HighMassVHbbAnalysis.Heppy.vhbbobj import *
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 


import logging
logging.basicConfig(level=logging.ERROR)

import os

cfg.Analyzer.nosubdir = True

treeProducer= cfg.Analyzer(
	class_object=AutoFillTreeProducer, 
	defaultFloatType = "F",
	verbose=False, 
	vectorTree = True,
        globalVariables	= [
                 NTupleVariable("json", lambda ev : getattr(ev,"json",True), help="Passing json selection"),
                 NTupleVariable("nPU0", lambda ev : [bx.nPU() for bx in  ev.pileUpInfo if bx.getBunchCrossing()==0][0], help="nPU in BX=0",mcOnly=True),
                 NTupleVariable("nPVs", lambda ev : len(ev.goodVertices), help="total number of good PVs"),
#		 NTupleVariable("Vtype", lambda ev : ev.Vtype, help="Event classification"),
#		 NTupleVariable("VtypeSim", lambda ev : ev.VtypeSim, help="Event classification",mcOnly=True),
#		 NTupleVariable("VMt", lambda ev : ev.V.goodMt, help="Transverse mass of the vector boson"),
#		 NTupleVariable("HVdPhi", lambda ev : deltaPhi(ev.V.phi(),ev.H.phi()), help="Delta phi between Higgs and Z/W"),
#		 NTupleVariable("fakeMET_sumet", lambda ev : ev.fakeMET.sumet, help="Fake SumET from Zmumu events removing muons"),
		 NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
#		 NTupleVariable("deltaR_jj",  lambda ev: deltaR(ev.hJets[0].eta(),ev.hJets[0].phi(),ev.hJets[1].eta(),ev.hJets[1].phi()) if len(ev.hJets) > 1 else -1, float, help="deltaR higgsJets"),
		 NTupleVariable("lheNj",  lambda ev: ev.lheNj, float,mcOnly=True, help="number of jets at LHE level"),
                 NTupleVariable("lheNb",  lambda ev: ev.lheNb, float,mcOnly=True, help="number of b-jets at LHE level"),
                 NTupleVariable("lheNc",  lambda ev: ev.lheNc, float,mcOnly=True, help="number of c-jets at LHE level"),
                 NTupleVariable("lheNg",  lambda ev: ev.lheNg, float,mcOnly=True, help="number of gluon jets at LHE level"),
                 NTupleVariable("lheNl",  lambda ev: ev.lheNl, float,mcOnly=True, help="number of light(uds) jets at LHE level"),
		 NTupleVariable("lheV_pt",  lambda ev: ev.lheV_pt, float,mcOnly=True, help="Vector pT at LHE level"),
                 NTupleVariable("lheHT",  lambda ev: ev.lheHT, float,mcOnly=True, help="HT at LHE level"),
#                 NTupleVariable("genTTHtoTauTauDecayMode", lambda ev: ev.genTTHtoTauTauDecayMode, int,mcOnly=True, help="gen level ttH, H -> tautau decay mode"),        
		#Soft Activity vars
##                 NTupleVariable("totSoftActivityJets2", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 2 ] ), int, help="number of jets from soft activity with pt>2Gev"),
# #                NTupleVariable("totSoftActivityJets5", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 5 ] ), int, help="number of jets from soft activity with pt>5Gev"),
#  #               NTupleVariable("totSoftActivityJets10", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 10 ] ), int, help="number of jets from soft activity with pt>10Gev"),
#                 NTupleVariable("ttCls",  lambda ev: getattr(ev, "ttbarCls", -1), float,mcOnly=True, help="ttbar classification via GeNHFHadronMatcher"),
		 NTupleVariable("mhtJet30",  lambda ev : ev.mhtJet30, help="mht with jets30"),
		 NTupleVariable("mhtPhiJet30",  lambda ev : ev.mhtPhiJet30, help="mht phi with jets30"),
		 NTupleVariable("htJet30",  lambda ev : ev.htJet30, help="ht  with jets30"),
		 NTupleVariable("met_rawpt",  lambda ev : ev.met.uncorPt(), help="raw met"),
		 NTupleVariable("metPuppi_pt",  lambda ev : ev.metPuppi.pt(), help="met from Puppi"),
		 NTupleVariable("metPuppi_phi",  lambda ev : ev.metPuppi.phi(), help="met phi from Puppi"),
		 NTupleVariable("metPuppi_rawpt",  lambda ev : ev.metPuppi.uncorPt(), help="raw met from Puppi"),
		 NTupleVariable("metNoHF_pt",  lambda ev : ev.metNoHF.pt(), help="met from NoHF"),
		 NTupleVariable("metNoHF_phi",  lambda ev : ev.metNoHF.phi(), help="met phi from NoHF"),
		 NTupleVariable("metNoHF_rawpt",  lambda ev : ev.metNoHF.uncorPt(), help="raw met from NoHF"),
		 NTupleVariable("metType1p2_pt",  lambda ev : ev.met.shiftedPt(12,2), help="type1type2met"),
		 NTupleVariable("metNoPU_pt",  lambda ev : ev.metNoPU.pt(), help="PFnoPU E_{T}^{miss}"),
		 NTupleVariable("metNoPU_phi",  lambda ev : ev.metNoPU.phi(), help="phi of PFnoPU E_{T}^{miss}"),
		 NTupleVariable("tkMet_pt",  lambda ev : ev.tkMet.pt(), help="E_{T}^{miss} from tracks"),
		 NTupleVariable("tkMet_phi",  lambda ev : ev.tkMet.phi(), help="phi of E_{T}^{miss} from tracks"),
		 NTupleVariable("tkMetPVchs_pt",  lambda ev : ev.tkMetPVchs.pt(), help="E_{T}^{miss} from tracks"),
		 NTupleVariable("tkMetPVchs_phi",  lambda ev : ev.tkMetPVchs.phi(), help="phi of E_{T}^{miss} from tracks"),
#		 NTupleVariable("isrJetVH",  lambda ev : ev.isrJetVH, help="Index of ISR jet in VH"),
		 NTupleVariable("Flag_hbheIsoFilter",  lambda ev : ev.hbheFilterIso, help="hbheFilterIso, after rerun"),
		 NTupleVariable("Flag_hbheFilterNew",  lambda ev : ev.hbheFilterNew, help="hbheFilterIso, after rerun"),
		 NTupleVariable("simPrimaryVertex_z", lambda ev: ev.genvertex, float,mcOnly=True, help="z coordinate of the simulated primary vertex"),
		 NTupleVariable("genHiggsDecayMode", lambda ev: ev.genHiggsDecayMode, float, mcOnly=True, help="decay mode of the Higgs boson"),
	],
	globalObjects = {
          "met"    : NTupleObject("met",     metType, help="PF E_{T}^{miss}, after default type 1 corrections"),
#          "fakeMET"    : NTupleObject("fakeMET", fourVectorType, help="fake MET in Zmumu event obtained removing the muons"),
#          "H"    : NTupleObject("H", fourVectorType, help="higgs"),
#          "HCSV"    : NTupleObject("HCSV", fourVectorType, help="higgs CSV selection"),
#          "H_reg"    : NTupleObject("H_reg", fourVectorType, help="regressed higgs"),
#          "HCSV_reg"    : NTupleObject("HCSV_reg", fourVectorType, help="regresses higgs CSV selection"),
#          "HaddJetsdR08"    : NTupleObject("HaddJetsdR08", fourVectorType, help="higgs with cen jets added if dR<0.8 from hJetsCSV selection"),
#          "V"    : NTupleObject("V", fourVectorType, help="z or w"),
#          "softActivityJets"    : NTupleObject("softActivity", softActivityType, help="VBF soft activity variables"),
#          "softActivityVHJets"    : NTupleObject("softActivityVH", softActivityType, help="VH soft activity variables"),
        },
	collections = {
		#standard dumping of objects
   	        "selectedLeptons" : NTupleCollection("selLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
 	        "inclusiveLeptons" : NTupleCollection("incLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
		#old style stuff, can be removed at some point
#   	        "vLeptons" : NTupleCollection("vLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
#   	        "aLeptons" : NTupleCollection("aLeptons", leptonTypeVHbb, 8, help="Additional leptons, not passing the preselection"),
## now store only indices, this lines are left commented for possible debugging
##	        "hJets"       : NTupleCollection("hJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Higgs jets"),
##	        "aJets"       : NTupleCollection("aJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Additional jets"),
#                "hjidx"       : NTupleCollection("hJidx",    objectInt, 2,help="Higgs jet indices"),
#                "hjidxDiJetPtByCSV"       : NTupleCollection("hJidx_sortcsv",    objectInt, 2,help="Higgs jet indices within hJets with CSV sorting "),
#                "ajidx"       : NTupleCollection("aJidx",    objectInt, 8,help="additional jet indices"),
#                "hjidxCSV"       : NTupleCollection("hJCidx",    objectInt, 2,help="Higgs jet indices CSV"),
#                "ajidxCSV"       : NTupleCollection("aJCidx",    objectInt, 8,help="additional jet indices CSV"),
                "cleanJetsAll"       : NTupleCollection("Jet",     jetTypeVHbb, 25, help="Cental+fwd jets after full selection and cleaning, sorted by b-tag"),
#treeProducer.collections["cleanAK08JetsAll"] = NTupleCollection("FatjetCleanAK08ungroomed",  ak8FatjetType,  10,                                                                                    
                "cleanAK08JetsAll" : NTupleCollection("FatjetCleanAK08ungroomed",     ak8FatjetType, 10, help="AK08 ungroomed jets, lepton cleaning applied"),
		"cleanAK08prunedJetsAll" : NTupleCollection("FatjetCleanAK08pruned",     fourVectorType, 10, help="AK08 pruned jets, lepton cleaning applied"),
		"cleanAK08prunedcalJetsAll" : NTupleCollection("FatjetCleanAK08prunedcal",     fourVectorType, 10, help="AK08 pruned cal jets, lepton cleaning applied"),
		"cleanAK08prunedregJetsAll" : NTupleCollection("FatjetCleanAK08prunedcalreg",    fourVectorType, 10, help="AK08 pruned cal jets, lepton cleaning applied"),
		"cleanAK08prunedsubJetsAll" : NTupleCollection("SubjetCleanAK08pruned",    subjetType, 10, help="AK08 pruned sub jets, lepton cleaning applied"),
             


#                "hjidxaddJetsdR08"       : NTupleCollection("hjidxaddJetsdR08",    objectInt, 5,help="Higgs jet indices with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),
#                "ajidxaddJetsdR08"       : NTupleCollection("ajidxaddJetsdR08",    objectInt, 8,help="additional jet indices with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),
#		"dRaddJetsdR08"       : NTupleCollection("dRaddJetsdR08",    objectFloat, 5,help="dR of add jet with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),        
#                "discardedJets"       : NTupleCollection("DiscardedJet",     jetTypeVHbb, 15, help="jets that were discarded"),
                "inclusiveTaus"  : NTupleCollection("TauGood", tauTypeVHbb, 25, help="Taus after the preselection"),
#                "softActivityJets"    : NTupleCollection("softActivityJets", fourVectorType, 5, help="jets made for soft activity"),
#                "softActivityVHJets"    : NTupleCollection("softActivityVHJets", fourVectorType, 5, help="jets made for soft activity VH version"),
                "goodVertices"    : NTupleCollection("primaryVertices", primaryVertexType, 4, help="first four PVs"),

		#dump of gen objects
                #"generatorSummary"    : NTupleCollection("GenSummary", genParticleWithLinksType, 30, help="Generator summary, see description in Heppy GeneratorAnalyzer",mcOnly=True),
                "genJets"    : NTupleCollection("GenJet",  genJetType, 15, help="Generated jets with hadron matching, sorted by pt descending",filter=lambda x: x.pt() > 20,mcOnly=True),
                "ak8genJets"    : NTupleCollection("ak08GenJet",   genParticleType, 15, help="Generated jets (without hadron matching, but with neutrinos), sorted by pt descending",filter=lambda x: x.pt() > 20,mcOnly=True), 
                "genHiggsSisters"    : NTupleCollection("GenHiggsSisters",     genParticleType, 4, help="Sisters of the Higgs bosons"),
                "gentopquarks"    : NTupleCollection("GenTop",     genParticleType, 4, help="Generated top quarks from hard scattering"),
                "genallstatus2bhadrons"    : NTupleCollection("GenStatus2bHad",     genParticleType, 15, help="Generated Status 2 b Hadrons"),
#                "genallstatus2bhadronsv2"    : NTupleCollection("GenStatus2bHad_v2",     genParticleType, 15, help="Generated Status 2 b Hadrons"),
                "gennus"          : NTupleCollection("GenNu", genParticleType, 4, help="Generated neutrinos"), #added
                "genMyTaus"       : NTupleCollection("GenTau", genParticleType, 4, help="Generated taus"), #added
                "gennusFromTop"    : NTupleCollection("GenNuFromTop",     genParticleType, 4, help="Generated neutrino from t->W decay"),
                "genbquarksFromH"      : NTupleCollection("GenBQuarkFromH",  genParticleType, 4, help="Generated bottom quarks from Higgs decays"),
                "genbquarksFromTop"      : NTupleCollection("GenBQuarkFromTop",  genParticleType, 4, help="Generated bottom quarks from top decays"),
                "genbquarksFromHafterISR"      : NTupleCollection("GenBQuarkFromHafterISR",  genParticleType, 4, help="Generated bottom quarks from Higgs decays"),
                "genwzquarks"     : NTupleCollection("GenWZQuark",   genParticleType, 6, help="Generated quarks from W/Z decays"),
                "genlepsFromWPrime": NTupleCollection("GenLepFromWPrime", genParticleType, 4, help="Generated leptons from W'"), 
                "genleps"         : NTupleCollection("GenLep",     genParticleType, 4, help="Generated leptons from W/Z decays"),
                "gentauleps"         : NTupleCollection("GenTaus",     genParticleType, 4, help="Generated taus"),
                "genlepsFromTop"         : NTupleCollection("GenLepFromTop",     genParticleType, 4, help="Generated leptons from t->W decays"),
                "gentauleps"      : NTupleCollection("GenLepFromTau", genParticleType, 6, help="Generated leptons from decays of taus from W/Z decays"),
		"genHiggsBosons"   : NTupleCollection("GenHiggsBoson", genParticleType, 4, help="Generated Higgs boson "),
		#"genZbosonsToLL"  : NTupleCollection("GenZbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		#"genWbosonsToLL"  : NTupleCollection("GenWbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		"genvbosons"       : NTupleCollection("GenVbosons", genParticleType, 6, help="Generated W or Z bosons, mass > 30"),
		"pileUpVertex_z"       : NTupleCollection("pileUpVertex_z",    objectFloat, 5,help="z position of hardest pile-up collisions"),        
		"pileUpVertex_ptHat"       : NTupleCollection("pileUpVertex_ptHat",    objectFloat, 5,help="z position of hardest pile-up collisions"),        
	}
	)

#Create shifted MET Ntuples
metNames={y.__get__(ROOT.pat.MET):x for x,y in ROOT.pat.MET.__dict__.items() if  (x[-2:]=="Up" or x[-4:]=="Down")}
print "met Names", metNames
shifted_met_keys = ["met_shifted_{0}".format(n) for n in range(12)] #we do not need noShift I gueess
shifted_met_names = ["met_shifted_%s"%metNames[n] for n in range(12)] #we do not need noShift I gueess
shifted_mets = {mk: NTupleObject(nm, shiftedMetType, help="PF E_{T}^{miss}, after default type 1 corrections, shifted with %s" %mk) for mk,nm in zip(shifted_met_keys,shifted_met_names)}
treeProducer.globalObjects.update(shifted_mets)

btag_weights = {}
for syst in ["JES", "LF", "HF", "Stats1", "Stats2", "cErr1", "cErr2"]:
	for sdir in ["Up", "Down"]:
		name = "bTagWeight"+syst+sdir
		btag_weights[name] = NTupleVariable("bTagWeight_" + syst + sdir,
			lambda ev, sname=syst+sdir: bweightcalc.calcEventWeight(
				ev.cleanJetsAll, kind="final", systematic=sname
			), float, mcOnly=True, help="b-tag CSV weight, variating "+syst+" "+sdir
		)
btag_weights["bTagWeight"] = NTupleVariable("bTagWeight",
	lambda ev: bweightcalc.calcEventWeight(
		ev.cleanJetsAll, kind="final", systematic="nominal"
	), float ,mcOnly=True, help="b-tag CSV weight, nominal"
)
#print list(btag_weights.values())
treeProducer.globalVariables += list(btag_weights.values())

# Lepton Analyzer, take its default config and fix loose iso consistent with tight definition
  #removed
from PhysicsTools.Heppy.analyzers.objects.LeptonAnalyzer import LeptonAnalyzer
LepAna = LeptonAnalyzer.defaultConfig
LepAna.mu_isoCorr  = "deltaBeta"
LepAna.loose_muon_isoCut = lambda muon : muon.relIso04 < 0.4 
LepAna.ele_isoCorr = "rhoArea"
LepAna.loose_electron_isoCut = lambda electron : electron.relIso03 < 0.4 

#added
#from PhysicsTools.Heppy.analyzers.objects.LeptonAnalyzer import LeptonAnalyzer
#LepAna = cfg.Analyzer(
#      LeptonAnalyzer, name="leptonAnalyzer",
#      # input collections
#      muons='slimmedMuons',
#      electrons='slimmedElectrons',
#      rhoMuon= 'fixedGridRhoFastjetAll',
#      rhoElectron = 'fixedGridRhoFastjetAll',
#      # energy scale corrections and ghost muon suppression (off by default)
#      doMuonScaleCorrections=False,
#      doElectronScaleCorrections=False, # "embedded" in 5.18 for regression
#      doSegmentBasedMuonCleaning=False,
#      # inclusive very loose muon selection
#      inclusive_muon_id  = "POG_ID_Loose",
#      inclusive_muon_pt  = 3,
#      inclusive_muon_eta = 2.4,
#      inclusive_muon_dxy = 1000,
#      inclusive_muon_dz  = 1000,
#      muon_dxydz_track = "innerTrack",
#      # veto muon selection
#      loose_muon_id     = "POG_ID_Loose",
#      loose_muon_pt     = 10,
#      loose_muon_eta    = 2.4,
#      loose_muon_dxy    = 1000,
#      loose_muon_dz     = 1000,
#      loose_muon_isoCut = (lambda mu : ( mu.relIso04 <= 0.4)), # this is not to apply loose_muon_relIso which is on DR=0.3
#      loose_muon_relIso = 0.4,
#      # inclusive very loose electron selection
#      inclusive_electron_id  = "",
#      inclusive_electron_pt  = 5,
#      inclusive_electron_eta = 2.5,
#      inclusive_electron_dxy = 0.5,
#      inclusive_electron_dz  = 1.0,
#      inclusive_electron_lostHits = 5.0,
#      # veto electron selection
#      loose_electron_id     = "POG_Cuts_ID_SPRING15_25ns_v1_ConvVetoDxyDz_Veto_full5x5",
#      loose_electron_pt     = 10,
#      loose_electron_eta    = 2.5,
#      loose_electron_dxy    = 0.5,
#      loose_electron_dz     = 1.0,
#      loose_electron_relIso = 1.0,
#      loose_electron_lostHits = 5.0,
#      # muon isolation correction method (can be "rhoArea" or "deltaBeta")
#      mu_isoCorr = "deltaBeta" ,
#      mu_effectiveAreas = "Spring15_25ns_v1", #(can be 'Data2012' or 'Phys14_25ns_v1' or 'Spring15_25ns_v1')
#      # electron isolation correction method (can be "rhoArea" or "deltaBeta")
#      ele_isoCorr = "rhoArea" ,
#      ele_effectiveAreas = "Spring15_25ns_v1" , #(can be 'Data2012' or 'Phys14_25ns_v1' or 'Spring15_25ns_v1' or 'Spring15_50ns_v1')
#      ele_tightId = "Cuts_2012" ,
#      # Mini-isolation, with pT dependent cone: will fill in the miniRelIso, miniRelIsoCharged, miniRelIsoNeutral variables of the leptons (see https://indico.cern.ch/event/368826/ )
#      doMiniIsolation = False, # off by default since it requires access to all PFCandidates 
#      packedCandidates = 'packedPFCandidates',
#      miniIsolationPUCorr = 'rhoArea', # Allowed options: 'rhoArea' (EAs for 03 cone scaled by R^2), 'deltaBeta', 'raw' (uncorrected), 'weights' (delta beta weights; not validated)
#      miniIsolationVetoLeptons = None, # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
#      # minimum deltaR between a loose electron and a loose muon (on overlaps, discard the electron)
#      min_dr_electron_muon = 0.05,
#      # do MC matching 
#      do_mc_match = True, # note: it will in any case try it only on MC, not on data
#      match_inclusiveLeptons = False, # match to all inclusive leptons
#      do_mc_match_photons = False, # do not do MC matching of electrons to photons
#)



from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer import VertexAnalyzer
VertexAna = VertexAnalyzer.defaultConfig
VertexAna.keepFailingEvents = True

from PhysicsTools.Heppy.analyzers.objects.PhotonAnalyzer import PhotonAnalyzer
PhoAna = PhotonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.TauAnalyzer import TauAnalyzer
TauAna = TauAnalyzer.defaultConfig
TauAna.inclusive_ptMin = 18.
TauAna.inclusive_etaMax = 2.5
TauAna.inclusive_dxyMax = 1000.
TauAna.inclusive_dzMax = 0.4
TauAna.inclusive_vetoLeptons = False
TauAna.inclusive_leptonVetoDR = 0.4
TauAna.inclusive_decayModeID = "decayModeFindingNewDMs"
TauAna.inclusive_tauID = "decayModeFindingNewDMs"
TauAna.inclusive_vetoLeptonsPOG = False
TauAna.inclusive_tauAntiMuonID = ""
TauAna.inclusive_tauAntiElectronID = ""
TauAna.inclusive_tauLooseID = "decayModeFindingNewDMs"

from PhysicsTools.Heppy.analyzers.objects.JetAnalyzer import JetAnalyzer
JetAna = JetAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.gen.LHEAnalyzer import LHEAnalyzer 
LHEAna = LHEAnalyzer.defaultConfig

#from PhysicsTools.Heppy.analyzers.gen.GeneratorAnalyzer import GeneratorAnalyzer 
from HighMassVHbbAnalysis.Heppy.GeneratorAnalyzer import GeneratorAnalyzer
GenAna = GeneratorAnalyzer.defaultConfig
from HighMassVHbbAnalysis.Heppy.VHGeneratorAnalyzer import GeneratorAnalyzer as  VHGeneratorAnalyzer
VHGenAna = VHGeneratorAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.METAnalyzer import METAnalyzer
METAna = METAnalyzer.defaultConfig
METAna.doTkMet = True
METAna.doMetNoPU = True
METAna.doTkGenMet = False
METAna.includeTkMetPVLoose = False
METAna.includeTkMetPVTight = False

METPuppiAna = copy.copy(METAna)
METPuppiAna.metCollection     = "slimmedMETsPuppi"
METPuppiAna.doMetNoPU = False
METPuppiAna.recalibrate = False
METPuppiAna.collectionPostFix = "Puppi"
METPuppiAna.copyMETsByValue = True
METPuppiAna.doTkMet = False
METPuppiAna.doMetNoPU = False


METNoHFAna = copy.copy(METPuppiAna)
METNoHFAna.metCollection     = "slimmedMETsNoHF"
METNoHFAna.doMetNoPU = False
METNoHFAna.recalibrate = False
METNoHFAna.collectionPostFix = "NoHF"
METNoHFAna.copyMETsByValue = True
METNoHFAna.doTkMet = False
METNoHFAna.doMetNoPU = False




from PhysicsTools.Heppy.analyzers.core.PileUpAnalyzer import PileUpAnalyzer
PUAna = PileUpAnalyzer.defaultConfig

from HighMassVHbbAnalysis.Heppy.VHbbAnalyzer import VHbbAnalyzer

from HighMassVHbbAnalysis.Heppy.HighMassVHbbAnalyzer import HighMassVHbbAnalyzer

JetAna.jetPt = 15
JetAna.jetEta = 4.7
JetAna.doQG=False
JetAna.QGpath=os.environ['CMSSW_BASE']+"/src/PhysicsTools/Heppy/data/pdfQG_AK4chs_13TeV_v2b.root"
JetAna.recalibrateJets=True
JetAna.jecPath=os.environ['CMSSW_BASE']+"/src/HighMassVHbbAnalysis/Heppy/data/jec"
JetAna.mcGT="74X_mcRun2_asymptotic_v4"
JetAna.dataGT = "Summer15_25nsV6_DATA"
JetAna.addJECShifts=True
JetAna.addJERShifts=True

#mu_pfRelIso04 = lambda mu : (mu.pfIsolationR04().sumChargedHadronPt + max( mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumPUPt,0.0)) / mu.pt()

VHbb = cfg.Analyzer(
    verbose=False,
    class_object=VHbbAnalyzer,
#    wEleSelection = lambda x : x.pt() > 30 and x.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
#    wMuSelection = lambda x : x.pt() > 30 and x.muonID("POG_ID_Tight"), #and mu_pfRelIso04(x) < 0.15,
    wEleSelection = lambda x : x.pt() > 30 and x.electronID(""),
    wMuSelection = lambda x : x.pt() > 30 and x.muonID("POG_ID_HighPt"), #and mu_pfRelIso04(x) < 0.15,
    zEleSelection = lambda x : x.pt() > 15 and x.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    zMuSelection = lambda x : x.pt() > 10 and x.muonID("POG_ID_Loose"), #and mu_pfRelIso04(x) < 0.25,
    zLeadingElePt = 20,
    zLeadingMuPt = 20,
    higgsJetsPreSelection = lambda x: (( x.puJetId() > 0 and x.jetID('POG_PFID_Loose')) or abs(x.eta())>3.0 ) and x.pt() >  20 ,
    boostedhiggsJetsPreSelection = lambda x:  x.pt() >  200 ,
    higgsJetsPreSelectionVBF = lambda x: x.pt() >  20 ,
#    higgsJetsPreSelectionVBF = lambda x: (( x.puJetId() > 0 and x.jetID('POG_PFID_Loose')) or abs(x.eta())>3.0 ) and x.pt() >  20,
    passall=False,
    doSoftActivityVH=False,
    doVBF=False,
    regressions = [
        {"weight":"Zll-spring15.weights.xml", "name":"jet0Regression_zll", "vtypes":[0,1]},
        {"weight":"Wln-spring15.weights.xml", "name":"jet0Regression_wln", "vtypes":[2,3]},
        {"weight":"Znn-spring15.weights.xml", "name":"jet0Regression_znn", "vtypes":[4,5,-1]}
    ],
    regressionVBF = [ 		
	{"weight":"VBF-spring15.weights.xml", "name":"jet0Regression_vbf", "vtypes":[0,1,2,3,4,5,-1]}
    ],
    VBFblikelihood = {"weight":"TMVA_blikelihood_vbf_singlebtag_v13_id.xml", "name":"BDGT"}
)




HighMassVHbb = cfg.Analyzer(
    verbose=False,
    class_object=HighMassVHbbAnalyzer,
#    higgsJetsPreSelectionVBF = lambda x: (( x.puJetId() > 0 and x.jetID('POG_PFID_Loose')) or abs(x.eta())>3.0 ) and x.pt() >  20,
    passall=False,

)




from HighMassVHbbAnalysis.Heppy.TTHtoTauTauAnalyzer import TTHtoTauTauAnalyzer
TTHtoTauTau = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauAnalyzer,
)
from HighMassVHbbAnalysis.Heppy.TTHtoTauTauGeneratorAnalyzer import TTHtoTauTauGeneratorAnalyzer
TTHtoTauTauGen = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauGeneratorAnalyzer,
)

#from HighMassHighMassHighMassVHbbAnalysis.Heppy.HeppyShell import HeppyShell
#sh = cfg.Analyzer( class_object=HeppyShell)

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer

from HighMassVHbbAnalysis.Heppy.TriggerTable import triggerTable
from HighMassVHbbAnalysis.Heppy.TriggerTableData import triggerTable as triggerTableData

TrigAna = cfg.Analyzer(
    verbose = False,
    class_object = TriggerBitAnalyzer,
    triggerBits = triggerTable,  #default is MC, use the triggerTableData in -data.py files
#   processName = 'HLT',
#   outprefix = 'HLT'
   )

from PhysicsTools.HeppyCore.framework.services.tfile import TFileService 
output_service = cfg.Service(
      TFileService,
      'outputfile',
      name="outputfile",
      fname='tree.root',
      option='recreate'
    )

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer
FlagsAna = TriggerBitAnalyzer.defaultEventFlagsConfig

from HighMassVHbbAnalysis.Heppy.hbheAnalyzer import *
hbheAna = hbheAnalyzer.defaultConfig



from PhysicsTools.Heppy.analyzers.gen.PDFWeightsAnalyzer import PDFWeightsAnalyzer
PdfAna = cfg.Analyzer(PDFWeightsAnalyzer,
    PDFWeights = [],
    doPDFVars = doPDFVars
)

if doPDFVars:
    treeProducer.globalVariables += [
        NTupleVariable("pdf_x1",  lambda ev: ev.pdf_x1, float,mcOnly=True, help="PDF energy fraction of first parton"),
        NTupleVariable("pdf_x2",  lambda ev: ev.pdf_x2, float,mcOnly=True, help="PDF energy fraction of second parton"),
        NTupleVariable("pdf_id1",  lambda ev: ev.pdf_id1, int,mcOnly=True, help="PDF id of first parton"),
        NTupleVariable("pdf_id2",  lambda ev: ev.pdf_id2, int,mcOnly=True, help="PDF id of second parton"),
        NTupleVariable("pdf_scale",  lambda ev: ev.pdf_scale, float,mcOnly=True, help="PDF scale"),
    ]

TrigAna.unrollbits=True

from PhysicsTools.Heppy.analyzers.core.JSONAnalyzer import JSONAnalyzer
jsonAna = cfg.Analyzer(JSONAnalyzer,
      passAll=True
      )
#jsonAna.json= "json.txt"
#sequence = [jsonAna,LHEAna,FlagsAna, hbheAna, GenAna,VHGenAna,PUAna,TrigAna,VertexAna,LepAna,PhoAna,TauAna,JetAna,METAna, METPuppiAna, METNoHFAna, PdfAna, VHbb,TTHtoTauTau,TTHtoTauTauGen,treeProducer]#,sh]

sequence = [jsonAna,LHEAna,FlagsAna, hbheAna, GenAna,VHGenAna,PUAna,TrigAna,VertexAna,LepAna,PhoAna,TauAna,JetAna,METAna, METPuppiAna, METNoHFAna, PdfAna, HighMassVHbb , treeProducer]#,sh]


from PhysicsTools.Heppy.utils.miniAodFiles import miniAodFiles
sample = cfg.MCComponent(
    files = [
    #"root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/WprimeToWhToWlephbb_narrow_M-4000_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/E4EB03F2-4170-E511-B5C4-00259073E466.root"
    "root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/WprimeToMuNu_M-3000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/10EFFE5D-4374-E511-905E-003048F0113E.root"
    #"root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/WprimeToWhToWlephbb_narrow_M-1200_13TeV-madgraph/MINIAODSIM/Asympt50ns_74X_mcRun2_asymptotic50ns_v0-v1/30000/6248D38C-8D76-E511-B243-20CF305B058C.root"
#'file:68790C43-4971-E511-AFD7-782BCB20E307.root'
#		"root://cms-xrd-global.cern.ch//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/009AE141-CA6D-E511-A060-002590A3716C.root"

#'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/WprimeToWhToWlephbb_narrow_M-1200_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/DE7A5FB2-DC71-E511-9539-001517F8083C.root'
#'file:DE7A5FB2-DC71-E511-9539-001517F8083C.root'
#"root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/ZprimeToZhToZinvhbb_narrow_M-1200_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/8C77EFF8-4771-E511-97BD-B083FED40671.root"
#"/store/mc/RunIISpring15MiniAODv2/ZprimeToZhToZinvhbb_narrow_M-1200_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/8C77EFF8-4771-E511-97BD-B083FED40671.root"
#		"/scratch/arizzi/00B6C8DE-E76E-E511-AEDE-008CFA000BB8.root" ##ttbar
],

    #files = ["226BB247-A565-E411-91CF-00266CFF0AF4.root"],
    name="ZHLL125", isEmbed=False,
    puFileMC="puMC.root",
    puFileData="puData.root", 
    splitFactor = 5
    )
sample.isMC=True

# the following is declared in case this cfg is used in input to the heppy.py script
selectedComponents = [sample]
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence, 
		     services = [output_service],
                     events_class = Events)

class TestFilter(logging.Filter):
    def filter(self, record):
        print record

# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 1, nEvents = 1)

    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
