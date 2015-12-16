from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HBB_HEPPY_D14_003'
config.General.workArea = 'crab_projects_D14_003'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
                             'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                              '../vhbb_combined.py',
                             'TMVAClassification_BDT.weights.xml',
                             'puData.root',
                             'puMC.root',
                              'json.txt',
                              "../Zll-spring15.weights.xml",
                              "../Wln-spring15.weights.xml",
                              "../Znn-spring15.weights.xml",
                              "../VBF-spring15.weights.xml",
                             '../TMVA_blikelihood_vbf_singlebtag_v13_id.xml',
                              '../regAK08.weights.xml' 
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/ZH_HToBB_ZToLL_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 100
config.Data.outLFNDirBase = '/store/group/cmst3/user/degrutto/WprimeHBBHeppyD14_TEST/'
config.Data.publication = True
config.Data.outputDatasetTag = 'HBB_HEPPY_D14'

config.section_("Site")
config.Site.storageSite = "T2_IT_Rome"
#config.Site.whitelist = ["T2_CH_CERN"]


#config.Data.ignoreLocality = True
