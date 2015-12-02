# HighMassVHbbAnalysis

Checkout:

	cmsrel CMSSW_7_4_15
  cd CMSSW_7_4_15/src/
  cmsenv
  git cms-merge-topic vhbb:vhbbHeppy74X
  git cms-merge-topic gkasieczka:dev-ak8calib  --> needed only for the meantime
  git clone git@github.com:degrutto/HighMassVHbbAnalysis.git 
  scram b -j 20 
  cd  HighMassVHbbAnalysis/test
  python vhbb_combined.py

To run with crab:

	cd crab 
	source launchall.sh datasetWprime.txt 
