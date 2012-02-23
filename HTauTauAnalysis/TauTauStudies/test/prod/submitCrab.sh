#!/bin/sh

multicrab -create -cfg multicrab_skim_ElecTau_17Feb2012_MC.cfg

crab -submit 1-500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim
crab -submit 501-1000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 1001-1500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 1501-2000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 2001-2500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 2501-3000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 3001-3500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 3501-4000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 4001-4500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 4501-5000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 5001-5500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 5501-6000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 6001-6500 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 6501-7000 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
crab -submit 7001-7157 -c QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v3_skim 
