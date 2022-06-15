#!/bin/bash

emin=10
emax=600
etamin=1.7
etamax=2.7
pid=130
n_events=500

seed=`shuf -i1-99999999 -n1`

echo "Start preparing environment!"
cmsrel  CMSSW_12_4_0_pre2
cd CMSSW_12_4_0_pre2/src
cmsenv
git cms-init
git cms-merge-topic Abhirikshma:TICLHackathon
git checkout Abhirikshma/TICLHackathon
scram b -j 

echo "Clone Repository"
git clone ssh://git@gitlab.cern.ch:7999/phzehetn/ticl_hackathon.git

echo "Run step 1"
cmsRun -n 4 ticl_hackathon/production/step1.py filename=step1_$1_$2.root Nevents=$n_events PID=$pid Emin=$emin Emax=$emax Etamin=$etamin Etamax=$etamax seed=$seed

echo "Run step 2"
cmsRun -n 4 ticl_hackathon/production/step2.py sourcename=step1_$1_$2.root filename=step2_$1_$2.root Nevents=$n_events

echo "Run step 3"
cmsRun -n 4 ticl_hackathon/production/step3.py sourcename=step2_$1_$2.root filename=step3_$1_$2.root Nevents=$n_events

echo "Moving files"
mv step1_${1}_${2}.root /eos/cms/store/group/dpg_hgcal/comm_hgcal/hackathon/samples/close_by_single_kaon/production/step1_${1}_${2}.root
mv step2_${1}_${2}.root /eos/cms/store/group/dpg_hgcal/comm_hgcal/hackathon/samples/close_by_single_kaon/production/step2_${1}_${2}.root
mv step3_${1}_${2}.root /eos/cms/store/group/dpg_hgcal/comm_hgcal/hackathon/samples/close_by_single_kaon/production/step3_${1}_${2}.root
