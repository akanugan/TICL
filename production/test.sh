#!/bin/bash

emin=10
emax=600
etamin=1.7
etamax=2.7
pid=22
n_events=2

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
cmsRun -n 4 ticl_hackathon/production/step1.py filename=step1_$1_$2.root Nevents=10 PID=$pid Emin=$emin Emax=$emax Etamin=$etamin Etamax=$etamax
echo "Saving step 1 to step1_${1}_${2}.root"
ls step1_${1}_${2}.root

echo "Run step 2"
cmsRun -n 4 ticl_hackathon/production/step2.py sourcename=step1_$1_$2.root filename=step2_$1_$2.root Nevents=$n_events
echo "Saving step 2 to step2_${1}_${2}.root"
ls step2_${1}_${2}.root

echo "Run step 3"
cmsRun -n 4 ticl_hackathon/production/step3.py sourcename=step2_$1_$2.root filename=step3_$1_$2.root Nevents=$n_events
        echo "Running step 4"
        cmsRun -n 8 $codedir/step4.py sourcename=step3_${i}.root filename=step4_${i}.root Nevents=$n_events > log/log4_${i}.out 2> log/log4_${i}.err
echo "Saving step 3 to step3_${1}_${2}.root"
mv step3_${1}_${2}.root /eos/cms/store/group/dpg_hgcal/comm_hgcal/hackathon/samples/condor_test/step3_${1}_${2}.root
