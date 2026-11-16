#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/phenix/u/${LOGNAME}

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=~/install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


# final MC pass: totalcalo_et/vz reweighting + low-stat gating + truth pT1-pT2 prior reweighting
# (needs hists/prior_reweight.root from getPriorReweighting.C, run on the primer pass from runmcrw.sh)
inname="/sphenix/tg/tg01/jets/jpark4/Run25OO/TTrees/SkimmedTrees/outree_skimjet_MC_inclusive_merged.root"
outname="/sphenix/user/vbailey/ooanalysis/hists/histMC_reweight.root"
root -b -q /sphenix/user/vbailey/ooanalysis/getDijets.C\(\"$inname\",\"$outname\",1,1,5,1\)
