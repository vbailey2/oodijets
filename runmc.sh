#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/phenix/u/${LOGNAME}

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=~/install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


inname="/sphenix/tg/tg01/jets/jpark4/Run25OO/TTrees/SkimmedTrees/outree_skimjet_MC_inclusive_merged.root"
outname="/sphenix/user/vbailey/ooanalysis/hists/histMC.root"
root -b -q /sphenix/user/vbailey/ooanalysis/getDijets.C\(\"$inname\",\"$outname\",1\)
