#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/phenix/u/${LOGNAME}

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=~/install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#NEW_PREFIX="/sphenix/user/vbailey/ooanalysis/analysis/JS-Jet/OODijets/offline/hists/hist"
#base="$1"

#if [[ "$base" =~ (-[0-9]+-part[0-9]+\.root)$ ]]; then
#    echo "${NEW_PREFIX}${BASH_REMATCH[1]}"
#    outname="${NEW_PREFIX}${BASH_REMATCH[1]}"
#    root -b -q /sphenix/user/vbailey/ooanalysis/analysis/JS-Jet/OODijets/offline/getDijets.C\(\"$base\",\"$outname\"\)
#else
#    echo "Error: '$base' doesn't match expected pattern" >&2
#    exit 1
#fi

inname="/sphenix/tg/tg01/jets/jpark4/Run25OO/TTrees/SkimmedTrees/outree_skimjet_data_trigger_jet12_merged.root"
outname="/sphenix/user/vbailey/ooanalysis/analysis/JS-Jet/OODijets/offline/hists/hist-full.root"
root -b -q /sphenix/user/vbailey/ooanalysis/analysis/JS-Jet/OODijets/offline/getDijets.C\(\"$inname\",\"$outname\",0\)
