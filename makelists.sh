#! /bin/bash

while read -r num; do
    [ -z "$num" ] && continue
    #ls /sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/TTrees/Data/Pro001_pcdb001_v001_VtxCut0_dst_calofitting-000"$num"-part*/* > "filelists/${num}.list"
    hadd /sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/TTrees/Data/Pro001_pcdb001_v001_VtxCut0_dst_calofitting-000"$num"-part*/*
done < runlist_oo.list
