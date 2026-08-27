#! /bin/bash

for iter in 1 2 3 4 5 6 7 8
do
	root -b -q unfoldData.C\($iter\)
	root -b -q projectDijets.C\(\"hists/hist-unfoldedData_iter$iter.root\",\"hists/projections_unfold_iter$iter.root\"\)
done
